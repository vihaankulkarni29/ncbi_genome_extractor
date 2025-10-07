#!/usr/bin/env python3
"""
NCBI Genome Extractor - A robust tool for downloading genome sequences from NCBI
"""

import argparse
import csv
import json
import logging
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from typing import List, Optional, Dict, Any
from urllib.parse import parse_qs, urlparse

import requests
from Bio import SeqIO
from tqdm import tqdm

from config import (
    NCBI_EMAIL, NCBI_API_KEY, DEFAULT_MAX_RESULTS, DEFAULT_RETRIES,
    DEFAULT_DELAY, DEFAULT_OUTPUT_DIR, NCBI_SEARCH_URL, NCBI_FETCH_URL, NCBI_BASE_URL
)


class NCBIExtractor:
    """Main class for NCBI genome extraction"""

    def __init__(self, email: str = NCBI_EMAIL, api_key: str = NCBI_API_KEY):
        self.email = email
        self.api_key = api_key
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': f'NCBIExtractor/1.0 (mailto:{self.email})'
        })

    def search_genomes(self, query: str, max_results: int = DEFAULT_MAX_RESULTS,
                      select_best: bool = True) -> List[str]:
        """Search NCBI for genome IDs matching the query"""
        if select_best:
            # Get more genomes initially to select the best ones
            initial_results = max_results * 5  # Get 5x more to have good selection
            initial_results = min(initial_results, 500)  # Cap at 500 to avoid API limits
        else:
            initial_results = max_results

        params = {
            'db': 'nuccore',
            'term': query,
            'retmax': initial_results,
            'retmode': 'xml',
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        # Remove None values
        params = {k: v for k, v in params.items() if v is not None}

        response = self._make_request(NCBI_SEARCH_URL, params)
        if not response:
            return []

        # Parse XML response
        root = ET.fromstring(response.text)
        ids = []
        for id_element in root.findall('.//Id'):
            ids.append(id_element.text)

        if select_best and len(ids) > max_results:
            # Select the best genomes based on metadata quality
            ids = self._select_best_genomes(ids, max_results)

        return ids

    def download_fasta(self, genome_id: str, output_dir: str = DEFAULT_OUTPUT_DIR,
                      retries: int = DEFAULT_RETRIES, accession: Optional[str] = None) -> tuple[bool, str]:
        """Download FASTA sequence for a single genome ID"""
        os.makedirs(output_dir, exist_ok=True)

        # Get accession number if not provided
        if not accession:
            accession = self._get_accession_from_id(genome_id)
            if not accession:
                accession = genome_id  # Fallback to genome_id if accession not found

        params = {
            'db': 'nuccore',
            'id': genome_id,
            'rettype': 'fasta',
            'retmode': 'text',
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        params = {k: v for k, v in params.items() if v is not None}

        for attempt in range(retries + 1):
            try:
                response = self._make_request(NCBI_FETCH_URL, params, stream=True)
                if not response:
                    continue

                filename = f"{accession}.fasta"
                filepath = os.path.join(output_dir, filename)

                # Download with progress bar for large files
                total_size = int(response.headers.get('content-length', 0))
                with open(filepath, 'wb') as f, tqdm(
                    desc=f"Downloading {accession}",
                    total=total_size,
                    unit='B',
                    unit_scale=True,
                    unit_divisor=1024,
                ) as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))

                logging.info(f"Successfully downloaded {accession}")
                return True, accession

            except Exception as e:
                logging.warning(f"Attempt {attempt + 1} failed for {accession}: {e}")
                if attempt < retries:
                    time.sleep(DEFAULT_DELAY * (2 ** attempt))  # Exponential backoff

        logging.error(f"Failed to download {accession} after {retries + 1} attempts")
        return False, accession

    def _get_accession_from_id(self, genome_id: str) -> Optional[str]:
        """Get accession number from genome ID"""
        try:
            summary_url = NCBI_BASE_URL + "esummary.fcgi"
            params = {
                'db': 'nuccore',
                'id': genome_id,
                'retmode': 'xml',
                'email': self.email,
                'api_key': self.api_key if self.api_key else None
            }

            params = {k: v for k, v in params.items() if v is not None}
            response = self._make_request(summary_url, params)

            if response:
                root = ET.fromstring(response.text)
                docsum = root.find('.//DocSum')
                if docsum is not None:
                    for item in docsum.findall('Item'):
                        if item.get('Name') == 'Caption':
                            return item.text

        except Exception as e:
            logging.warning(f"Failed to get accession for {genome_id}: {e}")

        return None

    def download_multiple_genomes(self, genome_ids: List[str], output_dir: str = DEFAULT_OUTPUT_DIR,
                                max_concurrent: int = 3) -> List[str]:
        """Download multiple genomes with progress tracking"""
        successful = []
        failed = []

        # Get accession numbers for all genomes first
        logging.info("Retrieving accession numbers for all genomes...")
        genome_accessions = {}
        for genome_id in genome_ids:
            accession = self._get_accession_from_id(genome_id)
            genome_accessions[genome_id] = accession or genome_id

        with tqdm(total=len(genome_ids), desc="Overall Progress") as pbar:
            for genome_id in genome_ids:
                success, accession = self.download_fasta(genome_id, output_dir, accession=genome_accessions[genome_id])
                if success:
                    successful.append(accession)
                else:
                    failed.append(accession)
                pbar.update(1)

        logging.info(f"Downloaded {len(successful)} genomes successfully")
        if failed:
            logging.warning(f"Failed to download {len(failed)} genomes: {failed}")

        return successful

    def parse_search_url(self, url: str) -> Optional[str]:
        """Extract search query from NCBI URL"""
        parsed = urlparse(url)
        if 'ncbi.nlm.nih.gov' not in parsed.netloc:
            return None

        query_params = parse_qs(parsed.query)
        if 'term' in query_params:
            return ' '.join(query_params['term'])

        # Try to extract from path
        path_match = re.search(r'/([^/]+)$', parsed.path)
        if path_match:
            return path_match.group(1).replace('+', ' ')

        return None

    def extract_metadata(self, genome_ids: List[str]) -> List[Dict[str, Any]]:
        """Extract comprehensive metadata for genome IDs"""
        metadata_list = []

        # Process in batches to avoid API limits
        batch_size = 100
        for i in range(0, len(genome_ids), batch_size):
            batch_ids = genome_ids[i:i + batch_size]
            batch_metadata = self._extract_metadata_batch(batch_ids)
            metadata_list.extend(batch_metadata)

        return metadata_list

    def _extract_metadata_batch(self, genome_ids: List[str]) -> List[Dict[str, Any]]:
        """Extract metadata for a batch of genome IDs"""
        metadata_list = []

        # Get document summaries from NCBI
        summary_url = NCBI_BASE_URL + "esummary.fcgi"
        params = {
            'db': 'nuccore',
            'id': ','.join(genome_ids),
            'retmode': 'xml',
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        params = {k: v for k, v in params.items() if v is not None}
        response = self._make_request(summary_url, params)

        if not response:
            logging.warning(f"Failed to get metadata for batch: {genome_ids[:3]}...")
            return [self._create_empty_metadata(gid) for gid in genome_ids]

        try:
            root = ET.fromstring(response.text)
            for docsum in root.findall('.//DocSum'):
                metadata = self._parse_docsum_metadata(docsum)
                metadata_list.append(metadata)
        except ET.ParseError as e:
            logging.error(f"Failed to parse metadata XML: {e}")
            return [self._create_empty_metadata(gid) for gid in genome_ids]

        return metadata_list

    def _select_best_genomes(self, genome_ids: List[str], max_results: int) -> List[str]:
        """Select the best genomes based on metadata quality"""
        logging.info(f"Evaluating {len(genome_ids)} genomes for metadata quality...")

        # Extract metadata for all genomes
        metadata_list = self.extract_metadata(genome_ids)

        # Score each genome based on metadata completeness
        scored_genomes = []
        for metadata in metadata_list:
            score = self._calculate_metadata_score(metadata)
            scored_genomes.append({
                'genome_id': metadata['genome_id'],
                'score': score,
                'metadata': metadata
            })

        # Sort by score (highest first) and select top genomes
        scored_genomes.sort(key=lambda x: x['score'], reverse=True)
        selected_genomes = scored_genomes[:max_results]

        selected_ids = [genome['genome_id'] for genome in selected_genomes]

        # Log selection results
        high_quality = sum(1 for genome in selected_genomes if genome['score'] >= 7)
        medium_quality = sum(1 for genome in selected_genomes if 4 <= genome['score'] < 7)
        low_quality = sum(1 for genome in selected_genomes if genome['score'] < 4)

        logging.info(f"Selected {len(selected_ids)} best genomes:")
        logging.info(f"  - High quality (score ≥7): {high_quality}")
        logging.info(f"  - Medium quality (score 4-6): {medium_quality}")
        logging.info(f"  - Low quality (score <4): {low_quality}")

        return selected_ids

    def _calculate_metadata_score(self, metadata: Dict[str, Any]) -> int:
        """Calculate a quality score for genome metadata (0-10 scale)"""
        score = 0

        # MIC data (most important for AMR research) - 3 points
        if metadata.get('mic_data') and len(metadata['mic_data']) > 0:
            mic_count = len(metadata['mic_data'])
            score += min(mic_count * 2, 3)  # Max 3 points for MIC data

        # BioSample data - 2 points
        if metadata.get('biosample'):
            score += 2

        # BioProject data - 1 point
        if metadata.get('bioproject'):
            score += 1

        # Collection date - 1 point
        if metadata.get('collection_date'):
            score += 1

        # Geographic data (country/isolation source) - 1 point
        if metadata.get('country') or metadata.get('isolation_source'):
            score += 1

        # Host information - 1 point
        if metadata.get('host'):
            score += 1

        # Antibiotic resistance data - 1 point
        if metadata.get('antibiotic_resistance') and len(metadata['antibiotic_resistance']) > 0:
            score += 1

        return min(score, 10)  # Cap at 10

    def _parse_docsum_metadata(self, docsum: ET.Element) -> Dict[str, Any]:
        """Parse metadata from NCBI DocSum XML"""
        metadata = {
            'genome_id': None,
            'accession': None,
            'title': None,
            'organism': None,
            'genus': None,
            'species': None,
            'biosample': None,
            'bioproject': None,
            'collection_date': None,
            'country': None,
            'host': None,
            'isolation_source': None,
            'mic_data': [],
            'resistance_phenotype': [],
            'antibiotic_resistance': []
        }

        # Extract basic information
        id_elem = docsum.find('Id')
        if id_elem is not None:
            metadata['genome_id'] = id_elem.text

        # Parse items
        for item in docsum.findall('Item'):
            name = item.get('Name')
            if name == 'Caption':
                metadata['accession'] = item.text
            elif name == 'Title':
                metadata['title'] = item.text
                # Extract organism from title
                if item.text:
                    self._extract_organism_from_title(metadata, item.text)
            elif name == 'Extra':
                # Extract BioSample and BioProject IDs
                if item.text:
                    self._extract_ids_from_extra(metadata, item.text)
            elif name == 'CreateDate':
                metadata['create_date'] = item.text
            elif name == 'UpdateDate':
                metadata['update_date'] = item.text

        # Try to get additional metadata from BioSample if available
        if metadata['biosample']:
            biosample_metadata = self._get_biosample_metadata(metadata['biosample'])
            if biosample_metadata:
                # Only update if we got actual data
                for key, value in biosample_metadata.items():
                    if value and (key not in metadata or not metadata[key]):
                        metadata[key] = value

        # Enhanced metadata validation and enrichment
        metadata = self._enhance_metadata_quality(metadata)

        return metadata

    def _extract_organism_from_title(self, metadata: Dict[str, Any], title: str):
        """Extract organism name from sequence title"""
        # Common patterns for organism names in titles
        patterns = [
            r'([A-Z][a-z]+ [a-z]+(?: subsp\.? [a-z]+)?(?: str\.? .+?)?)',
            r'([A-Z][a-z]+ [a-z]+)',
            r'([A-Z][a-z]+ sp\.)'
        ]

        for pattern in patterns:
            match = re.search(pattern, title, re.IGNORECASE)
            if match:
                organism = match.group(1).strip()
                metadata['organism'] = organism

                # Split into genus and species
                parts = organism.split()
                if len(parts) >= 2:
                    metadata['genus'] = parts[0]
                    metadata['species'] = ' '.join(parts[1:])
                break

    def _extract_ids_from_extra(self, metadata: Dict[str, Any], extra: str):
        """Extract BioSample and BioProject IDs from extra field"""
        # Look for BioSample
        biosample_match = re.search(r'SAM[NED]?\d+', extra, re.IGNORECASE)
        if biosample_match:
            metadata['biosample'] = biosample_match.group(0)

        # Look for BioProject
        bioproject_match = re.search(r'PRJ[NED][A-Z]\d+', extra, re.IGNORECASE)
        if bioproject_match:
            metadata['bioproject'] = bioproject_match.group(0)

    def _get_biosample_metadata(self, biosample_id: str) -> Dict[str, Any]:
        """Get additional metadata from BioSample database"""
        metadata = {}

        try:
            # Query BioSample database
            biosample_url = NCBI_BASE_URL + "efetch.fcgi"
            params = {
                'db': 'biosample',
                'id': biosample_id,
                'retmode': 'xml',
                'email': self.email,
                'api_key': self.api_key if self.api_key else None
            }

            params = {k: v for k, v in params.items() if v is not None}
            response = self._make_request(biosample_url, params)

            if response:
                # Parse BioSample XML for additional metadata
                metadata.update(self._parse_biosample_xml(response.text))

        except Exception as e:
            logging.warning(f"Failed to get BioSample metadata for {biosample_id}: {e}")

        return metadata

    def _parse_biosample_xml(self, xml_content: str) -> Dict[str, Any]:
        """Parse BioSample XML for AMR-related metadata"""
        metadata = {}

        try:
            root = ET.fromstring(xml_content)

            # Extract collection date (try multiple attribute names)
            for attr_name in ['collection_date', 'collection date', 'date']:
                collection_date = root.find(f'.//Attribute[@attribute_name="{attr_name}"]')
                if collection_date is not None and collection_date.text:
                    metadata['collection_date'] = collection_date.text
                    break

            # Extract country/geographic location (try multiple attribute names)
            for attr_name in ['geo_loc_name', 'geographic location', 'country', 'location']:
                country = root.find(f'.//Attribute[@attribute_name="{attr_name}"]')
                if country is not None and country.text:
                    metadata['country'] = country.text
                    break

            # Extract host (try multiple attribute names)
            for attr_name in ['host', 'host organism', 'source host']:
                host = root.find(f'.//Attribute[@attribute_name="{attr_name}"]')
                if host is not None and host.text:
                    metadata['host'] = host.text
                    break

            # Extract isolation source (try multiple attribute names)
            for attr_name in ['isolation_source', 'isolation source', 'source', 'sample type']:
                isolation = root.find(f'.//Attribute[@attribute_name="{attr_name}"]')
                if isolation is not None and isolation.text:
                    metadata['isolation_source'] = isolation.text
                    break

            # Extract additional metadata from structured attributes
            for attr in root.findall('.//Attribute'):
                attr_name = attr.get('attribute_name', '').lower()
                attr_value = attr.text or ''

                if not attr_value.strip():
                    continue

                # MIC data patterns (expanded)
                if any(keyword in attr_name for keyword in ['mic', 'minimum inhibitory concentration', 'mic_', 'mic50', 'mic90']):
                    antibiotic_name = attr_name.replace('mic', '').replace('minimum inhibitory concentration', '').replace('_', ' ').strip()
                    if not antibiotic_name:
                        antibiotic_name = 'unknown'
                    mic_data_entry = {
                        'antibiotic': antibiotic_name,
                        'value': attr_value,
                        'unit': self._extract_mic_unit(attr_value)
                    }
                    if 'mic_data' not in metadata:
                        metadata['mic_data'] = []
                    metadata['mic_data'].append(mic_data_entry)

                # Resistance phenotype (expanded)
                elif any(keyword in attr_name for keyword in ['resistance', 'phenotype', 'susceptibility', 'antibiotic resistance']):
                    if 'resistance_phenotype' not in metadata:
                        metadata['resistance_phenotype'] = []
                    metadata['resistance_phenotype'].append(attr_value)

                # Specific antibiotic resistance (expanded list)
                elif any(abx in attr_name for abx in [
                    'ampicillin', 'tetracycline', 'ciprofloxacin', 'gentamicin', 'trimethoprim',
                    'sulfamethoxazole', 'cefotaxime', 'ceftriaxone', 'meropenem', 'amikacin',
                    'azithromycin', 'clindamycin', 'erythromycin', 'vancomycin', 'linezolid',
                    'daptomycin', 'colistin', 'polymyxin', 'carbapenem', 'esbl', 'mrsa'
                ]):
                    resistance_entry = {
                        'antibiotic': attr_name,
                        'resistance': attr_value
                    }
                    if 'antibiotic_resistance' not in metadata:
                        metadata['antibiotic_resistance'] = []
                    metadata['antibiotic_resistance'].append(resistance_entry)

                # Additional useful metadata
                elif attr_name in ['strain', 'serotype', 'serovar', 'pathotype', 'virulence']:
                    metadata[attr_name] = attr_value

        except ET.ParseError as e:
            logging.warning(f"Failed to parse BioSample XML: {e}")
        except Exception as e:
            logging.warning(f"Error parsing BioSample metadata: {e}")

        return metadata

    def _extract_mic_unit(self, mic_value: str) -> str:
        """Extract MIC unit from MIC value string"""
        # Common MIC units
        units = ['ug/ml', 'μg/ml', 'mg/l', 'mcg/ml', 'µg/ml']
        for unit in units:
            if unit in mic_value.lower():
                return unit
        return 'unknown'

    def _enhance_metadata_quality(self, metadata: Dict[str, Any]) -> Dict[str, Any]:
        """Enhance and validate metadata quality"""
        # Try to extract additional information from title if missing
        if metadata.get('title') and not metadata.get('organism'):
            temp_metadata = metadata.copy()
            self._extract_organism_from_title(temp_metadata, metadata['title'])
            for key in ['organism', 'genus', 'species']:
                if temp_metadata.get(key) and not metadata.get(key):
                    metadata[key] = temp_metadata[key]

        # Validate and clean collection date
        if metadata.get('collection_date'):
            metadata['collection_date'] = self._normalize_date(metadata['collection_date'])

        # Ensure lists are properly initialized
        for list_field in ['mic_data', 'resistance_phenotype', 'antibiotic_resistance']:
            if list_field not in metadata:
                metadata[list_field] = []
            elif not isinstance(metadata[list_field], list):
                metadata[list_field] = [metadata[list_field]]

        # Add metadata quality score
        metadata['quality_score'] = self._calculate_metadata_score(metadata)

        return metadata

    def _normalize_date(self, date_str: str) -> str:
        """Normalize date strings to consistent format"""
        if not date_str:
            return date_str

        # Try to parse various date formats
        try:
            # Handle YYYY-MM-DD, YYYY/MM/DD, etc.
            if '-' in date_str or '/' in date_str:
                return date_str.strip()
            # Handle just year
            elif len(date_str) == 4 and date_str.isdigit():
                return f"{date_str}-01-01"
            else:
                return date_str.strip()
        except:
            return date_str.strip()

    def _create_empty_metadata(self, genome_id: str) -> Dict[str, Any]:
        """Create empty metadata structure for failed extractions"""
        return {
            'genome_id': genome_id,
            'accession': None,
            'title': None,
            'organism': None,
            'genus': None,
            'species': None,
            'biosample': None,
            'bioproject': None,
            'collection_date': None,
            'country': None,
            'host': None,
            'isolation_source': None,
            'mic_data': [],
            'resistance_phenotype': [],
            'antibiotic_resistance': [],
            'quality_score': 0
        }

    def save_metadata_to_json(self, metadata_list: List[Dict[str, Any]], output_file: str):
        """Save metadata to JSON format"""
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(metadata_list, f, indent=2, ensure_ascii=False)

    def save_metadata_to_csv(self, metadata_list: List[Dict[str, Any]], output_file: str):
        """Save metadata to CSV format"""
        if not metadata_list:
            return

        # Flatten nested structures for CSV
        flattened_data = []
        for item in metadata_list:
            flat_item = {}
            for key, value in item.items():
                if isinstance(value, list):
                    if key in ['mic_data', 'resistance_phenotype', 'antibiotic_resistance']:
                        # Convert lists to semicolon-separated strings
                        if value and isinstance(value[0], dict):
                            # Handle list of dictionaries (MIC data, antibiotic resistance)
                            flat_item[key] = '; '.join([f"{v.get('antibiotic', '')}: {v.get('value', '')}" for v in value])
                        else:
                            # Handle list of strings
                            flat_item[key] = '; '.join(value)
                    else:
                        flat_item[key] = str(value)
                else:
                    flat_item[key] = value
            flattened_data.append(flat_item)

        fieldnames = flattened_data[0].keys() if flattened_data else []
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(flattened_data)

    def _make_request(self, url: str, params: dict, stream: bool = False) -> Optional[requests.Response]:
        """Make HTTP request with error handling"""
        try:
            response = self.session.get(url, params=params, stream=stream, timeout=30)
            response.raise_for_status()
            return response
        except requests.RequestException as e:
            logging.error(f"Request failed: {e}")
            return None


def main():
    """Main CLI function"""
    parser = argparse.ArgumentParser(description="Download genome sequences and metadata from NCBI")
    parser.add_argument('--query', help='NCBI search query')
    parser.add_argument('--url', help='NCBI search URL')
    parser.add_argument('--output_dir', default=DEFAULT_OUTPUT_DIR,
                       help='Output directory for FASTA files')
    parser.add_argument('--max_results', type=int, default=DEFAULT_MAX_RESULTS,
                       help='Maximum number of genomes to download')
    parser.add_argument('--retries', type=int, default=DEFAULT_RETRIES,
                       help='Number of retry attempts')
    parser.add_argument('--delay', type=float, default=DEFAULT_DELAY,
                       help='Delay between requests')
    parser.add_argument('--log_level', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level')
    parser.add_argument('--metadata_format', choices=['json', 'csv'], default='json',
                       help='Format for metadata output (default: json)')
    parser.add_argument('--no_metadata', action='store_true',
                       help='Skip metadata extraction (not recommended for AMR research)')
    parser.add_argument('--no_quality_selection', action='store_true',
                       help='Disable quality-based genome selection (use first N results)')

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Validate arguments
    if not args.query and not args.url:
        parser.error("Either --query or --url must be provided")

    # Initialize extractor
    extractor = NCBIExtractor()

    # Get search query
    if args.url:
        query = extractor.parse_search_url(args.url)
        if not query:
            logging.error("Could not parse query from URL")
            sys.exit(1)
    else:
        query = args.query

    logging.info(f"Searching for genomes with query: {query}")

    # Search for genomes
    select_best = not args.no_quality_selection
    genome_ids = extractor.search_genomes(query, args.max_results, select_best)
    if not genome_ids:
        logging.warning("No genomes found for the query")
        sys.exit(0)

    logging.info(f"Found {len(genome_ids)} genomes")

    # Download genomes
    successful = extractor.download_multiple_genomes(
        genome_ids, args.output_dir
    )

    # Extract metadata by default (unless --no_metadata is specified)
    if successful and not args.no_metadata:
        logging.info("Extracting metadata for downloaded genomes...")

        # Extract metadata for successful downloads
        metadata_list = extractor.extract_metadata(successful)

        # Determine metadata output file
        base_name = os.path.basename(args.output_dir.rstrip(os.sep))
        metadata_file = f"{base_name}_metadata.{args.metadata_format}"
        metadata_path = os.path.join(args.output_dir, metadata_file)

        # Save metadata in requested format
        if args.metadata_format == 'json':
            extractor.save_metadata_to_json(metadata_list, metadata_path)
            logging.info(f"Metadata saved to {metadata_path}")
        elif args.metadata_format == 'csv':
            extractor.save_metadata_to_csv(metadata_list, metadata_path)
            logging.info(f"Metadata saved to {metadata_path}")

        # Log summary of extracted metadata
        mic_count = sum(len(m.get('mic_data', [])) for m in metadata_list)
        biosample_count = sum(1 for m in metadata_list if m.get('biosample'))
        resistance_count = sum(len(m.get('antibiotic_resistance', [])) for m in metadata_list)

        logging.info(f"Metadata extraction summary:")
        logging.info(f"  - Genomes with BioSample IDs: {biosample_count}")
        logging.info(f"  - MIC data entries: {mic_count}")
        logging.info(f"  - Antibiotic resistance entries: {resistance_count}")

    logging.info(f"Process completed. Downloaded {len(successful)} genomes to {args.output_dir}")


if __name__ == '__main__':
    main()