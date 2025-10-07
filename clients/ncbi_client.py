#!/usr/bin/env python3
"""
NCBI Client - Handles all interactions with NCBI Entrez API
"""

import logging
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from typing import List, Optional, Dict, Any

import requests
from Bio import SeqIO
from tqdm import tqdm

from config import (
    NCBI_EMAIL, NCBI_API_KEY, DEFAULT_MAX_RESULTS, DEFAULT_RETRIES,
    DEFAULT_DELAY, NCBI_SEARCH_URL, NCBI_FETCH_URL, NCBI_BASE_URL
)


class NCBIClient:
    """Client for NCBI Entrez API interactions"""

    def __init__(self, email: str = NCBI_EMAIL, api_key: str = NCBI_API_KEY):
        self.email = email
        self.api_key = api_key
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': f'FederatedGenomeHarvester/1.0 (mailto:{self.email})'
        })

    def fetch_genomes(self, query: str, max_results: int = DEFAULT_MAX_RESULTS) -> List[Dict[str, Any]]:
        """
        Search NCBI for genomes and return raw records

        Args:
            query: NCBI search query
            max_results: Maximum number of results to return

        Returns:
            List of raw genome records from NCBI
        """
        logging.info(f"Searching NCBI with query: {query}")

        # Search for genome IDs
        genome_ids = self._search_genomes(query, max_results)
        if not genome_ids:
            logging.warning("No genomes found for the query")
            return []

        logging.info(f"Found {len(genome_ids)} genomes")

        # Extract metadata for all genomes
        logging.info("Extracting metadata for genomes...")
        raw_records = self._extract_metadata_batch(genome_ids)

        return raw_records

    def download_fasta(self, accession: str, output_dir: str,
                      retries: int = DEFAULT_RETRIES) -> bool:
        """Download FASTA sequence for a genome accession"""
        os.makedirs(output_dir, exist_ok=True)

        params = {
            'db': 'nuccore',
            'id': accession,
            'rettype': 'fasta',
            'retmode': 'text',
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        params = {k: v for k, v in params.items() if v is not None}

        filename = f"{accession}.fasta"
        filepath = os.path.join(output_dir, filename)

        for attempt in range(retries + 1):
            try:
                response = self._make_request(NCBI_FETCH_URL, params, stream=True)
                if not response:
                    continue

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
                return True

            except Exception as e:
                logging.warning(f"Attempt {attempt + 1} failed for {accession}: {e}")
                if attempt < retries:
                    time.sleep(DEFAULT_DELAY * (2 ** attempt))  # Exponential backoff

        logging.error(f"Failed to download {accession} after {retries + 1} attempts")
        return False

    def _search_genomes(self, query: str, max_results: int) -> List[str]:
        """Search NCBI for genome IDs matching the query"""
        params = {
            'db': 'nuccore',
            'term': query,
            'retmax': max_results,
            'retmode': 'xml',
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        params = {k: v for k, v in params.items() if v is not None}

        response = self._make_request(NCBI_SEARCH_URL, params)
        if not response:
            return []

        # Parse XML response
        root = ET.fromstring(response.text)
        ids = []
        for id_element in root.findall('.//Id'):
            ids.append(id_element.text)

        return ids

    def _extract_metadata_batch(self, genome_ids: List[str]) -> List[Dict[str, Any]]:
        """Extract metadata for a batch of genome IDs"""
        metadata_list = []

        # Process in batches to avoid API limits
        batch_size = 100
        for i in range(0, len(genome_ids), batch_size):
            batch_ids = genome_ids[i:i + batch_size]
            batch_metadata = self._get_metadata_batch(batch_ids)
            metadata_list.extend(batch_metadata)

        return metadata_list

    def _get_metadata_batch(self, genome_ids: List[str]) -> List[Dict[str, Any]]:
        """Get metadata for a batch of genome IDs"""
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

        # Add quality score
        metadata['quality_score'] = self._calculate_metadata_score(metadata)

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

            # Extract MIC data and resistance information (expanded)
            mic_data = []
            resistance_phenotype = []
            antibiotic_resistance = []

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
                    mic_data.append(mic_data_entry)

                # Resistance phenotype (expanded)
                elif any(keyword in attr_name for keyword in ['resistance', 'phenotype', 'susceptibility', 'antibiotic resistance']):
                    resistance_phenotype.append(attr_value)

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
                    antibiotic_resistance.append(resistance_entry)

            if mic_data:
                metadata['mic_data'] = mic_data
            if resistance_phenotype:
                metadata['resistance_phenotype'] = resistance_phenotype
            if antibiotic_resistance:
                metadata['antibiotic_resistance'] = antibiotic_resistance

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

    def _make_request(self, url: str, params: dict, stream: bool = False) -> Optional[requests.Response]:
        """Make HTTP request with error handling"""
        try:
            response = self.session.get(url, params=params, stream=stream, timeout=30)
            response.raise_for_status()
            return response
        except requests.RequestException as e:
            logging.error(f"Request failed: {e}")
            return None