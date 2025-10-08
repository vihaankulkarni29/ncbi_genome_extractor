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
        """Search NCBI for genome IDs matching the query, focusing on complete genomes"""
        # Enhance query to prioritize complete genomes and chromosome-level assemblies
        enhanced_query = f"({query}) AND (complete genome[TI] OR chromosome[TITL] OR \"complete sequence\"[TI]) NOT (scaffold OR contig OR plasmid[TITL])"

        params = {
            'db': 'nuccore',
            'term': enhanced_query,
            'retmax': max_results * 2,  # Get more to filter for quality
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

        # Filter for high-quality complete genomes
        filtered_ids = self._filter_complete_genomes(ids, max_results)
        return filtered_ids[:max_results]

    def _filter_complete_genomes(self, genome_ids: List[str], max_results: int) -> List[str]:
        """Filter genome IDs to prioritize complete genomes with good metadata"""
        if not genome_ids:
            return []

        # Get metadata for all candidates to assess quality
        logging.info(f"Evaluating {len(genome_ids)} genomes for completeness and metadata quality...")

        # Process in smaller batches for quality assessment
        batch_size = 50
        quality_genomes = []

        for i in range(0, len(genome_ids), batch_size):
            batch_ids = genome_ids[i:i + batch_size]
            batch_metadata = self._get_metadata_batch(batch_ids)

            for metadata in batch_metadata:
                # Calculate quality score
                quality_score = self._calculate_metadata_score(metadata)

                # Check if it's a complete genome (not scaffold/contig)
                accession = metadata.get('accession', '')
                title = metadata.get('title', '').lower()

                # Prioritize complete genomes
                is_complete = (
                    'chromosome' in title or
                    'complete genome' in title or
                    'complete sequence' in title or
                    accession.startswith('CP') or  # Complete genome prefix
                    accession.startswith('NC_') or # RefSeq complete
                    accession.startswith('NZ_CP')  # Complete bacterial genomes
                )

                # Avoid scaffolds and contigs
                is_not_fragment = not (
                    'scaffold' in title or
                    'contig' in title or
                    'plasmid' in title or
                    accession.startswith('NZ_') and not accession.startswith('NZ_CP')
                )

                if is_complete and is_not_fragment:
                    quality_genomes.append({
                        'id': metadata.get('genome_id'),
                        'accession': accession,
                        'quality_score': quality_score,
                        'has_biosample': bool(metadata.get('biosample')),
                        'has_collection_data': bool(metadata.get('collection_date') or metadata.get('country')),
                        'metadata': metadata
                    })

        # Sort by quality: complete genomes with metadata first
        quality_genomes.sort(key=lambda x: (
            x['quality_score'],  # Primary: quality score
            x['has_biosample'],  # Secondary: has BioSample
            x['has_collection_data'],  # Tertiary: has collection data
            x['accession'].startswith('CP')  # Prefer complete genome accessions
        ), reverse=True)

        selected_ids = [genome['id'] for genome in quality_genomes[:max_results]]

        high_quality = sum(1 for g in quality_genomes[:max_results] if g['quality_score'] >= 7)
        medium_quality = sum(1 for g in quality_genomes[:max_results] if 4 <= g['quality_score'] < 7)
        low_quality = sum(1 for g in quality_genomes[:max_results] if g['quality_score'] < 4)

        logging.info(f"Selected {len(selected_ids)} high-quality complete genomes:")
        logging.info(f"  - High quality (score ≥7): {high_quality}")
        logging.info(f"  - Medium quality (score 4-6): {medium_quality}")
        logging.info(f"  - Low quality (score <4): {low_quality}")

        return selected_ids

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
                # Update metadata with BioSample data
                metadata.update(biosample_metadata)
                logging.debug(f"Enhanced metadata with BioSample data for {metadata['accession']}")

        # If no BioSample, try to extract metadata from title and other fields
        if not metadata.get('country') and not metadata.get('collection_date'):
            self._extract_metadata_from_title(metadata)

        # Try to get assembly information which might have BioSample data
        if not metadata.get('biosample'):
            assembly_biosample = self._get_biosample_from_assembly(metadata.get('accession'))
            if assembly_biosample:
                metadata['biosample'] = assembly_biosample
                # Try to get BioSample metadata again
                biosample_metadata = self._get_biosample_metadata(assembly_biosample)
                if biosample_metadata:
                    metadata.update(biosample_metadata)

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
        """Extract BioSample and BioProject IDs from extra field with enhanced patterns"""
        if not extra:
            return

        # Look for BioSample - comprehensive patterns
        biosample_patterns = [
            r'SAM[NED]?\d+',  # Standard format: SAMN123, SAME123, SAMD123
            r'SAM[NED][A-Z]?\d+',  # With optional letter: SAMNA123, SAMEA123
            r'SAM[NED]\d+[A-Z]?',  # Number then optional letter: SAMN123A
            r'SAM[NED][A-Z]+\d+',  # Letters then number: SAMNA123
            r'SAM[NED]\d+[A-Z]+\d*',  # Complex: SAMN123A456
        ]

        for pattern in biosample_patterns:
            biosample_match = re.search(pattern, extra, re.IGNORECASE)
            if biosample_match:
                biosample_id = biosample_match.group(0).upper()
                # Validate format
                if re.match(r'SAM[NED][A-Z]*\d+', biosample_id):
                    metadata['biosample'] = biosample_id
                    break

        # Look for BioProject - comprehensive patterns
        bioproject_patterns = [
            r'PRJ[NED][A-Z]\d+',  # Standard: PRJNA123, PRJEB123
            r'PRJ[NED]\d+',       # Without letter: PRJN123
            r'PRJ[A-Z]+\d+',      # Any letters: PRJNA123, PRJEB123
            r'PRJ[A-Z]*\d+',      # Flexible: PRJ123, PRJNA123
        ]

        for pattern in bioproject_patterns:
            bioproject_match = re.search(pattern, extra, re.IGNORECASE)
            if bioproject_match:
                bioproject_id = bioproject_match.group(0).upper()
                # Validate format
                if re.match(r'PRJ[A-Z]*\d+', bioproject_id):
                    metadata['bioproject'] = bioproject_id
                    break

        # Look for SRA/ENA accession which might link to BioSample
        sra_patterns = [
            r'SRR\d+',  # SRA run: SRR123456
            r'SRX\d+',  # SRA experiment: SRX123456
            r'SRS\d+',  # SRA sample: SRS123456
            r'ERP\d+',  # ENA project: ERP123456
            r'ERS\d+',  # ENA sample: ERS123456
            r'ERR\d+',  # ENA run: ERR123456
        ]

        for pattern in sra_patterns:
            sra_match = re.search(pattern, extra, re.IGNORECASE)
            if sra_match:
                metadata['sra_accession'] = sra_match.group(0).upper()
                break

        # Also try to find IDs in other fields if not found in Extra
        if not metadata.get('biosample'):
            # Check title for BioSample/Project IDs
            title = metadata.get('title', '')
            if title:
                for pattern in biosample_patterns:
                    match = re.search(pattern, title, re.IGNORECASE)
                    if match:
                        biosample_id = match.group(0).upper()
                        if re.match(r'SAM[NED][A-Z]*\d+', biosample_id):
                            metadata['biosample'] = biosample_id
                            break

                if not metadata.get('bioproject'):
                    for pattern in bioproject_patterns:
                        match = re.search(pattern, title, re.IGNORECASE)
                        if match:
                            bioproject_id = match.group(0).upper()
                            if re.match(r'PRJ[A-Z]*\d+', bioproject_id):
                                metadata['bioproject'] = bioproject_id
                                break

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
        """Parse BioSample XML for comprehensive AMR-related metadata"""
        metadata = {}

        try:
            root = ET.fromstring(xml_content)

            # Extract basic BioSample information
            biosample_elem = root.find('.//BioSample')
            if biosample_elem is not None:
                # Get BioSample accession
                accession = biosample_elem.get('accession')
                if accession:
                    metadata['biosample_accession'] = accession

                # Get BioSample description
                description = biosample_elem.find('.//Description')
                if description is not None and description.text:
                    metadata['biosample_description'] = description.text.strip()

            # Extract collection date (comprehensive patterns)
            collection_date_attrs = [
                'collection_date', 'collection date', 'date', 'collection-date',
                'sampling_date', 'sampling date', 'isolation_date', 'isolation date'
            ]

            for attr_name in collection_date_attrs:
                collection_date = root.find(f'.//Attribute[@attribute_name="{attr_name}"]')
                if collection_date is not None and collection_date.text:
                    date_text = collection_date.text.strip()
                    # Try to standardize date format
                    metadata['collection_date'] = self._standardize_date(date_text)
                    break

            # Extract country/geographic location (comprehensive patterns)
            location_attrs = [
                'geo_loc_name', 'geographic location', 'country', 'location',
                'geographic_location', 'geo_loc', 'isolation_country', 'origin_country'
            ]

            for attr_name in location_attrs:
                country = root.find(f'.//Attribute[@attribute_name="{attr_name}"]')
                if country is not None and country.text:
                    location_text = country.text.strip()
                    metadata['country'] = self._standardize_location(location_text)
                    break

            # Extract host information (comprehensive patterns)
            host_attrs = [
                'host', 'host organism', 'source host', 'host_species',
                'host_name', 'host_taxon', 'host_organism'
            ]

            for attr_name in host_attrs:
                host = root.find(f'.//Attribute[@attribute_name="{attr_name}"]')
                if host is not None and host.text:
                    metadata['host'] = host.text.strip()
                    break

            # Extract isolation source (comprehensive patterns)
            source_attrs = [
                'isolation_source', 'isolation source', 'source', 'sample type',
                'source_material', 'material', 'sample_type', 'specimen_type'
            ]

            for attr_name in source_attrs:
                isolation = root.find(f'.//Attribute[@attribute_name="{attr_name}"]')
                if isolation is not None and isolation.text:
                    metadata['isolation_source'] = isolation.text.strip()
                    break

            # Extract clinical/disease information
            clinical_attrs = [
                'disease', 'clinical_disease', 'diagnosis', 'condition',
                'health_status', 'disease_state', 'clinical_condition'
            ]

            for attr_name in clinical_attrs:
                disease = root.find(f'.//Attribute[@attribute_name="{attr_name}"]')
                if disease is not None and disease.text:
                    metadata['disease'] = disease.text.strip()
                    break

            # Extract MIC data and resistance information (highly expanded)
            mic_data = []
            resistance_phenotype = []
            antibiotic_resistance = []

            for attr in root.findall('.//Attribute'):
                attr_name = attr.get('attribute_name', '').lower().strip()
                attr_value = attr.text or ''

                if not attr_value.strip():
                    continue

                attr_value_clean = attr_value.strip()

                # MIC data patterns (highly expanded)
                mic_keywords = [
                    'mic', 'minimum inhibitory concentration', 'mic_', 'mic50', 'mic90',
                    'mic_value', 'mic_result', 'minimum_inhibitory_concentration'
                ]

                if any(keyword in attr_name for keyword in mic_keywords):
                    antibiotic_name = self._extract_antibiotic_from_mic_name(attr_name)
                    if antibiotic_name:
                        mic_data_entry = {
                            'antibiotic': antibiotic_name,
                            'value': attr_value_clean,
                            'unit': self._extract_mic_unit(attr_value_clean),
                            'source': 'biosample'
                        }
                        mic_data.append(mic_data_entry)

                # Resistance phenotype (highly expanded)
                resistance_keywords = [
                    'resistance', 'phenotype', 'susceptibility', 'antibiotic resistance',
                    'resistance_phenotype', 'susceptibility_profile', 'resistance_profile',
                    'antibiotic_susceptibility', 'drug_resistance'
                ]

                if any(keyword in attr_name for keyword in resistance_keywords):
                    resistance_phenotype.append(attr_value_clean)

                # Specific antibiotic resistance (comprehensive list)
                antibiotic_list = [
                    # Beta-lactams
                    'ampicillin', 'amoxicillin', 'penicillin', 'oxacillin', 'methicillin',
                    'cefazolin', 'cefotaxime', 'ceftriaxone', 'ceftazidime', 'cefepime',
                    'meropenem', 'imipenem', 'ertapenem', 'doripenem', 'aztreonam',

                    # Aminoglycosides
                    'gentamicin', 'tobramycin', 'amikacin', 'kanamycin', 'streptomycin',
                    'netilmicin', 'neomycin',

                    # Fluoroquinolones
                    'ciprofloxacin', 'levofloxacin', 'moxifloxacin', 'norfloxacin',
                    'ofloxacin', 'gemifloxacin',

                    # Tetracyclines
                    'tetracycline', 'doxycycline', 'minocycline', 'tigecycline',

                    # Sulfonamides
                    'trimethoprim', 'sulfamethoxazole', 'trimethoprim_sulfamethoxazole',

                    # Macrolides
                    'erythromycin', 'azithromycin', 'clarithromycin', 'telithromycin',

                    # Lincosamides
                    'clindamycin', 'lincomycin',

                    # Glycopeptides
                    'vancomycin', 'teicoplanin',

                    # Oxazolidinones
                    'linezolid',

                    # Lipopeptides
                    'daptomycin',

                    # Polymyxins
                    'colistin', 'polymyxin_b',

                    # Others
                    'chloramphenicol', 'rifampin', 'fusidic_acid', 'fosfomycin',
                    'nitrofurantoin', 'tigecycline',

                    # Resistance markers
                    'esbl', 'mrsa', 'vre', 'carbapenemase', 'extended_spectrum_beta_lactamase',
                    'methicillin_resistant', 'vancomycin_resistant', 'carbapenem_resistant'
                ]

                # Check if attribute name contains any antibiotic
                for antibiotic in antibiotic_list:
                    if antibiotic.replace('_', ' ') in attr_name or antibiotic in attr_name:
                        resistance_entry = {
                            'antibiotic': antibiotic.replace('_', ' '),
                            'resistance': attr_value_clean,
                            'source': 'biosample'
                        }
                        antibiotic_resistance.append(resistance_entry)
                        break

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

    def _extract_metadata_from_title(self, metadata: Dict[str, Any]):
        """Extract additional metadata from sequence title when BioSample is not available"""
        title = metadata.get('title', '')
        if not title:
            return

        # Extract collection year from title
        year_patterns = [
            r'isolat\w*\s+(\d{4})',  # isolate 2020
            r'(\d{4})[\s\-]isolat',  # 2020-isolate
            r'strain.*(\d{4})',      # strain ... 2020
            r'(\d{4})[\s\-]strain',  # 2020-strain
        ]

        for pattern in year_patterns:
            match = re.search(pattern, title, re.IGNORECASE)
            if match:
                year = match.group(1)
                if 1900 <= int(year) <= 2030:  # Reasonable year range
                    metadata['collection_date'] = f"{year}-01-01"  # Default to Jan 1st
                    break

        # Extract country/location information from title
        location_patterns = [
            r'(\b\w+),\s*India',  # City, India
            r'India[\s\-](\w+)',  # India-City
            r'(\w+)\s*India',     # City India
        ]

        for pattern in location_patterns:
            match = re.search(pattern, title, re.IGNORECASE)
            if match:
                location = match.group(1)
                if len(location) > 2:  # Avoid abbreviations
                    metadata['country'] = f"{location}, India"
                    break

        # If we found India but no specific location
        if 'india' in title.lower() and not metadata.get('country'):
            metadata['country'] = 'India'

    def _get_biosample_from_assembly(self, accession: str) -> Optional[str]:
        """Try to get BioSample ID from assembly information"""
        if not accession:
            return None

        try:
            # Query assembly database for this accession
            assembly_url = NCBI_BASE_URL + "elink.fcgi"
            params = {
                'dbfrom': 'nuccore',
                'db': 'assembly',
                'id': accession,
                'email': self.email,
                'api_key': self.api_key if self.api_key else None
            }

            params = {k: v for k, v in params.items() if v is not None}
            response = self._make_request(assembly_url, params)

            if response:
                # Parse the elink response to get assembly IDs
                root = ET.fromstring(response.text)
                assembly_ids = []
                for link in root.findall('.//Link/Id'):
                    assembly_ids.append(link.text)

                if assembly_ids:
                    # Get assembly details which might include BioSample
                    summary_url = NCBI_BASE_URL + "esummary.fcgi"
                    summary_params = {
                        'db': 'assembly',
                        'id': ','.join(assembly_ids[:5]),  # Limit to first 5
                        'email': self.email,
                        'api_key': self.api_key if self.api_key else None
                    }

                    summary_response = self._make_request(summary_url, dict(summary_params))
                    if summary_response:
                        summary_root = ET.fromstring(summary_response.text)
                        for docsum in summary_root.findall('.//DocSum'):
                            # Look for BioSample in assembly metadata
                            for item in docsum.findall('Item'):
                                name = item.get('Name')
                                if name == 'BioSampleAccn' or name == 'BioSampleId':
                                    biosample_id = item.text
                                    if biosample_id:
                                        return biosample_id

        except Exception as e:
            logging.debug(f"Failed to get BioSample from assembly for {accession}: {e}")

        return None

    def _extract_antibiotic_from_mic_name(self, mic_attr_name: str) -> str:
        """Extract antibiotic name from MIC attribute name"""
        # Remove common MIC prefixes/suffixes
        clean_name = mic_attr_name.lower()
        clean_name = re.sub(r'^mic_?', '', clean_name)
        clean_name = re.sub(r'_?mic$', '', clean_name)
        clean_name = re.sub(r'minimum inhibitory concentration', '', clean_name)
        clean_name = re.sub(r'mic\d+', '', clean_name)  # Remove mic50, mic90, etc.
        clean_name = clean_name.strip('_').strip()

        # Convert underscores to spaces
        antibiotic_name = clean_name.replace('_', ' ')

        # If we end up with empty or just numbers, return unknown
        if not antibiotic_name or antibiotic_name.isdigit():
            return 'unknown'

        return antibiotic_name

    def _standardize_date(self, date_text: str) -> str:
        """Standardize date format to YYYY-MM-DD when possible"""
        if not date_text:
            return date_text

        # Try to parse various date formats
        date_patterns = [
            (r'(\d{4})[\-/](\d{1,2})[\-/](\d{1,2})', r'\1-\2-\3'),  # YYYY-MM-DD or YYYY/MM/DD
            (r'(\d{1,2})[\-/](\d{1,2})[\-/](\d{4})', r'\3-\1-\2'),  # MM-DD-YYYY or MM/DD/YYYY
            (r'(\d{4})[\-/](\d{1,2})', r'\1-\2-01'),               # YYYY-MM (assume 1st)
            (r'(\d{4})', r'\1-01-01'),                             # YYYY (assume Jan 1st)
        ]

        for pattern, replacement in date_patterns:
            match = re.search(pattern, date_text)
            if match:
                try:
                    formatted_date = re.sub(pattern, replacement, date_text)
                    # Validate the date
                    parts = formatted_date.split('-')
                    if len(parts) == 3:
                        year, month, day = map(int, parts)
                        if 1900 <= year <= 2030 and 1 <= month <= 12 and 1 <= day <= 31:
                            return formatted_date
                except (ValueError, IndexError):
                    continue

        # If we can't standardize, return original
        return date_text

    def _standardize_location(self, location_text: str) -> str:
        """Standardize location/country names"""
        if not location_text:
            return location_text

        # Common location standardization
        location_lower = location_text.lower().strip()

        # Handle "India" variations
        if 'india' in location_lower:
            return 'India'

        # Handle common country name variations
        country_mappings = {
            'usa': 'USA',
            'united states': 'USA',
            'united states of america': 'USA',
            'uk': 'United Kingdom',
            'united kingdom': 'United Kingdom',
            'china': 'China',
            'japan': 'Japan',
            'germany': 'Germany',
            'france': 'France',
            'italy': 'Italy',
            'spain': 'Spain',
            'canada': 'Canada',
            'australia': 'Australia',
            'brazil': 'Brazil',
            'mexico': 'Mexico',
            'south korea': 'South Korea',
            'korea': 'South Korea',
        }

        for key, value in country_mappings.items():
            if key in location_lower:
                return value

        # If no mapping found, return title case
        return location_text.strip().title()

    def _extract_mic_unit(self, mic_value: str) -> str:
        """Extract MIC unit from MIC value string"""
        # Common MIC units
        units = ['ug/ml', 'μg/ml', 'mg/l', 'mcg/ml', 'µg/ml', 'μg/l', 'mcg/l']
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