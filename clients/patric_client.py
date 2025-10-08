#!/usr/bin/env python3
"""
PATRIC Client - Pathosystems Resource Integration Center
PATRIC provides comprehensive bacterial genomics data with AMR information
"""

import logging
import time
from typing import List, Dict, Any, Optional

import requests
from tqdm import tqdm

from config import PATRIC_API_KEY, PATRIC_BASE_URL


class PATRICClient:
    """Client for PATRIC API - comprehensive bacterial genomics"""

    def __init__(self, api_key: Optional[str] = PATRIC_API_KEY):
        self.api_key = api_key
        self.base_url = PATRIC_BASE_URL or "https://www.patricbrc.org/api"
        self.session = requests.Session()

        # PATRIC API authentication
        if self.api_key:
            self.session.headers.update({
                'Authorization': f'Bearer {self.api_key}',
                'Accept': 'application/json',
                'Content-Type': 'application/json'
            })
            logging.info("Using PATRIC API key for enhanced access")
        else:
            logging.warning("No PATRIC API key provided - PATRIC requires API key for programmatic access")
            logging.info("Get API key from: https://www.patricbrc.org/")
            logging.info("Without API key, PATRIC searches will be limited or unavailable")

    def fetch_genomes(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """
        Search PATRIC for bacterial genomes

        Args:
            query: Search query (taxon, genome name, etc.)
            max_results: Maximum number of results

        Returns:
            List of genome records with metadata
        """
        logging.info(f"Searching PATRIC with query: {query}")

        try:
            # PATRIC genome search endpoint
            search_url = f"{self.base_url}/genome"

            # PATRIC uses POST requests with JSON payload
            payload = {
                "q": query,
                "limit": max_results,
                "offset": 0,
                "fields": "genome_id,genome_name,organism_name,taxon_id,genome_status,strain,isolation_country,isolation_date,host_name,isolation_source,collection_date,sequencing_centers,assembly_method,genome_length,contigs,n50,gc_content,mlst,amr_profile"
            }

            response = self.session.post(search_url, json=payload, timeout=30)
            response.raise_for_status()

            data = response.json()

            # Parse PATRIC response
            genomes = []
            if isinstance(data, list):
                genomes = data
            elif isinstance(data, dict) and 'response' in data:
                response_data = data['response']
                if isinstance(response_data, list):
                    genomes = response_data
                elif isinstance(response_data, dict) and 'docs' in response_data:
                    genomes = response_data['docs']

            # Convert to harmonized format
            harmonized_genomes = []
            for genome in genomes:
                harmonized = self._parse_patric_genome(genome)
                if harmonized:
                    harmonized_genomes.append(harmonized)

            logging.info(f"Found {len(harmonized_genomes)} genomes in PATRIC")
            return harmonized_genomes

        except Exception as e:
            logging.error(f"Error searching PATRIC: {e}")
            return []

    def download_fasta(self, genome_id: str, output_dir: str, retries: int = 3) -> tuple[bool, str]:
        """Download FASTA sequence from PATRIC"""
        import os
        os.makedirs(output_dir, exist_ok=True)

        try:
            # PATRIC download endpoint
            download_url = f"{self.base_url}/download"

            payload = {
                "genomes": [genome_id],
                "file_format": "fasta",
                "file_type": "contigs"
            }

            filename = f"{genome_id}.fasta"
            filepath = os.path.join(output_dir, filename)

            for attempt in range(retries + 1):
                try:
                    response = self.session.post(download_url, json=payload, stream=True, timeout=60)
                    response.raise_for_status()

                    # Download with progress bar
                    total_size = int(response.headers.get('content-length', 0))
                    with open(filepath, 'wb') as f, tqdm(
                        desc=f"Downloading {genome_id}",
                        total=total_size,
                        unit='B',
                        unit_scale=True,
                        unit_divisor=1024,
                    ) as pbar:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                pbar.update(len(chunk))

                    logging.info(f"Successfully downloaded {genome_id}")
                    return True, genome_id

                except Exception as e:
                    logging.warning(f"Attempt {attempt + 1} failed for {genome_id}: {e}")
                    if attempt < retries:
                        time.sleep(2 ** attempt)

            logging.error(f"Failed to download {genome_id} after {retries + 1} attempts")
            return False, genome_id

        except Exception as e:
            logging.error(f"Error downloading {genome_id}: {e}")
            return False, genome_id

    def _parse_patric_genome(self, genome: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Parse PATRIC genome record into harmonized format"""
        try:
            harmonized = {
                'accession': genome.get('genome_id', ''),
                'organism': genome.get('organism_name', ''),
                'strain': genome.get('strain', ''),
                'database': 'PATRIC',
                'source_url': f"https://www.patricbrc.org/view/Genome/{genome.get('genome_id', '')}",
                'metadata_quality': 'high'
            }

            # Extract location data
            if genome.get('isolation_country'):
                harmonized['country'] = genome['isolation_country']

            # Extract dates
            if genome.get('isolation_date'):
                harmonized['collection_date'] = genome['isolation_date']
            elif genome.get('collection_date'):
                harmonized['collection_date'] = genome['collection_date']

            # Extract host information
            if genome.get('host_name'):
                harmonized['host'] = genome['host_name']

            # Extract isolation source
            if genome.get('isolation_source'):
                harmonized['isolation_source'] = genome['isolation_source']

            # Extract assembly statistics
            if genome.get('genome_length'):
                harmonized['genome_length'] = genome['genome_length']
            if genome.get('contigs'):
                harmonized['contigs'] = genome['contigs']
            if genome.get('n50'):
                harmonized['n50'] = genome['n50']
            if genome.get('gc_content'):
                harmonized['gc_content'] = genome['gc_content']

            # Extract MLST data
            if genome.get('mlst'):
                mlst_data = genome['mlst']
                if isinstance(mlst_data, dict):
                    harmonized['mlst_st'] = mlst_data.get('ST', '')
                    harmonized['mlst_scheme'] = mlst_data.get('scheme', '')

            # Extract AMR profile (PATRIC has excellent AMR data)
            if genome.get('amr_profile'):
                amr_data = genome['amr_profile']
                harmonized.update(self._parse_patric_amr(amr_data))

            # Calculate quality score
            harmonized['quality_score'] = self._calculate_patric_quality_score(harmonized)

            return harmonized

        except Exception as e:
            logging.warning(f"Error parsing PATRIC genome: {e}")
            return None

    def _parse_patric_amr(self, amr_profile: Any) -> Dict[str, Any]:
        """Parse PATRIC AMR profile"""
        amr_data = {
            'mic_data': [],
            'antibiotic_resistance': [],
            'resistance_phenotype': []
        }

        try:
            if isinstance(amr_profile, list):
                for amr_entry in amr_profile:
                    if isinstance(amr_entry, dict):
                        antibiotic = amr_entry.get('antibiotic', '')
                        mic = amr_entry.get('mic', '')
                        resistance = amr_entry.get('resistance', '')

                        if mic:
                            mic_entry = {
                                'antibiotic': antibiotic,
                                'value': str(mic),
                                'unit': amr_entry.get('unit', 'mg/L'),
                                'source': 'PATRIC'
                            }
                            amr_data['mic_data'].append(mic_entry)

                        if resistance:
                            resistance_entry = {
                                'antibiotic': antibiotic,
                                'resistance': resistance,
                                'source': 'PATRIC'
                            }
                            amr_data['antibiotic_resistance'].append(resistance_entry)

                            # Add to phenotype
                            if resistance.lower() in ['resistant', 'intermediate']:
                                amr_data['resistance_phenotype'].append(f"{antibiotic} {resistance}")

        except Exception as e:
            logging.warning(f"Error parsing PATRIC AMR data: {e}")

        return amr_data

    def _calculate_patric_quality_score(self, genome: Dict[str, Any]) -> int:
        """Calculate quality score for PATRIC genomes"""
        score = 6  # Base score for PATRIC (good quality database)

        # AMR data
        if genome.get('mic_data') and len(genome['mic_data']) > 0:
            score += 2
        if genome.get('antibiotic_resistance') and len(genome['antibiotic_resistance']) > 0:
            score += 1

        # Metadata completeness
        if genome.get('country'):
            score += 1
        if genome.get('collection_date'):
            score += 1

        return min(score, 10)

    def search_amr_genomes(self, antibiotic: Optional[str] = None,
                          resistance: Optional[str] = None,
                          max_results: int = 50) -> List[Dict[str, Any]]:
        """Search for genomes with specific AMR patterns"""
        try:
            query_parts = []

            if antibiotic:
                query_parts.append(f"antibiotic:{antibiotic}")
            if resistance:
                query_parts.append(f"resistance:{resistance}")

            query = " AND ".join(query_parts) if query_parts else "genome_status:complete"

            return self.fetch_genomes(query, max_results)

        except Exception as e:
            logging.error(f"Error searching AMR genomes: {e}")
            return []