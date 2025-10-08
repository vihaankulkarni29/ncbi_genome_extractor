#!/usr/bin/env python3
"""
EnteroBase Client - Specialized for E. coli and Salmonella genomics
EnteroBase provides rich AMR metadata and strain typing data
"""

import logging
import time
from typing import List, Dict, Any, Optional

import requests
from tqdm import tqdm

from config import ENTEROBASE_API_KEY, ENTEROBASE_BASE_URL


class EnteroBaseClient:
    """Client for EnteroBase API - specialized for E. coli and Salmonella"""

    def __init__(self, api_key: Optional[str] = ENTEROBASE_API_KEY):
        self.api_key = api_key
        self.base_url = ENTEROBASE_BASE_URL or "https://enterobase.warwick.ac.uk"
        self.session = requests.Session()

        # EnteroBase API key in query parameters
        self.session.headers.update({
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        })

        if not self.api_key:
            logging.warning("No EnteroBase API key provided - EnteroBase requires API key for programmatic access")
            logging.info("Get API key from: https://enterobase.warwick.ac.uk/enterobase/")
            logging.info("Without API key, EnteroBase searches will be limited or unavailable")

    def fetch_genomes(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """
        Search EnteroBase for E. coli/Salmonella genomes

        Args:
            query: Search query (strain names, ST types, etc.)
            max_results: Maximum number of results

        Returns:
            List of genome records with rich metadata
        """
        logging.info(f"Searching EnteroBase with query: {query}")

        try:
            # EnteroBase search endpoint
            search_url = f"{self.base_url}/search"

            # Build search parameters
            params = {
                'query': query,
                'species': 'Escherichia coli',  # Focus on E. coli
                'limit': max_results,
                'offset': 0,
                'fields': 'strain_name,lab_contact,owner,created,modified,version,download_accession,serotype,mlst,ngstar,core_genome_mlST,amr_profile,assembly_stats,experimental_data'
            }

            response = self.session.get(search_url, params=params, timeout=30)
            response.raise_for_status()

            data = response.json()

            # Parse EnteroBase response
            genomes = []
            if 'results' in data:
                for result in data['results']:
                    genome = self._parse_enterobase_genome(result)
                    if genome:
                        genomes.append(genome)

            logging.info(f"Found {len(genomes)} genomes in EnteroBase")
            return genomes

        except Exception as e:
            logging.error(f"Error searching EnteroBase: {e}")
            return []

    def download_fasta(self, genome_id: str, output_dir: str, retries: int = 3) -> tuple[bool, str]:
        """Download FASTA sequence from EnteroBase"""
        import os
        os.makedirs(output_dir, exist_ok=True)

        try:
            # EnteroBase download endpoint
            download_url = f"{self.base_url}/download"

            params = {
                'accession': genome_id,
                'type': 'assembly',
                'format': 'fasta'
            }

            filename = f"{genome_id}.fasta"
            filepath = os.path.join(output_dir, filename)

            for attempt in range(retries + 1):
                try:
                    response = self.session.get(download_url, params=params, stream=True, timeout=60)
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

    def _parse_enterobase_genome(self, result: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Parse EnteroBase genome record into harmonized format"""
        try:
            genome = {
                'accession': result.get('download_accession', ''),
                'organism': 'Escherichia coli',  # EnteroBase focuses on specific species
                'strain': result.get('strain_name', ''),
                'database': 'EnteroBase',
                'source_url': f"https://enterobase.warwick.ac.uk/species/index/ecoli",
                'metadata_quality': 'high'  # EnteroBase has rich metadata
            }

            # Extract collection/lab information
            if 'lab_contact' in result:
                genome['lab_contact'] = result['lab_contact']

            if 'owner' in result:
                genome['owner'] = result['owner']

            # Extract dates
            if 'created' in result:
                genome['created_date'] = result['created']

            if 'modified' in result:
                genome['modified_date'] = result['modified']

            # Extract serotype information
            if 'serotype' in result:
                genome['serotype'] = result['serotype']

            # Extract MLST information
            if 'mlst' in result:
                mlst_data = result['mlst']
                if isinstance(mlst_data, dict):
                    genome['mlst_st'] = mlst_data.get('ST', '')
                    genome['mlst_profile'] = mlst_data.get('profile', '')

            # Extract AMR profile (this is the gold standard!)
            if 'amr_profile' in result:
                amr_data = result['amr_profile']
                genome.update(self._parse_amr_profile(amr_data))

            # Extract assembly statistics
            if 'assembly_stats' in result:
                assembly_stats = result['assembly_stats']
                if isinstance(assembly_stats, dict):
                    genome['genome_length'] = assembly_stats.get('total_length', 0)
                    genome['contigs'] = assembly_stats.get('contigs', 0)
                    genome['n50'] = assembly_stats.get('N50', 0)
                    genome['gc_content'] = assembly_stats.get('gc_content', 0)

            # Calculate quality score
            genome['quality_score'] = self._calculate_enterobase_quality_score(genome)

            return genome

        except Exception as e:
            logging.warning(f"Error parsing EnteroBase genome: {e}")
            return None

    def _parse_amr_profile(self, amr_profile: Dict[str, Any]) -> Dict[str, Any]:
        """Parse EnteroBase AMR profile into standardized format"""
        amr_data = {
            'mic_data': [],
            'antibiotic_resistance': [],
            'resistance_phenotype': []
        }

        try:
            if isinstance(amr_profile, dict):
                # EnteroBase provides detailed MIC and resistance data
                for antibiotic, data in amr_profile.items():
                    if isinstance(data, dict):
                        # MIC value
                        if 'mic' in data:
                            mic_entry = {
                                'antibiotic': antibiotic,
                                'value': str(data['mic']),
                                'unit': data.get('unit', 'mg/L'),
                                'source': 'EnteroBase'
                            }
                            amr_data['mic_data'].append(mic_entry)

                        # Resistance interpretation
                        if 'resistance' in data:
                            resistance_entry = {
                                'antibiotic': antibiotic,
                                'resistance': data['resistance'],
                                'source': 'EnteroBase'
                            }
                            amr_data['antibiotic_resistance'].append(resistance_entry)

                            # Add to phenotype list
                            if data['resistance'].lower() in ['resistant', 'intermediate']:
                                amr_data['resistance_phenotype'].append(f"{antibiotic} {data['resistance']}")

        except Exception as e:
            logging.warning(f"Error parsing AMR profile: {e}")

        return amr_data

    def _calculate_enterobase_quality_score(self, genome: Dict[str, Any]) -> int:
        """Calculate quality score for EnteroBase genomes (0-10 scale)"""
        score = 5  # Base score for EnteroBase (high quality database)

        # AMR data (most important)
        if genome.get('mic_data') and len(genome['mic_data']) > 0:
            score += 2
        if genome.get('antibiotic_resistance') and len(genome['antibiotic_resistance']) > 0:
            score += 1

        # MLST data
        if genome.get('mlst_st'):
            score += 1

        # Assembly quality
        if genome.get('genome_length', 0) > 4000000:  # Complete E. coli genome
            score += 1

        return min(score, 10)

    def search_amr_genomes(self, antibiotic: Optional[str] = None,
                          resistance: Optional[str] = None,
                          max_results: int = 50) -> List[Dict[str, Any]]:
        """Search for genomes with specific AMR patterns"""
        try:
            query_parts = ["species:Escherichia coli"]

            if antibiotic:
                query_parts.append(f"amr:{antibiotic}")
            if resistance:
                query_parts.append(f"resistance:{resistance}")

            query = " AND ".join(query_parts)

            return self.fetch_genomes(query, max_results)

        except Exception as e:
            logging.error(f"Error searching AMR genomes: {e}")
            return []