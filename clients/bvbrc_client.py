#!/usr/bin/env python3
"""
BV-BRC Client - Handles all interactions with BV-BRC (formerly PATRIC) API
"""

import logging
import time
from typing import List, Dict, Any, Optional

import requests
from tqdm import tqdm

from config import BV_BRC_BASE_URL, BV_BRC_API_KEY


class BVBRCClient:
    """Client for BV-BRC API interactions"""

    def __init__(self, api_key: Optional[str] = BV_BRC_API_KEY):
        self.api_key = api_key
        self.base_url = BV_BRC_BASE_URL
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'FederatedGenomeHarvester/1.0',
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        })

        # Add API key to headers if provided
        if self.api_key:
            self.session.headers.update({'Authorization': f'Bearer {self.api_key}'})

    def fetch_genomes(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """
        Search BV-BRC for genomes and return raw records

        Args:
            query: BV-BRC search query
            max_results: Maximum number of results to return

        Returns:
            List of raw genome records from BV-BRC
        """
        logging.info(f"Searching BV-BRC with query: {query}")

        try:
            # BV-BRC uses a POST request with JSON payload for complex queries
            search_payload = {
                "query": query,
                "fields": [
                    "genome_id", "genome_name", "organism", "taxon_id", "genome_status",
                    "strain", "serovar", "biovar", "pathovar", "isolation_country",
                    "isolation_date", "host_name", "host_group", "isolation_source",
                    "collection_date", "completion_date", "sequencing_centers",
                    "assembly_method", "genome_length", "contigs", "n50",
                    "gc_content", "refseq_cds", "patric_cds"
                ],
                "sort": [{"field": "genome_name", "direction": "asc"}],
                "limit": max_results
            }

            response = self._make_request("genome", search_payload)
            if not response:
                return []

            data = response.json()
            genomes = data.get('response', {}).get('docs', [])

            logging.info(f"Found {len(genomes)} genomes in BV-BRC")
            return genomes

        except Exception as e:
            logging.error(f"Error searching BV-BRC: {e}")
            return []

    def download_fasta(self, genome_id: str, output_dir: str,
                      retries: int = 3) -> tuple[bool, str]:
        """Download FASTA sequence for a BV-BRC genome"""
        import os
        os.makedirs(output_dir, exist_ok=True)

        try:
            # For BV-BRC, genome_id is typically the accession, so use it directly
            accession = genome_id

            # BV-BRC provides direct download URLs for genome sequences
            download_url = f"{self.base_url}/download/genome/{genome_id}.fna"

            filename = f"{accession}.fasta"
            filepath = os.path.join(output_dir, filename)

            for attempt in range(retries + 1):
                try:
                    response = self.session.get(download_url, stream=True, timeout=30)
                    response.raise_for_status()

                    # Download with progress bar
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
                        time.sleep(1 * (2 ** attempt))  # Exponential backoff

            logging.error(f"Failed to download {accession} after {retries + 1} attempts")
            return False, accession

        except Exception as e:
            logging.error(f"Error downloading {genome_id}: {e}")
            return False, genome_id

    def _make_request(self, endpoint: str, payload: Optional[Dict[str, Any]] = None,
                     method: str = "POST") -> Optional[requests.Response]:
        """Make HTTP request to BV-BRC API"""
        url = f"{self.base_url}/{endpoint}"

        try:
            if method.upper() == "POST":
                response = self.session.post(url, json=payload, timeout=30)
            else:
                response = self.session.get(url, params=payload, timeout=30)

            response.raise_for_status()
            return response

        except requests.RequestException as e:
            logging.error(f"BV-BRC API request failed: {e}")
            return None
        except Exception as e:
            logging.error(f"Unexpected error in BV-BRC request: {e}")
            return None

    def get_genome_amr_data(self, genome_id: str) -> Dict[str, Any]:
        """Get AMR data for a specific genome"""
        try:
            payload = {
                "query": f"genome_id:{genome_id}",
                "fields": [
                    "genome_id", "genome_name", "antibiotic", "resistant_phenotype",
                    "measurement", "measurement_value", "measurement_sign",
                    "laboratory_typing_method", "testing_standard"
                ]
            }

            response = self._make_request("genome_amr", payload)
            if response:
                data = response.json()
                return data.get('response', {}).get('docs', [])
            return {}

        except Exception as e:
            logging.warning(f"Failed to get AMR data for {genome_id}: {e}")
            return {}

    def search_amr_genomes(self, antibiotic: Optional[str] = None, resistance: Optional[str] = None,
                          max_results: int = 50) -> List[Dict[str, Any]]:
        """Search for genomes with specific AMR patterns"""
        try:
            query_parts = []
            if antibiotic:
                query_parts.append(f"antibiotic:{antibiotic}")
            if resistance:
                query_parts.append(f"resistant_phenotype:{resistance}")

            query = " AND ".join(query_parts) if query_parts else "antibiotic:*"

            payload = {
                "query": query,
                "fields": [
                    "genome_id", "genome_name", "organism", "antibiotic",
                    "resistant_phenotype", "measurement_value"
                ],
                "sort": [{"field": "genome_name", "direction": "asc"}],
                "limit": max_results
            }

            response = self._make_request("genome_amr", payload)
            if response:
                data = response.json()
                docs = data.get('response', {}).get('docs', [])
                logging.info(f"Found {len(docs)} AMR genome records")
                return docs
            return []

        except Exception as e:
            logging.error(f"Error searching AMR genomes: {e}")
            return []