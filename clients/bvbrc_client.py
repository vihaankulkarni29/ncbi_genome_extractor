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
        # BV-BRC API endpoints - trying the correct service endpoints
        self.base_url = "https://www.bv-brc.org"
        self.api_base = "https://www.bv-brc.org/api"
        self.services_base = "https://www.bv-brc.org/services"

        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'FederatedGenomeHarvester/1.0',
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        })

        if not self.api_key:
            logging.info("No BV-BRC API key provided - using public access (limited functionality)")
        else:
            logging.info("Using BV-BRC API key for enhanced access")

    def fetch_genomes(self, query: str, max_results: int = 100) -> List[Dict[str, Any]]:
        """
        Search BV-BRC for genomes using multiple API endpoint strategies

        Args:
            query: BV-BRC search query
            max_results: Maximum number of results to return

        Returns:
            List of raw genome records from BV-BRC
        """
        logging.info(f"Searching BV-BRC with query: {query}")

        # Try multiple API endpoint strategies
        strategies = [
            self._try_solr_api,
            self._try_service_api,
            self._try_rest_api
        ]

        for strategy in strategies:
            try:
                genomes = strategy(query, max_results)
                if genomes:
                    logging.info(f"Successfully found {len(genomes)} genomes using {strategy.__name__}")
                    return genomes
            except Exception as e:
                logging.debug(f"Strategy {strategy.__name__} failed: {e}")
                continue

        logging.warning("All BV-BRC API strategies failed")
        return []

    def _try_solr_api(self, query: str, max_results: int) -> List[Dict[str, Any]]:
        """Try Solr-based API endpoint"""
        search_url = f"{self.api_base}/genome"

        params = {
            'q': query,
            'fl': 'genome_id,genome_name,organism_name,taxon_id,genome_status,strain,isolation_country,isolation_date,host_name,isolation_source,collection_date',
            'sort': 'genome_name asc',
            'rows': str(max_results),
            'wt': 'json'
        }

        if self.api_key:
            params['auth_token'] = self.api_key

        response = self.session.get(search_url, params=params, timeout=30)
        response.raise_for_status()

        data = response.json()
        return self._extract_genomes_from_response(data)

    def _try_service_api(self, query: str, max_results: int) -> List[Dict[str, Any]]:
        """Try service-based API endpoint"""
        search_url = f"{self.services_base}/Genome"

        payload = {
            "method": "query",
            "params": [{
                "q": query,
                "limit": max_results,
                "offset": 0
            }],
            "id": 1,
            "jsonrpc": "2.0"
        }

        if self.api_key:
            payload['params'][0]['auth_token'] = self.api_key

        response = self.session.post(search_url, json=payload, timeout=30)
        response.raise_for_status()

        data = response.json()
        return self._extract_genomes_from_response(data)

    def _try_rest_api(self, query: str, max_results: int) -> List[Dict[str, Any]]:
        """Try REST API endpoint"""
        search_url = f"{self.base_url}/api/genome/search"

        payload = {
            "query": query,
            "limit": max_results,
            "fields": ["genome_id", "genome_name", "organism_name", "isolation_country"]
        }

        if self.api_key:
            payload['auth_token'] = self.api_key

        response = self.session.post(search_url, json=payload, timeout=30)
        response.raise_for_status()

        data = response.json()
        return self._extract_genomes_from_response(data)

    def _extract_genomes_from_response(self, data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Extract genome list from various response formats"""
        genomes = []

        # Try different response structures
        if isinstance(data, list):
            genomes = data
        elif isinstance(data, dict):
            # Check for various possible structures
            for key in ['result', 'response', 'docs', 'genomes', 'data']:
                if key in data:
                    candidate = data[key]
                    if isinstance(candidate, list):
                        genomes = candidate
                        break
                    elif isinstance(candidate, dict):
                        # Check nested structures
                        for subkey in ['docs', 'genomes', 'result']:
                            if subkey in candidate and isinstance(candidate[subkey], list):
                                genomes = candidate[subkey]
                                break
                        if genomes:
                            break

        return genomes

    def download_fasta(self, genome_id: str, output_dir: str,
                      retries: int = 3) -> tuple[bool, str]:
        """Download FASTA sequence for a BV-BRC genome"""
        import os
        os.makedirs(output_dir, exist_ok=True)

        try:
            # For BV-BRC, genome_id is typically the accession, so use it directly
            accession = genome_id

            # BV-BRC provides direct download URLs for genome sequences
            download_url = f"{self.base_url}/download/genome/{genome_id}.fasta"

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

    def _make_request(self, endpoint: str, params: Optional[Dict[str, Any]] = None,
                     method: str = "GET") -> Optional[requests.Response]:
        """Make HTTP request to BV-BRC API"""
        # BV-BRC API base URL construction
        if endpoint:
            url = f"{self.base_url}/{endpoint}"
        else:
            url = self.base_url.rstrip('/')

        try:
            if method.upper() == "POST":
                # For POST requests, use JSON payload
                response = self.session.post(url, json=params, timeout=30)
            else:
                # For GET requests, use query parameters
                response = self.session.get(url, params=params, timeout=30)

            response.raise_for_status()
            return response

        except requests.RequestException as e:
            logging.error(f"BV-BRC API request failed: {e}")
            logging.error(f"URL: {url}")
            if params:
                logging.error(f"Params: {params}")
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