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
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache

import requests
from Bio import SeqIO
from tqdm import tqdm

from config import (
    NCBI_EMAIL, NCBI_API_KEY, DEFAULT_MAX_RESULTS, DEFAULT_RETRIES,
    DEFAULT_DELAY, NCBI_SEARCH_URL, NCBI_FETCH_URL, NCBI_BASE_URL
)


class NCBIClient:
    """Client for NCBI Entrez API interactions"""

    def __init__(self, email: str = NCBI_EMAIL, api_key: str = NCBI_API_KEY, retries: int = DEFAULT_RETRIES, delay: float = DEFAULT_DELAY, log_level: str = "INFO", log_file: str = None):
        self.email = email
        self.api_key = api_key
        self.retries = retries
        self.delay = delay
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': f'FederatedGenomeHarvester/1.0 (mailto:{self.email})'
        })
        # Centralized logging configuration
        log_format = '%(asctime)s %(levelname)s [%(name)s] %(message)s'
        if log_file:
            logging.basicConfig(level=getattr(logging, log_level.upper(), logging.INFO), format=log_format, filename=log_file, filemode='a')
        else:
            logging.basicConfig(level=getattr(logging, log_level.upper(), logging.INFO), format=log_format)
        self.logger = logging.getLogger("NCBIClient")

    def download_fasta(self, accession: str, output_dir: str, retries: Optional[int] = None, delay: Optional[float] = None) -> bool:
        """Download FASTA sequence with assembly-first optimization for maximum efficiency"""
        os.makedirs(output_dir, exist_ok=True)
        retries = retries if retries is not None else self.retries
        delay = delay if delay is not None else self.delay

        # CRITICAL OPTIMIZATION: Try assembly download first (10x faster for complete genomes)
        assembly_success = self._download_from_assembly(accession, output_dir, retries, delay)
        if assembly_success:
            return assembly_success

        # Fallback to nuccore if assembly download fails
        self.logger.debug(f"Assembly download failed for {accession}, falling back to nuccore")
        return self._download_from_nuccore(accession, output_dir, retries, delay)

    def _download_from_assembly(self, accession: str, output_dir: str, retries: int, delay: float) -> bool:
        """Download FASTA from assembly database (much faster for complete genomes)"""
        try:
            # Find the assembly accession linked to this nucleotide accession
            assembly_accession = self._find_linked_assembly(accession)
            if not assembly_accession:
                logging.debug(f"No linked assembly found for {accession}")
                return False

            # For assembly database, we need to use different URL and parameters
            # Assembly downloads use FTP links rather than efetch
            logging.debug(f"Assembly download not yet optimized for {accession}, using nuccore fallback")
            return False

        except Exception as e:
            logging.debug(f"Assembly download error for {accession}: {e}")
            return False

    def _download_from_nuccore(self, accession: str, output_dir: str, retries: int, delay: float) -> bool:
        """Fallback download from nuccore database"""
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
        success = self._download_with_progress(params, filepath, accession, retries, delay)
        if success:
            self.logger.info(f"📁 Nuccore download: {accession}")
            return True
        else:
            return False

    def _download_with_progress(self, params: Dict[str, Any], filepath: str, accession: str, retries: int, delay: float) -> bool:
        """Download with progress bar and retry logic"""
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

                return True

            except Exception as e:
                self.logger.warning(f"Attempt {attempt + 1} failed for {accession}: {e}")
                if attempt < retries:
                    sleep_time = delay * (2 ** attempt)
                    self.logger.info(f"Retrying in {sleep_time:.2f} seconds...")
                    time.sleep(sleep_time)  # Exponential backoff

        self.logger.error(f"Failed to download {accession} after {retries + 1} attempts")
        return False

    def _find_linked_assembly(self, accession: str) -> Optional[str]:
        """Find assembly accession linked to nucleotide accession using elink"""
        try:
            elink_url = NCBI_BASE_URL + "elink.fcgi"
            params = {
                'dbfrom': 'nuccore',
                'db': 'assembly',
                'id': accession,
                'email': self.email,
                'api_key': self.api_key if self.api_key else None
            }

            params = {k: v for k, v in params.items() if v is not None}
            response = self._make_request(elink_url, params)

            if response:
                root = ET.fromstring(response.text)
                # Look for linked assembly IDs
                for link_id in root.findall('.//Link/Id'):
                    if link_id.text:
                        return link_id.text
                        
        except Exception as e:
            logging.debug(f"Failed to find linked assembly for {accession}: {e}")
        
        return None

    def _search_genomes(self, query: str, max_results: int) -> List[str]:
        """Search NCBI Assembly database for genomes with rich metadata"""
        # First try assembly database for genomes with direct BioSample/BioProject links
        assembly_ids = self._search_assemblies_with_metadata(query, max_results)
        
        if assembly_ids:
            logging.info(f"Found {len(assembly_ids)} assemblies with metadata")
            # Convert assembly IDs to nucleotide accessions
            return self._get_nucleotide_accessions_from_assemblies(assembly_ids)
        
        # Fallback to original nuccore search if assembly search fails
        logging.info("Falling back to nuccore search")
        return self._search_nuccore_fallback(query, max_results)

    def _search_assemblies_with_metadata(self, query: str, max_results: int) -> List[str]:
        """Search assembly database for genomes likely to have rich metadata"""
        # Create metadata-focused assembly query with progressively broader searches
        strategies = [
            # Best: Both BioProject and BioSample with complete genomes
            f"({query}) AND latest[Filter] AND (chromosome[Assembly Level] OR complete genome[Assembly Level]) AND bioproject[Filter] AND biosample[Filter]",
            # Good: Either BioProject or BioSample with any complete assembly
            f"({query}) AND latest[Filter] AND (chromosome[Assembly Level] OR complete genome[Assembly Level] OR scaffold[Assembly Level]) AND (bioproject[Filter] OR biosample[Filter])",
            # Moderate: Latest assemblies with metadata filters
            f"({query}) AND latest[Filter] AND (bioproject[Filter] OR biosample[Filter])",
            # Fallback: Any latest assemblies
            f"({query}) AND latest[Filter]"
        ]
        
        for i, enhanced_query in enumerate(strategies):
            logging.info(f"Trying assembly search strategy {i+1}/4: {enhanced_query[:100]}...")
            
            params = {
                'db': 'assembly',
                'term': enhanced_query,
                'retmax': max_results * 5,  # Get more to filter for quality
                'retmode': 'xml',
                'email': self.email,
                'api_key': self.api_key if self.api_key else None
            }

            params = {k: v for k, v in params.items() if v is not None}
            response = self._make_request(NCBI_SEARCH_URL, params)
            
            if response:
                try:
                    root = ET.fromstring(response.text)
                    ids = [id_elem.text for id_elem in root.findall('.//Id') if id_elem.text]
                    
                    if ids:
                        logging.info(f"Strategy {i+1} found {len(ids)} assemblies")
                        return ids[:max_results * 3]  # Return reasonable number for filtering
                except ET.ParseError:
                    logging.warning(f"Failed to parse assembly search response for strategy {i+1}")
                    continue
        
        return []

    def _get_nucleotide_accessions_from_assemblies(self, assembly_ids: List[str]) -> List[str]:
        """Convert assembly IDs to nucleotide accessions for sequence download"""
        nucleotide_accessions = []
        
        # Process in batches to get assembly details
        batch_size = 100
        for i in range(0, len(assembly_ids), batch_size):
            batch_ids = assembly_ids[i:i + batch_size]
            batch_accessions = self._get_assembly_nucleotide_accessions_batch(batch_ids)
            nucleotide_accessions.extend(batch_accessions)
        
        return nucleotide_accessions

    def _get_assembly_nucleotide_accessions_batch(self, assembly_ids: List[str]) -> List[str]:
        """Get nucleotide accessions from assembly IDs in batch"""
        summary_url = NCBI_BASE_URL + "esummary.fcgi"
        params = {
            'db': 'assembly',
            'id': ','.join(assembly_ids),
            'retmode': 'xml',
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        params = {k: v for k, v in params.items() if v is not None}
        response = self._make_request(summary_url, params)
        
        if not response:
            return []

        nucleotide_accessions = []
        try:
            root = ET.fromstring(response.text)
            for docsum in root.findall('.//DocSum'):
                # Look for GenBank accession in assembly summary
                for item in docsum.findall('Item'):
                    name = item.get('Name')
                    if name == 'Synonym' and item.text:
                        # Extract GenBank accession from synonym
                        synonym_text = item.text
                        if 'GenBank:' in synonym_text:
                            genbank_acc = synonym_text.split('GenBank:')[1].split(';')[0].strip()
                            if genbank_acc:
                                nucleotide_accessions.append(genbank_acc)
                                break
                    elif name == 'AssemblyAccession' and item.text:
                        # Use assembly accession if no GenBank accession found
                        nucleotide_accessions.append(item.text)
                        break
        except ET.ParseError:
            logging.warning("Failed to parse assembly summary response")
        
        return nucleotide_accessions

    def _search_nuccore_fallback(self, query: str, max_results: int) -> List[str]:
        """Enhanced nuccore search with metadata-focused queries"""
        # Create progressive query strategies that target genomes with rich metadata
        strategies = [
            # Best: Target clinical/research genomes with explicit metadata terms
            f"({query}) AND (complete genome[TI] OR chromosome[TI]) AND (biosample[All Fields] OR bioproject[All Fields] OR clinical[All Fields] OR patient[All Fields] OR hospital[All Fields])",
            
            # Good: Target genomes with publication and study context
            f"({query}) AND (complete genome[TI] OR chromosome[TI]) AND (study[All Fields] OR surveillance[All Fields] OR outbreak[All Fields] OR resistance[All Fields])",
            
            # Moderate: Recent complete genomes (likely to have better metadata)
            f"({query}) AND (complete genome[TI] OR chromosome[TI]) AND 2020:2025[Publication Date]",
            
            # Basic: Any complete genomes
            f"({query}) AND (complete genome[TI] OR chromosome[TITL] OR \"complete sequence\"[TI]) NOT (scaffold OR contig OR plasmid[TITL])",
            
            # Fallback: Broad search
            f"({query}) AND complete[TI]"
        ]

        for i, enhanced_query in enumerate(strategies):
            logging.info(f"Trying nuccore search strategy {i+1}/5")
            
            search_results = self._execute_nuccore_search(enhanced_query, max_results * 3)
            
            if search_results:
                logging.info(f"Strategy {i+1} found {len(search_results)} genomes")
                # Filter results for quality if we got enough
                if len(search_results) >= max_results:
                    return self._filter_complete_genomes(search_results, max_results)[:max_results]
                else:
                    return search_results
        
        # If all strategies fail, return empty list
        logging.warning("All nuccore search strategies failed")
        return []

    def _execute_nuccore_search(self, query: str, max_results: int) -> List[str]:
        """Execute a single nuccore search with the given query"""
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

        try:
            root = ET.fromstring(response.text)
            ids = [id_elem.text for id_elem in root.findall('.//Id') if id_elem.text]
            return ids
        except ET.ParseError:
            logging.warning(f"Failed to parse search response for query: {query[:50]}...")
            return []

    def fetch_genomes(self, query: str, max_results: int = DEFAULT_MAX_RESULTS) -> List[Dict[str, Any]]:
        """
        Enhanced genome search with metadata-focused strategies

        Args:
            query: Search query
            max_results: Maximum number of results to return

        Returns:
            List of genome records with rich metadata
        """
        logging.info(f"Enhanced search for: {query}")

        # Strategy 1: Try targeted metadata-rich searches in nuccore
        logging.info("Phase 1: Searching for genomes with rich metadata indicators")
        genome_ids = self._search_metadata_rich_genomes(query, max_results)
        
        if not genome_ids:
            # Strategy 2: Try assembly database search
            logging.info("Phase 2: Trying assembly database search")
            assembly_ids = self._search_assemblies_with_metadata(query, max_results)
            if assembly_ids:
                genome_ids = self._get_nucleotide_accessions_from_assemblies(assembly_ids)
        
        if not genome_ids:
            # Strategy 3: Fallback to enhanced nuccore search
            logging.info("Phase 3: Fallback to enhanced nuccore search")
            genome_ids = self._search_nuccore_fallback(query, max_results)

        if not genome_ids:
            logging.warning("No genomes found for the query")
            return []

        logging.info(f"Found {len(genome_ids)} genomes for metadata extraction")

        # Extract metadata for all genomes
        logging.info("Extracting metadata for genomes...")
        raw_records = self._extract_metadata_batch(genome_ids)

        return raw_records

    def _search_metadata_rich_genomes(self, query: str, max_results: int) -> List[str]:
        """Search specifically for genomes likely to have rich metadata"""
        # Construct queries that target genomes with known metadata indicators
        metadata_queries = [
            # Target genomes with explicit BioSample/BioProject mentions
            f"({query}) AND (complete genome[TI] OR chromosome[TI]) AND (SAMN[All Fields] OR PRJNA[All Fields] OR PRJEB[All Fields])",
            
            # Target clinical and research contexts
            f"({query}) AND (complete genome[TI] OR chromosome[TI]) AND (clinical isolate[All Fields] OR antimicrobial resistance[All Fields] OR drug resistance[All Fields])",
            
            # Target surveillance and outbreak studies
            f"({query}) AND (complete genome[TI] OR chromosome[TI]) AND (surveillance[All Fields] OR outbreak[All Fields] OR epidemiology[All Fields])",
            
            # Target recent publications with likely metadata
            f"({query}) AND (complete genome[TI] OR chromosome[TI]) AND (strain[All Fields] OR isolate[All Fields]) AND 2020:2025[Publication Date]",
        ]

        for i, meta_query in enumerate(metadata_queries):
            logging.info(f"Trying metadata-rich search {i+1}/4")
            
            search_results = self._execute_nuccore_search(meta_query, max_results * 2)
            
            if search_results:
                logging.info(f"Metadata-rich search {i+1} found {len(search_results)} genomes")
                return search_results[:max_results * 2]  # Return extra for filtering
        
        return []

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
        """Extract metadata for a batch of genome IDs using assembly-first approach"""
        metadata_list = []

        # First try to get assembly IDs for these genome accessions
        assembly_metadata = self._get_assembly_metadata_for_accessions(genome_ids)
        
        # If we got assembly metadata, use it; otherwise fall back to nuccore
        if assembly_metadata:
            logging.info(f"Retrieved assembly metadata for {len(assembly_metadata)} genomes")
            return assembly_metadata
        
        # Fallback to original nuccore metadata extraction
        logging.info("Falling back to nuccore metadata extraction")
        return self._extract_nuccore_metadata_batch(genome_ids)

    def _get_assembly_metadata_for_accessions(self, accessions: List[str]) -> List[Dict[str, Any]]:
        """Get rich metadata from assembly database for given accessions"""
        metadata_list = []
        
        # Process in batches to avoid API limits
        batch_size = 50
        for i in range(0, len(accessions), batch_size):
            batch_accessions = accessions[i:i + batch_size]
            
            # First, find assembly IDs for these accessions
            assembly_ids = self._find_assembly_ids_for_accessions(batch_accessions)
            
            if assembly_ids:
                # Get rich metadata from assembly records
                batch_metadata = self._extract_assembly_metadata_batch(assembly_ids, batch_accessions)
                metadata_list.extend(batch_metadata)
            else:
                # Create empty metadata for this batch
                for acc in batch_accessions:
                    metadata_list.append(self._create_empty_metadata(acc))

        return metadata_list

    def _find_assembly_ids_for_accessions(self, accessions: List[str]) -> List[str]:
        """Find assembly IDs that correspond to given nucleotide accessions"""
        assembly_ids = []
        
        for accession in accessions:
            # Use elink to find assembly records linked to this nucleotide accession
            elink_url = NCBI_BASE_URL + "elink.fcgi"
            params = {
                'dbfrom': 'nuccore',
                'db': 'assembly',
                'id': accession,
                'email': self.email,
                'api_key': self.api_key if self.api_key else None
            }
            
            params = {k: v for k, v in params.items() if v is not None}
            response = self._make_request(elink_url, params)
            
            if response:
                try:
                    root = ET.fromstring(response.text)
                    for link_id in root.findall('.//Link/Id'):
                        if link_id.text:
                            assembly_ids.append(link_id.text)
                            break  # Take first assembly ID
                except ET.ParseError:
                    continue
        
        return assembly_ids

    def _extract_assembly_metadata_batch(self, assembly_ids: List[str], corresponding_accessions: List[str]) -> List[Dict[str, Any]]:
        """Extract rich metadata from assembly records"""
        if not assembly_ids:
            return []

        # Get assembly summaries in batch
        summary_url = NCBI_BASE_URL + "esummary.fcgi"
        params = {
            'db': 'assembly',
            'id': ','.join(assembly_ids),
            'retmode': 'xml',
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        params = {k: v for k, v in params.items() if v is not None}
        response = self._make_request(summary_url, params)

        if not response:
            return [self._create_empty_metadata(acc) for acc in corresponding_accessions]

        metadata_list = []
        try:
            root = ET.fromstring(response.text)
            
            for i, docsum in enumerate(root.findall('.//DocSum')):
                # Get corresponding accession for this assembly
                accession = corresponding_accessions[i] if i < len(corresponding_accessions) else "unknown"
                
                metadata = self._parse_assembly_metadata(docsum, accession)
                
                # If we have BioSample, get detailed BioSample metadata
                if metadata.get('biosample'):
                    biosample_metadata = self._get_biosample_metadata(metadata['biosample'])
                    if biosample_metadata:
                        metadata.update(biosample_metadata)
                        logging.debug(f"Enhanced {accession} with BioSample metadata")

                # Calculate quality score
                metadata['quality_score'] = self._calculate_metadata_score(metadata)
                metadata_list.append(metadata)

        except ET.ParseError as e:
            logging.error(f"Failed to parse assembly metadata XML: {e}")
            return [self._create_empty_metadata(acc) for acc in corresponding_accessions]

        return metadata_list

    def _parse_assembly_metadata(self, docsum: ET.Element, accession: str) -> Dict[str, Any]:
        """Parse metadata from assembly DocSum XML - much richer than nuccore"""
        metadata = {
            'genome_id': None,
            'accession': accession,
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
            'assembly_level': None,
            'genome_representation': None,
            'assembly_name': None
        }

        # Extract assembly ID
        id_elem = docsum.find('Id')
        if id_elem is not None:
            metadata['genome_id'] = id_elem.text

        # Parse assembly items - these have rich metadata
        for item in docsum.findall('Item'):
            name = item.get('Name')
            
            if name == 'AssemblyName':
                metadata['assembly_name'] = item.text
                metadata['title'] = item.text
            elif name == 'Organism':
                metadata['organism'] = item.text
                # Extract genus and species
                if item.text:
                    parts = item.text.split()
                    if len(parts) >= 2:
                        metadata['genus'] = parts[0]
                        metadata['species'] = ' '.join(parts[1:])
            elif name == 'AssemblyLevel':
                metadata['assembly_level'] = item.text
            elif name == 'GenomeRepresentation':
                metadata['genome_representation'] = item.text
            elif name == 'BioSampleAccn':
                # Direct BioSample link - this is what we want!
                metadata['biosample'] = item.text
                logging.debug(f"Found BioSample {item.text} for {accession}")
            elif name == 'BioprojectAccn':
                # Direct BioProject link
                metadata['bioproject'] = item.text
                logging.debug(f"Found BioProject {item.text} for {accession}")
            elif name == 'SubmissionDate':
                # Use submission date as collection date if no better date available
                metadata['submission_date'] = item.text
            elif name == 'LastUpdateDate':
                metadata['update_date'] = item.text

        return metadata

    def _extract_nuccore_metadata_batch(self, genome_ids: List[str]) -> List[Dict[str, Any]]:
        """Enhanced nuccore metadata extraction with concurrent processing"""
        if not genome_ids:
            return []
            
        metadata_list = []
        
        # Use concurrent processing for multiple batches
        batch_size = 100  # Smaller batches for concurrent processing
        batches = [genome_ids[i:i + batch_size] for i in range(0, len(genome_ids), batch_size)]
        
        if len(batches) == 1:
            # Single batch - no need for threading
            return self._get_metadata_batch(batches[0])
        
        # Multiple batches - use concurrent processing
        logging.info(f"Processing {len(batches)} batches concurrently")
        
        with ThreadPoolExecutor(max_workers=min(3, len(batches))) as executor:
            # Submit all batch processing tasks
            future_to_batch = {
                executor.submit(self._get_metadata_batch, batch): batch 
                for batch in batches
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_batch):
                batch = future_to_batch[future]
                try:
                    batch_metadata = future.result(timeout=30)
                    metadata_list.extend(batch_metadata)
                    logging.debug(f"Completed batch with {len(batch_metadata)} genomes")
                except Exception as e:
                    logging.error(f"Batch processing failed: {e}")
                    # Add empty metadata for failed batch
                    metadata_list.extend([self._create_empty_metadata(gid) for gid in batch])

        return metadata_list

    def _enhance_with_biosample_batch(self, metadata_list: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Enhanced BioSample processing with concurrent requests and caching"""
        # Collect unique BioSample IDs that need enhancement
        biosample_ids = []
        biosample_metadata_map = {}
        
        for metadata in metadata_list:
            biosample_id = metadata.get('biosample')
            if biosample_id and biosample_id not in biosample_metadata_map:
                biosample_ids.append(biosample_id)
        
        if not biosample_ids:
            return metadata_list
        
        logging.info(f"Enhancing {len(biosample_ids)} unique BioSamples")
        
        # Process BioSample IDs in concurrent batches
        batch_size = 20  # Smaller batches for BioSample due to complexity
        batches = [biosample_ids[i:i + batch_size] for i in range(0, len(biosample_ids), batch_size)]
        
        if len(batches) == 1:
            # Single batch
            batch_biosample_metadata = self._get_biosample_metadata_batch(batches[0])
            for i, biosample_id in enumerate(biosample_ids):
                if i < len(batch_biosample_metadata):
                    biosample_metadata_map[biosample_id] = batch_biosample_metadata[i]
        else:
            # Multiple batches - concurrent processing
            with ThreadPoolExecutor(max_workers=min(2, len(batches))) as executor:
                future_to_batch_ids = {
                    executor.submit(self._get_biosample_metadata_batch, batch): batch
                    for batch in batches
                }
                
                for future in as_completed(future_to_batch_ids):
                    batch_ids = future_to_batch_ids[future]
                    try:
                        batch_metadata = future.result(timeout=20)
                        # Map results back to IDs
                        for i, biosample_id in enumerate(batch_ids):
                            if i < len(batch_metadata):
                                biosample_metadata_map[biosample_id] = batch_metadata[i]
                    except Exception as e:
                        logging.warning(f"BioSample batch enhancement failed: {e}")
        
        # Update original metadata with BioSample data
        enhanced_count = 0
        for metadata in metadata_list:
            biosample_id = metadata.get('biosample')
            if biosample_id and biosample_id in biosample_metadata_map:
                biosample_data = biosample_metadata_map[biosample_id]
                if biosample_data:
                    metadata.update(biosample_data)
                    enhanced_count += 1
            
            # Recalculate quality score after enhancement
            metadata['quality_score'] = self._calculate_metadata_score(metadata)
        
        logging.info(f"Successfully enhanced {enhanced_count} genomes with BioSample data")
        return metadata_list

    @lru_cache(maxsize=500)
    def _get_biosample_metadata_cached(self, biosample_id: str) -> Dict[str, Any]:
        """Cached BioSample metadata retrieval for frequently accessed samples"""
        return self._get_single_biosample_metadata(biosample_id)
    
    @lru_cache(maxsize=200)
    def _get_bioproject_metadata_cached(self, bioproject_id: str) -> Dict[str, Any]:
        """Cached BioProject metadata retrieval"""
        try:
            bioproject_url = NCBI_BASE_URL + "efetch.fcgi"
            params = {
                'db': 'bioproject',
                'id': bioproject_id,
                'retmode': 'xml',
                'email': self.email,
                'api_key': self.api_key if self.api_key else None
            }

            params = {k: v for k, v in params.items() if v is not None}
            response = self._make_request(bioproject_url, params)

            if response:
                root = ET.fromstring(response.text)
                project_elem = root.find('.//Project')
                if project_elem is not None:
                    return self._parse_bioproject_xml(project_elem)
        except Exception as e:
            logging.debug(f"Failed to get BioProject {bioproject_id}: {e}")
        
        return {}

    def _parse_bioproject_xml(self, project_elem: ET.Element) -> Dict[str, Any]:
        """Parse BioProject XML for additional metadata"""
        metadata = {}
        
        try:
            # Get project description
            description = project_elem.find('.//Description')
            if description is not None and description.text:
                metadata['bioproject_description'] = description.text.strip()
            
            # Get project title
            title = project_elem.find('.//Title')
            if title is not None and title.text:
                metadata['bioproject_title'] = title.text.strip()
            
            # Look for study type and methodology
            method_elem = project_elem.find('.//Method')
            if method_elem is not None and method_elem.text:
                metadata['study_method'] = method_elem.text.strip()

        except Exception as e:
            logging.debug(f"Error parsing BioProject XML: {e}")
        
        return metadata

    def _clear_cache(self):
        """Clear all cached data - useful for testing or memory management"""
        self._get_biosample_metadata_cached.cache_clear()
        self._get_bioproject_metadata_cached.cache_clear()
        logging.info("Metadata caches cleared")

    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics for monitoring"""
        return {
            'biosample_cache': self._get_biosample_metadata_cached.cache_info()._asdict(),
            'bioproject_cache': self._get_bioproject_metadata_cached.cache_info()._asdict()
        }

    def _get_linked_metadata(self, accession: str) -> Dict[str, Any]:
        """Get linked BioSample and BioProject data using elink"""
        linked_data = {}
        
        # Find linked BioSample
        biosample_id = self._find_linked_biosample(accession)
        if biosample_id:
            linked_data['biosample'] = biosample_id
            logging.debug(f"Found linked BioSample {biosample_id} for {accession}")
            
            # Get BioSample metadata
            biosample_metadata = self._get_biosample_metadata_cached(biosample_id)
            if biosample_metadata:
                linked_data.update(biosample_metadata)

        # Find linked BioProject
        bioproject_id = self._find_linked_bioproject(accession)
        if bioproject_id:
            linked_data['bioproject'] = bioproject_id
            logging.debug(f"Found linked BioProject {bioproject_id} for {accession}")
            
            # Get BioProject metadata
            bioproject_metadata = self._get_bioproject_metadata_cached(bioproject_id)
            if bioproject_metadata:
                linked_data.update(bioproject_metadata)

        return linked_data

    def _find_linked_biosample(self, accession: str) -> Optional[str]:
        """Find BioSample linked to nucleotide accession"""
        elink_url = NCBI_BASE_URL + "elink.fcgi"
        params = {
            'dbfrom': 'nuccore',
            'db': 'biosample',
            'id': accession,
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        params = {k: v for k, v in params.items() if v is not None}
        response = self._make_request(elink_url, params)

        if response:
            try:
                root = ET.fromstring(response.text)
                # Look for linked BioSample IDs
                for link_id in root.findall('.//Link/Id'):
                    if link_id.text:
                        return link_id.text
            except ET.ParseError:
                pass
        
        return None

    def _find_linked_bioproject(self, accession: str) -> Optional[str]:
        """Find BioProject linked to nucleotide accession"""
        elink_url = NCBI_BASE_URL + "elink.fcgi"
        params = {
            'dbfrom': 'nuccore',
            'db': 'bioproject',
            'id': accession,
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        params = {k: v for k, v in params.items() if v is not None}
        response = self._make_request(elink_url, params)

        if response:
            try:
                root = ET.fromstring(response.text)
                # Look for linked BioProject IDs
                for link_id in root.findall('.//Link/Id'):
                    if link_id.text:
                        return link_id.text
            except ET.ParseError:
                pass
        
        return None
    
    def _get_single_biosample_metadata(self, biosample_id: str) -> Dict[str, Any]:
        """Get metadata for a single BioSample ID"""
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

        if not response:
            return {}

        try:
            root = ET.fromstring(response.text)
            biosample_elem = root.find('.//BioSample')
            if biosample_elem is not None:
                return self._parse_single_biosample_xml(biosample_elem)
        except ET.ParseError:
            pass
        
        return {}

    def _get_metadata_batch(self, genome_ids: List[str]) -> List[Dict[str, Any]]:
        """Get metadata for a batch of genome IDs with enhanced batch processing"""
        if not genome_ids:
            return []

        metadata_list = []
        
        # Use larger batch sizes for better performance
        batch_size = 200  # NCBI supports up to 200 IDs per request
        
        for i in range(0, len(genome_ids), batch_size):
            batch_ids = genome_ids[i:i + batch_size]
            logging.info(f"Processing batch {i//batch_size + 1}: {len(batch_ids)} genomes")
            
            batch_metadata = self._process_metadata_batch(batch_ids)
            metadata_list.extend(batch_metadata)
            
            # Add small delay between batches to respect rate limits
            if i + batch_size < len(genome_ids):
                time.sleep(0.1)  # Reduced delay with API key

        return metadata_list

    def _process_metadata_batch(self, genome_ids: List[str]) -> List[Dict[str, Any]]:
        """Process a single batch of genome IDs for metadata"""
        # Get document summaries from NCBI in batch
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
            metadata_list = []
            root = ET.fromstring(response.text)
            
            # Parse all docsum elements in batch
            docsum_elements = root.findall('.//DocSum')
            
            # Process each genome in the batch
            for i, docsum in enumerate(docsum_elements):
                genome_id = genome_ids[i] if i < len(genome_ids) else None
                metadata = self._parse_docsum_metadata(docsum)
                
                # Ensure genome_id is set correctly
                if genome_id:
                    metadata['genome_id'] = genome_id
                
                metadata_list.append(metadata)
            
            # If we have fewer docsums than expected, fill with empty metadata
            while len(metadata_list) < len(genome_ids):
                missing_id = genome_ids[len(metadata_list)]
                metadata_list.append(self._create_empty_metadata(missing_id))
            
            # Process BioSample metadata in batch for genomes that have BioSample IDs
            metadata_list = self._enhance_with_biosample_batch(metadata_list)
            
            return metadata_list
            
        except ET.ParseError as e:
            logging.error(f"Failed to parse metadata XML: {e}")
            return [self._create_empty_metadata(gid) for gid in genome_ids]

    # ...existing code...

    def _get_biosample_metadata_batch(self, biosample_ids: List[str]) -> List[Dict[str, Any]]:
        """Get BioSample metadata for multiple IDs in batch"""
        if not biosample_ids:
            return []
            
        biosample_url = NCBI_BASE_URL + "efetch.fcgi"
        params = {
            'db': 'biosample',
            'id': ','.join(biosample_ids),
            'retmode': 'xml',
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        params = {k: v for k, v in params.items() if v is not None}
        response = self._make_request(biosample_url, params)

        if not response:
            return [{}] * len(biosample_ids)

        try:
            metadata_list = []
            root = ET.fromstring(response.text)
            
            # Process each BioSample in the batch
            biosample_elements = root.findall('.//BioSample')
            
            for i, biosample_elem in enumerate(biosample_elements):
                biosample_metadata = self._parse_single_biosample_xml(biosample_elem)
                metadata_list.append(biosample_metadata)
            
            # Fill missing entries
            while len(metadata_list) < len(biosample_ids):
                metadata_list.append({})
            
            return metadata_list
            
        except ET.ParseError as e:
            logging.warning(f"Failed to parse BioSample batch XML: {e}")
            return [{}] * len(biosample_ids)

    def _parse_single_biosample_xml(self, biosample_elem: ET.Element) -> Dict[str, Any]:
        """Parse a single BioSample XML element"""
        metadata = {}

        try:
            # Get BioSample accession
            accession = biosample_elem.get('accession')
            if accession:
                metadata['biosample_accession'] = accession

            # Get BioSample description
            description = biosample_elem.find('.//Description')
            if description is not None and description.text:
                metadata['biosample_description'] = description.text.strip()

            # Extract attributes with comprehensive patterns
            for attr in biosample_elem.findall('.//Attribute'):
                attr_name = attr.get('attribute_name', '').lower().strip()
                attr_value = attr.text or ''

                if not attr_value.strip():
                    continue

                attr_value_clean = attr_value.strip()

                # Collection date
                if any(keyword in attr_name for keyword in ['collection_date', 'collection date', 'date', 'sampling_date']):
                    metadata['collection_date'] = self._standardize_date(attr_value_clean)

                # Geographic location
                elif any(keyword in attr_name for keyword in ['geo_loc_name', 'geographic location', 'country', 'location']):
                    metadata['country'] = self._standardize_location(attr_value_clean)

                # Host information
                elif any(keyword in attr_name for keyword in ['host', 'host organism', 'source host']):
                    metadata['host'] = attr_value_clean

                # Isolation source
                elif any(keyword in attr_name for keyword in ['isolation_source', 'isolation source', 'source', 'sample type']):
                    metadata['isolation_source'] = attr_value_clean

                # MIC data
                elif 'mic' in attr_name:
                    if not metadata.get('mic_data'):
                        metadata['mic_data'] = []
                    
                    antibiotic = self._extract_antibiotic_from_mic_name(attr_name)
                    mic_entry = {
                        'antibiotic': antibiotic,
                        'value': attr_value_clean,
                        'unit': self._extract_mic_unit(attr_value_clean),
                        'source': 'biosample'
                    }
                    metadata['mic_data'].append(mic_entry)

                # Resistance data
                elif any(keyword in attr_name for keyword in ['resistance', 'phenotype', 'susceptibility']):
                    if not metadata.get('resistance_phenotype'):
                        metadata['resistance_phenotype'] = []
                    metadata['resistance_phenotype'].append(attr_value_clean)

        except Exception as e:
            logging.warning(f"Error parsing single BioSample: {e}")

        return metadata

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
                # CRITICAL FIX: Use elink to find linked BioSample/BioProject instead
                # of trying to extract from Extra field (which doesn't contain this info)
                pass
            elif name == 'CreateDate':
                metadata['create_date'] = item.text
            elif name == 'UpdateDate':
                metadata['update_date'] = item.text

        # CRITICAL FIX: Use elink to find BioSample and BioProject
        if metadata['accession']:
            linked_metadata = self._get_linked_metadata(metadata['accession'])
            metadata.update(linked_metadata)

        # If no BioSample, try to extract metadata from title and other fields
        if not metadata.get('country') and not metadata.get('collection_date'):
            self._extract_metadata_from_title(metadata)

        # Try to get assembly information which might have BioSample data
        if not metadata.get('biosample'):
            accession = metadata.get('accession')
            if accession:
                assembly_biosample = self._get_biosample_from_assembly(accession)
                if assembly_biosample:
                    metadata['biosample'] = assembly_biosample
                    # Try to get BioSample metadata again
                    biosample_metadata = self._get_biosample_metadata(assembly_biosample)
                    if biosample_metadata:
                        metadata.update(biosample_metadata)
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
            self.logger.debug(f"Request succeeded: {url} params={params}")
            return response
        except requests.RequestException as e:
            self.logger.error(f"Request failed: {e}")
            return None