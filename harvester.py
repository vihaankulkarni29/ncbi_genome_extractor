#!/usr/bin/env python3
"""
Federated Genome Harvester - Multi-source genome data retrieval tool
"""

import argparse
import json
import logging
import os
import sys
import time
from typing import List, Dict, Any
from concurrent.futures import ThreadPoolExecutor, as_completed

from clients.ncbi_client import NCBIClient
from clients.bvbrc_client import BVBRCClient
from harmonizer import harmonize_data
from config import DEFAULT_OUTPUT_DIR


def download_genome_parallel(accession: str, database: str, output_dir: str, sources: List[str]) -> tuple:
    """Download a single genome - designed for parallel execution"""
    try:
        if database.lower() == 'ncbi' and 'ncbi' in sources:
            client = NCBIClient()
            result = client.download_fasta(accession, output_dir)
            if isinstance(result, tuple) and len(result) == 2:
                success, final_accession = result
                if success:
                    return True, final_accession, None
                else:
                    return False, accession, "Download failed"
            else:
                if result:
                    return True, accession, None
                else:
                    return False, accession, "Download failed"
        
        elif database.lower() == 'bvbrc' and 'bvbrc' in sources:
            client = BVBRCClient()
            result = client.download_fasta(accession, output_dir)
            if isinstance(result, tuple) and len(result) == 2:
                success, final_accession = result
                if success:
                    return True, final_accession, None
                else:
                    return False, accession, "Download failed"
            else:
                if result:
                    return True, accession, None
                else:
                    return False, accession, "Download failed"
        
        elif database.lower() == 'enterobase' and 'enterobase' in sources:
            from clients.enterobase_client import EnteroBaseClient
            client = EnteroBaseClient()
            result = client.download_fasta(accession, output_dir)
            if isinstance(result, tuple) and len(result) == 2:
                success, final_accession = result
                if success:
                    return True, final_accession, None
                else:
                    return False, accession, "Download failed"
            else:
                if result:
                    return True, accession, None
                else:
                    return False, accession, "Download failed"
        
        elif database.lower() == 'patric' and 'patric' in sources:
            from clients.patric_client import PATRICClient
            client = PATRICClient()
            result = client.download_fasta(accession, output_dir)
            if isinstance(result, tuple) and len(result) == 2:
                success, final_accession = result
                if success:
                    return True, final_accession, None
                else:
                    return False, accession, "Download failed"
            else:
                if result:
                    return True, accession, None
                else:
                    return False, accession, "Download failed"
        
        else:
            return False, accession, f"Database {database} not in sources {sources}"
            
    except Exception as e:
        return False, accession, str(e)


def download_genomes_parallel(records: List[Dict[str, Any]], output_dir: str, sources: List[str], max_workers: int = 4) -> List[str]:
    """Download multiple genomes in parallel"""
    logging.info(f"Starting parallel download of {len(records)} genomes with {max_workers} workers...")
    
    successful = []
    failed = []
    error_details = []
    start_time = time.time()
    
    # Prepare download tasks
    download_tasks = []
    for record in records:
        accession = record.get('accession', '')
        database = record.get('database', '').lower()
        if accession and database:
            download_tasks.append((accession, database, output_dir, sources))
    
    if not download_tasks:
        logging.warning("No valid accessions found for download")
        return successful
    
    # Execute downloads in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all download tasks
        future_to_accession = {
            executor.submit(download_genome_parallel, accession, database, output_dir, sources): accession
            for accession, database, output_dir, sources in download_tasks
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_accession):
            accession = future_to_accession[future]
            try:
                success, final_accession, error = future.result(timeout=300)  # 5 minute timeout per download
                if success:
                    successful.append(final_accession)
                    logging.debug(f"✓ Downloaded: {final_accession}")
                else:
                    failed.append(accession)
                    error_details.append({
                        'accession': accession,
                        'error': error,
                        'source': 'download',
                        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                    })
                    logging.warning(f"✗ Failed: {accession} - {error}")
            except Exception as e:
                failed.append(accession)
                error_details.append({
                    'accession': accession,
                    'error': str(e),
                    'source': 'exception',
                    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                })
                logging.error(f"✗ Exception downloading {accession}: {e}")
    
    end_time = time.time()
    total_time = end_time - start_time

    logging.info(f"Parallel download completed in {total_time:.2f} seconds")
    logging.info(f"✓ Successfully downloaded: {len(successful)}/{len(download_tasks)} genomes")
    logging.info(f"✗ Failed downloads: {len(failed)}/{len(download_tasks)} genomes")
    if successful:
        logging.info(f"⚡ Average speed: {total_time/len(successful):.2f} seconds per genome")

    # Save error report if there are failures
    if error_details:
        error_report_path = os.path.join(output_dir, 'failed_downloads_report.csv')
        save_error_report_csv(error_details, error_report_path)
        logging.info(f"Detailed error report saved to {error_report_path}")

    return successful
def save_error_report_csv(error_details: list, output_file: str):
    """Save error details to a CSV file for post-run analysis"""
    import csv
    if not error_details:
        return
    fieldnames = ['accession', 'error', 'source', 'timestamp']
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(error_details)


def main():
    """Main CLI function"""
    parser = argparse.ArgumentParser(
        description="Federated Genome Harvester - Multi-database genome data retrieval with quality filtering",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # GENERAL GENOME SEARCH - Any organism, any region
  python harvester.py --query "Escherichia coli" --max_results 50 --min_quality_score 3

  # STRICT QUALITY FILTERING - Only high-quality genomes with AMR data
  python harvester.py --query "Klebsiella pneumoniae AND carbapenem" --require_amr_data --require_location --min_quality_score 5

  # Search specific AMR-focused databases
  python harvester.py --source enterobase --query "Salmonella enterica" --max_results 25
  python harvester.py --source patric --query "Acinetobacter baumannii" --max_results 25

  # Search traditional databases
  python harvester.py --source ncbi --query "SARS-CoV-2[Organism]" --max_results 10
  python harvester.py --source bvbrc --query "Pseudomonas aeruginosa" --max_results 10

  # Download specific genomes from any source
  python harvester.py --query "CP186630,NZ_CP064825" --download_only

  # Geographic research - any region
  python harvester.py --query "Escherichia coli AND Europe" --max_results 100 --download --metadata_format csv --min_quality_score 4

  # AMR surveillance - any pathogen
  python harvester.py --query "Staphylococcus aureus AND MRSA" --require_amr_data --max_results 50
        """
    )

    parser.add_argument('--source', choices=['ncbi', 'bvbrc', 'enterobase', 'patric', 'all'],
                       default='all',
                       help='Data source to query (ncbi, bvbrc, enterobase, patric, or all for all sources, default: all)')
    parser.add_argument('--query', help='Search query for the selected source')
    parser.add_argument('--max_results', type=int, default=100,
                       help='Maximum number of results to return (default: 100)')
    parser.add_argument('--output_dir', default=DEFAULT_OUTPUT_DIR,
                       help=f'Output directory for files (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--download', action='store_true',
                       help='Download FASTA sequences for found genomes')
    parser.add_argument('--download_only', action='store_true',
                       help='Only download FASTA sequences (skip metadata extraction)')
    parser.add_argument('--parallel_downloads', type=int, default=4,
                       help='Number of parallel download workers (default: 4, max: 8)')
    parser.add_argument('--metadata_format', choices=['json', 'csv'], default='json',
                        help='Format for metadata output (default: json)')
    parser.add_argument('--min_quality_score', type=int, default=0,
                        help='Minimum quality score for genome inclusion (0-10, default: 0)')
    parser.add_argument('--require_amr_data', action='store_true',
                        help='Only download genomes with AMR resistance data')
    parser.add_argument('--require_location', action='store_true',
                        help='Only download genomes with geographic location data')
    parser.add_argument('--log_level', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level')

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Validate arguments
    if not args.query and not args.download_only:
        parser.error("Either --query or --download_only must be provided")

    # Validate parallel downloads argument
    if args.parallel_downloads < 1 or args.parallel_downloads > 8:
        parser.error("--parallel_downloads must be between 1 and 8")

    # Determine which sources to search
    if args.source == 'all':
        sources = ['ncbi', 'bvbrc', 'enterobase', 'patric']
    else:
        sources = [args.source]

    all_harmonized_records = []

    try:
        if args.download_only:
            # Download mode - requires specific genome IDs/accessions
            if not args.query:
                parser.error("--query must be provided with --download_only (genome IDs/accessions)")

            genome_ids = args.query.split(',')
            logging.info(f"Downloading {len(genome_ids)} genomes from {args.source.upper()} using {args.parallel_downloads} parallel workers...")

            # Convert genome IDs to records format for parallel download
            download_records = []
            for genome_id in genome_ids:
                genome_id = genome_id.strip()
                # For download-only mode, we'll try the first available source
                for source in sources:
                    download_records.append({
                        'accession': genome_id,
                        'database': source.upper()
                    })
                    break  # Only try the first source for each genome

            # Use parallel download
            successful = download_genomes_parallel(download_records, args.output_dir, sources, args.parallel_downloads)
            logging.info(f"Parallel download completed: {len(successful)} genomes successfully downloaded")

        else:
            # Search and retrieve mode
            for source in sources:
                logging.info(f"Searching {source.upper()} with query: {args.query}")

                # Initialize the appropriate client
                client = None
                if source == 'ncbi':
                    client = NCBIClient()
                elif source == 'bvbrc':
                    client = BVBRCClient()
                elif source == 'enterobase':
                    from clients.enterobase_client import EnteroBaseClient
                    client = EnteroBaseClient()
                elif source == 'patric':
                    from clients.patric_client import PATRICClient
                    client = PATRICClient()

                if client is None:
                    logging.error(f"Unknown source: {source}")
                    continue

                # Fetch raw genome records
                raw_records = client.fetch_genomes(args.query, args.max_results)

                if not raw_records:
                    logging.warning(f"No genomes found in {source.upper()} for the query")
                    continue

                # Harmonize the data
                harmonized_records = harmonize_data(raw_records, source)
                all_harmonized_records.extend(harmonized_records)
                logging.info(f"Found {len(harmonized_records)} genomes in {source.upper()}")

            if not all_harmonized_records:
                logging.warning("No genomes found in any database for the query")
                return

            logging.info(f"Total harmonized records: {len(all_harmonized_records)}")

            # Save harmonized metadata
            os.makedirs(args.output_dir, exist_ok=True)
            source_name = 'all' if args.source == 'all' else args.source
            metadata_file = f"{source_name}_harmonized_metadata.{args.metadata_format}"
            metadata_path = os.path.join(args.output_dir, metadata_file)

            if args.metadata_format == 'json':
                with open(metadata_path, 'w', encoding='utf-8') as f:
                    json.dump(all_harmonized_records, f, indent=2, ensure_ascii=False)
            else:
                save_metadata_to_csv(all_harmonized_records, metadata_path)

            logging.info(f"Metadata saved to {metadata_path}")

            # Apply quality filtering
            filtered_records = []
            for record in all_harmonized_records:
                quality_score = record.get('quality_score', 0)

                # Check minimum quality score
                if quality_score < args.min_quality_score:
                    continue

                # Check AMR data requirement
                if args.require_amr_data:
                    amr_phenotypes = record.get('amr_phenotypes', [])
                    mic_data = record.get('mic_data', [])
                    antibiotic_resistance = record.get('antibiotic_resistance', [])
                    if not (amr_phenotypes or mic_data or antibiotic_resistance):
                        continue

                # Check location requirement
                if args.require_location:
                    country = record.get('country')
                    if not country:
                        continue

                filtered_records.append(record)

            logging.info(f"After quality filtering: {len(filtered_records)}/{len(all_harmonized_records)} genomes meet criteria")

            if not filtered_records:
                logging.warning("No genomes met the quality criteria. Try lowering requirements.")
                return

            # Update harmonized records to filtered set
            all_harmonized_records = filtered_records

            # Save filtered metadata
            metadata_file = f"{source_name}_filtered_harmonized_metadata.{args.metadata_format}"
            metadata_path = os.path.join(args.output_dir, metadata_file)

            if args.metadata_format == 'json':
                with open(metadata_path, 'w', encoding='utf-8') as f:
                    json.dump(all_harmonized_records, f, indent=2, ensure_ascii=False)
            else:
                save_metadata_to_csv(all_harmonized_records, metadata_path)

            logging.info(f"Filtered metadata saved to {metadata_path}")

            # Download FASTA sequences if requested
            if args.download:
                logging.info(f"Downloading FASTA sequences for quality-filtered genomes using {args.parallel_downloads} parallel workers...")
                
                # Use parallel download for all quality-filtered genomes
                successful = download_genomes_parallel(all_harmonized_records, args.output_dir, sources, args.parallel_downloads)
                logging.info(f"Parallel download completed: {len(successful)} high-quality FASTA sequences downloaded")

            # Print summary
            print_summary(all_harmonized_records, args.source)

    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)


def save_metadata_to_csv(records: List[Dict[str, Any]], output_file: str):
    """Save metadata to CSV format with enhanced field handling"""
    import csv

    # Debug: log what we're receiving
    logging.info(f"CSV save called with {len(records)} records")
    
    if not records:
        logging.warning("No records provided to CSV save function")
        return

    # Flatten nested structures for CSV with enhanced field handling
    flattened_data = []
    for item in records:
        flat_item = {}
        for key, value in item.items():
            if isinstance(value, list):
                if key == 'amr_phenotypes':
                    flat_item[key] = '; '.join(value) if value else ''
                elif key == 'mic_data':
                    # Format MIC data as readable strings
                    mic_strings = []
                    for mic in value:
                        if isinstance(mic, dict):
                            antibiotic = mic.get('antibiotic', 'Unknown')
                            mic_value = mic.get('mic_value', 'Unknown')
                            mic_strings.append(f"{antibiotic}: {mic_value}")
                        else:
                            mic_strings.append(str(mic))
                    flat_item[key] = '; '.join(mic_strings) if mic_strings else ''
                elif key == 'antibiotic_resistance':
                    # Format resistance data
                    resistance_strings = []
                    for res in value:
                        if isinstance(res, dict):
                            antibiotic = res.get('antibiotic', 'Unknown')
                            resistance = res.get('resistance', 'Unknown')
                            resistance_strings.append(f"{antibiotic}: {resistance}")
                        else:
                            resistance_strings.append(str(res))
                    flat_item[key] = '; '.join(resistance_strings) if resistance_strings else ''
                else:
                    # Convert other lists to semicolon-separated strings
                    flat_item[key] = '; '.join(str(v) for v in value) if value else ''
            elif isinstance(value, dict):
                # Flatten dictionaries with key_subkey naming
                for subkey, subvalue in value.items():
                    flat_item[f"{key}_{subkey}"] = str(subvalue) if subvalue is not None else ''
            else:
                # Handle regular values, ensuring None becomes empty string
                flat_item[key] = str(value) if value is not None else ''
        flattened_data.append(flat_item)

    # Debug: log flattening results
    logging.info(f"Flattened {len(records)} records into {len(flattened_data)} flat records")
    if flattened_data:
        logging.info(f"First flattened record has {len(flattened_data[0])} fields")

    # Ensure consistent field ordering for better CSV readability
    if flattened_data:
        # Define preferred field order
        priority_fields = [
            'accession', 'organism', 'strain', 'title', 'biosample', 'bioproject',
            'collection_date', 'country', 'host', 'isolation_source',
            'amr_phenotypes', 'mic_data', 'antibiotic_resistance',
            'quality_score', 'database', 'genome_id'
        ]
        
        # Get all available fields
        all_fields = set()
        for item in flattened_data:
            all_fields.update(item.keys())
        
        # Order fields: priority first, then alphabetical for the rest
        fieldnames = []
        for field in priority_fields:
            if field in all_fields:
                fieldnames.append(field)
                all_fields.remove(field)
        
        # Add remaining fields alphabetically
        fieldnames.extend(sorted(all_fields))
    else:
        fieldnames = []

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(flattened_data)

    logging.info(f"CSV metadata saved with {len(fieldnames)} columns and {len(flattened_data)} rows")


def print_summary(records: List[Dict[str, Any]], source: str):
    """Print a summary of the retrieved data"""
    if not records:
        return

    print(f"\n{'='*50}")
    print(f"Genome Harvest Summary - {source.upper()}")
    print(f"{'='*50}")

    total_genomes = len(records)

    # Count organisms
    organisms = {}
    countries = {}
    amr_count = 0

    for record in records:
        # Count organisms
        organism = record.get('organism', 'Unknown')
        organisms[organism] = organisms.get(organism, 0) + 1

        # Count countries
        country = record.get('country')
        if country:
            countries[country] = countries.get(country, 0) + 1

        # Count AMR phenotypes
        amr_phenotypes = record.get('amr_phenotypes', [])
        if amr_phenotypes:
            amr_count += 1

    print(f"Total genomes: {total_genomes}")
    print(f"Genomes with AMR data: {amr_count}")

    if organisms:
        print(f"\nTop organisms:")
        sorted_orgs = sorted(organisms.items(), key=lambda x: x[1], reverse=True)[:5]
        for org, count in sorted_orgs:
            print(f"  - {org}: {count}")

    if countries:
        print(f"\nGeographic distribution:")
        sorted_countries = sorted(countries.items(), key=lambda x: x[1], reverse=True)[:5]
        for country, count in sorted_countries:
            print(f"  - {country}: {count}")

    print(f"\nMetadata saved to: {source}_harmonized_metadata.json")
    print(f"{'='*50}")


if __name__ == '__main__':
    main()