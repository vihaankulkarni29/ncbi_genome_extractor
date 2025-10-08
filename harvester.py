#!/usr/bin/env python3
"""
Federated Genome Harvester - Multi-source genome data retrieval tool
"""

import argparse
import json
import logging
import os
import sys
from typing import List, Dict, Any

from clients.ncbi_client import NCBIClient
from clients.bvbrc_client import BVBRCClient
from harmonizer import harmonize_data
from config import DEFAULT_OUTPUT_DIR


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
            logging.info(f"Downloading {len(genome_ids)} genomes from {args.source.upper()}...")

            successful = []
            for genome_id in genome_ids:
                genome_id = genome_id.strip()
                # Try to download from available sources
                for source in sources:
                    if source == 'ncbi':
                        client = NCBIClient()
                        result = client.download_fasta(genome_id, args.output_dir)
                        if isinstance(result, tuple) and len(result) == 2:
                            success, accession = result
                            if success:
                                successful.append(accession)
                                break  # Found it, no need to try other sources
                        else:
                            if result:
                                successful.append(genome_id)
                                break
                    elif source == 'bvbrc':
                        client = BVBRCClient()
                        result = client.download_fasta(genome_id, args.output_dir)
                        if isinstance(result, tuple) and len(result) == 2:
                            success, accession = result
                            if success:
                                successful.append(accession)
                                break
                        else:
                            if result:
                                successful.append(genome_id)
                                break
                    elif source == 'enterobase':
                        from clients.enterobase_client import EnteroBaseClient
                        client = EnteroBaseClient()
                        result = client.download_fasta(genome_id, args.output_dir)
                        if isinstance(result, tuple) and len(result) == 2:
                            success, accession = result
                            if success:
                                successful.append(accession)
                                break
                    elif source == 'patric':
                        from clients.patric_client import PATRICClient
                        client = PATRICClient()
                        result = client.download_fasta(genome_id, args.output_dir)
                        if isinstance(result, tuple) and len(result) == 2:
                            success, accession = result
                            if success:
                                successful.append(accession)
                                break

            logging.info(f"Downloaded {len(successful)} genomes successfully")

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
                logging.info("Downloading FASTA sequences for quality-filtered genomes...")
                successful = []

                for record in all_harmonized_records:
                    accession = record.get('accession', '')
                    database = record.get('database', '').lower()
                    if accession and database in sources:
                        if database == 'ncbi':
                            client = NCBIClient()
                            result = client.download_fasta(accession, args.output_dir)
                            if isinstance(result, tuple) and len(result) == 2:
                                success, final_accession = result
                                if success:
                                    successful.append(final_accession)
                            else:
                                if result:
                                    successful.append(accession)
                        elif database == 'bvbrc':
                            client = BVBRCClient()
                            result = client.download_fasta(accession, args.output_dir)
                            if isinstance(result, tuple) and len(result) == 2:
                                success, final_accession = result
                                if success:
                                    successful.append(final_accession)
                            else:
                                if result:
                                    successful.append(accession)
                        elif database == 'enterobase':
                            from clients.enterobase_client import EnteroBaseClient
                            client = EnteroBaseClient()
                            result = client.download_fasta(accession, args.output_dir)
                            if isinstance(result, tuple) and len(result) == 2:
                                success, final_accession = result
                                if success:
                                    successful.append(final_accession)
                        elif database == 'patric':
                            from clients.patric_client import PATRICClient
                            client = PATRICClient()
                            result = client.download_fasta(accession, args.output_dir)
                            if isinstance(result, tuple) and len(result) == 2:
                                success, final_accession = result
                                if success:
                                    successful.append(final_accession)

                logging.info(f"Downloaded {len(successful)} high-quality FASTA sequences")

            # Print summary
            print_summary(all_harmonized_records, args.source)

    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)


def save_metadata_to_csv(records: List[Dict[str, Any]], output_file: str):
    """Save metadata to CSV format"""
    import csv

    if not records:
        return

    # Flatten nested structures for CSV
    flattened_data = []
    for item in records:
        flat_item = {}
        for key, value in item.items():
            if isinstance(value, list):
                if key == 'amr_phenotypes':
                    flat_item[key] = '; '.join(value) if value else ''
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