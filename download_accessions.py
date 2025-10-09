#!/usr/bin/env python3
"""
Download genomes from accession list
"""

import os
import sys
from clients.ncbi_client import NCBIClient

def main():
    if len(sys.argv) != 2:
        print("Usage: python download_accessions.py <accession_file>")
        sys.exit(1)

    accession_file = sys.argv[1]
    output_dir = "./accession_list_genomes"

    # Read accessions
    with open(accession_file, 'r') as f:
        accessions = [line.strip() for line in f if line.strip()]

    print(f"Found {len(accessions)} accessions to download")

    # Initialize extractor
    extractor = NCBIClient()

    # Download in batches to avoid API limits
    batch_size = 100
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        print(f"Downloading batch {i//batch_size + 1}/{(len(accessions) + batch_size - 1)//batch_size}")

        # Download FASTA files
        for acc in batch:
            try:
                success = extractor.download_fasta(acc, output_dir)
                if success:
                    print(f"[OK] Downloaded {acc}")
                else:
                    print(f"[FAIL] Failed {acc}")
            except Exception as e:
                print(f"[ERROR] Error downloading {acc}: {e}")

    # Extract metadata for all downloaded genomes
    downloaded_files = [f for f in os.listdir(output_dir) if f.endswith('.fasta')]
    genome_ids = [f.replace('.fasta', '') for f in downloaded_files]

    print(f"Extracting metadata for {len(genome_ids)} genomes...")
    # Use fetch_genomes to get metadata (it will internally call extract_metadata)
    metadata_records = extractor.fetch_genomes(','.join(genome_ids), len(genome_ids))

    # Save metadata to JSON
    import json
    with open(os.path.join(output_dir, "accession_metadata.json"), 'w') as f:
        json.dump(metadata_records, f, indent=2)

    print(f"Downloaded {len(downloaded_files)} genomes to {output_dir}")

if __name__ == "__main__":
    main()