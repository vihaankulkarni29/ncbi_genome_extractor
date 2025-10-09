#!/usr/bin/env python3
"""
Analyze contigs for AMR genes and metadata
"""

import os
import json
import subprocess
import pandas as pd
from pathlib import Path

def run_abricate_on_contigs(contig_dir: str, output_file: str):
    """Run ABRicate on all contig FASTA files"""
    contig_files = list(Path(contig_dir).glob("*.fasta"))

    if not contig_files:
        print(f"No FASTA files found in {contig_dir}")
        return None

    print(f"Running ABRicate on {len(contig_files)} contig files...")

    # Run ABRicate on all files
    cmd = ["abricate", "--db", "card"] + [str(f) for f in contig_files]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        with open(output_file, 'w') as f:
            f.write(result.stdout)
        print(f"ABRicate results saved to {output_file}")
        return output_file
    else:
        print(f"ABRicate failed: {result.stderr}")
        return None

def extract_metadata_from_contigs(contig_dir: str, metadata_file: str):
    """Extract metadata from contig FASTA headers"""
    contig_files = list(Path(contig_dir).glob("*.fasta"))
    metadata_records = []

    for fasta_file in contig_files:
        accession = fasta_file.stem  # Remove .fasta extension

        # Read the FASTA file to get header information
        try:
            with open(fasta_file, 'r') as f:
                header = f.readline().strip()
                if header.startswith('>'):
                    # Extract information from FASTA header
                    description = header[1:]  # Remove '>' character

                    # Create metadata record
                    record = {
                        'accession': accession,
                        'fasta_header': description,
                        'file_path': str(fasta_file),
                        'organism': 'Unknown',  # Would need additional parsing
                        'contig_type': 'assembly_contig',
                        'source': 'NCBI_contig_download'
                    }
                    metadata_records.append(record)

        except Exception as e:
            print(f"Error reading {fasta_file}: {e}")

    # Save metadata
    with open(metadata_file, 'w') as f:
        json.dump(metadata_records, f, indent=2)

    print(f"Metadata extracted for {len(metadata_records)} contigs")
    return metadata_file

def combine_amr_and_metadata(abricate_file: str, metadata_file: str, output_file: str):
    """Combine ABRicate results with metadata"""
    try:
        # Read ABRicate results
        abricate_df = pd.read_csv(abricate_file, sep='\t')

        # Read metadata
        with open(metadata_file, 'r') as f:
            metadata_records = json.load(f)

        # Create metadata lookup
        metadata_dict = {record['accession']: record for record in metadata_records}

        # Add metadata columns to ABRicate results
        abricate_df['organism'] = abricate_df['#FILE'].apply(
            lambda x: metadata_dict.get(Path(x).stem, {}).get('organism', 'Unknown')
        )
        abricate_df['contig_type'] = abricate_df['#FILE'].apply(
            lambda x: metadata_dict.get(Path(x).stem, {}).get('contig_type', 'Unknown')
        )

        # Save combined results
        abricate_df.to_csv(output_file, index=False)
        print(f"Combined AMR and metadata results saved to {output_file}")

        # Summary statistics
        total_contigs = len(metadata_records)
        amr_positive_contigs = abricate_df['#FILE'].nunique()
        total_amr_hits = len(abricate_df)

        print("\n📊 Analysis Summary:")
        print(f"  - Total contigs analyzed: {total_contigs}")
        print(f"  - Contigs with AMR genes: {amr_positive_contigs}")
        print(f"  - Total AMR gene hits: {total_amr_hits}")

        if total_amr_hits > 0:
            print("  - Top AMR genes found:")
            top_genes = abricate_df['GENE'].value_counts().head(5)
            for gene, count in top_genes.items():
                print(f"    • {gene}: {count} hits")

        return output_file

    except Exception as e:
        print(f"Error combining results: {e}")
        return None

def main():
    # Configuration
    contig_dir = "./accession_list_genomes"  # Directory with contig FASTA files
    abricate_output = "./contig_amr_results.tab"
    metadata_output = "./contig_metadata.json"
    combined_output = "./contig_amr_with_metadata.csv"

    print("🔬 Analyzing contigs for AMR genes and metadata...")

    # Check if contig directory exists
    if not os.path.exists(contig_dir):
        print(f"❌ Contig directory {contig_dir} not found!")
        print("💡 To analyze contigs:")
        print("   1. Download contigs using download_accessions.py")
        print("   2. Run this script")
        return

    # Step 1: Run ABRicate analysis
    print("\n1️⃣ Running ABRicate AMR analysis...")
    abricate_result = run_abricate_on_contigs(contig_dir, abricate_output)

    # Step 2: Extract metadata from contigs
    print("\n2️⃣ Extracting metadata from contig headers...")
    metadata_result = extract_metadata_from_contigs(contig_dir, metadata_output)

    # Step 3: Combine results
    if abricate_result and metadata_result:
        print("\n3️⃣ Combining AMR and metadata results...")
        combine_amr_and_metadata(abricate_result, metadata_result, combined_output)

    print("\n✅ Contig analysis complete!")
    print("📁 Output files:")
    print(f"   • AMR results: {abricate_output}")
    print(f"   • Metadata: {metadata_output}")
    print(f"   • Combined analysis: {combined_output}")

if __name__ == "__main__":
    main()