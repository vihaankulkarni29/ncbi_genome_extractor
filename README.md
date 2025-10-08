# Federated Genome Harvester

A robust, multi-database tool for downloading genome sequences as FASTA files with **smart genome selection** and comprehensive metadata extraction. Supports antimicrobial resistance research and general genomic studies across NCBI, BV-BRC, EnteroBase, and PATRIC databases. Designed for bioinformaticians and researchers worldwide.

## Features

- **Multi-database genome retrieval**: Search across NCBI, BV-BRC, EnteroBase, and PATRIC
- **🎯 Smart genome selection**: Automatically selects highest-quality genomes based on metadata completeness
- **Quality filtering**: Configurable quality thresholds ensure only valuable genomes are downloaded
- Handle large files efficiently with streaming downloads
- Parse search URLs or custom queries from any supported database
- Batch download multiple genomes with progress tracking
- **Extract comprehensive metadata** including BioSample IDs, BioProject IDs, collection dates, organism information, and MIC data
- Export metadata in JSON or CSV formats for downstream analysis
- Robust error handling and automatic retries
- Command-line interface optimized for bioinformatics workflows
- **AMR-focused**: Special support for antimicrobial resistance research

## Installation

1. Clone this repository:
```bash
git clone https://github.com/vihaankulkarni/ncbi_genome_extractor.git
cd ncbi_genome_extractor
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. (Optional) Set up NCBI API key for faster downloads:
   - Get your API key from: https://www.ncbi.nlm.nih.gov/account/settings/
   - Update `NCBI_API_KEY` in `config.py`

## Usage

### Basic Usage

**Download genomes with automatic metadata extraction:**
```bash
python ncbi_genome_extractor.py --url "https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli+genome"
```

**Download from a custom query:**
```bash
python ncbi_genome_extractor.py --query "SARS-CoV-2 complete genome" --output_dir ./genomes
```

**Smart genome selection - get the BEST genomes with complete metadata:**
```bash
python ncbi_genome_extractor.py --query "Escherichia coli[Organism] AND complete genome" --max_results 10
```

**Clean, simple commands - metadata is extracted automatically!**

### 🎯 Smart Genome Selection

**The tool automatically selects the BEST genomes based on metadata quality:**

**Quality Scoring System (0-10 points):**
- ✅ **MIC Data**: 3 points (most valuable for AMR research)
- ✅ **BioSample ID**: 2 points
- ✅ **BioProject ID**: 1 point
- ✅ **Collection Date**: 1 point
- ✅ **Geographic Data**: 1 point (country/isolation source)
- ✅ **Host Information**: 1 point
- ✅ **Antibiotic Resistance**: 1 point

**How it works:**
1. Searches for 5x more genomes than requested (e.g., 50 genomes for `--max_results 10`)
2. Extracts metadata for all candidate genomes
3. Scores each genome based on metadata completeness
4. Selects the top N genomes with highest quality scores
5. Downloads only the highest-quality genomes

**Example Output:**
```
INFO - Evaluating 50 genomes for metadata quality...
INFO - Selected 10 best genomes:
INFO -   - High quality (score ≥7): 3
INFO -   - Medium quality (score 4-6): 5
INFO -   - Low quality (score <4): 2
```

**Disable smart selection (use first N results):**
```bash
python ncbi_genome_extractor.py --query "your query" --max_results 10 --no_quality_selection
```

### Handling Large Searches

**⚠️ Important for researchers:** NCBI searches can return thousands of genomes. Always use `--max_results` to control download size:

**For searches with 1000+ genomes:**
```bash
# Download only 50 genomes from a large search
python ncbi_genome_extractor.py --url "https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli" --max_results 50
```

**For focused research:**
```bash
# Get specific antibiotic-resistant strains
python ncbi_genome_extractor.py --query "Escherichia coli AND carbapenem resistance" --max_results 25
```

**Start small, then expand:**
```bash
# Test with 5 genomes first
python ncbi_genome_extractor.py --query "your complex query" --max_results 5

# Then scale up if needed
python ncbi_genome_extractor.py --query "your complex query" --max_results 100
```

### Advanced Options

**Essential for large searches:**
- `--max_results`: **Maximum number of genomes to download** (default: 100) - **CRITICAL** for controlling download size
- `--output_dir`: Directory to save FASTA files and metadata (default: ./output)

**Smart selection options:**
- `--no_quality_selection`: Disable quality-based genome selection (use first N results instead of best N)

**Metadata options:**
- `--metadata_format`: Format for metadata output (json or csv, default: json)
- `--no_metadata`: Skip metadata extraction (not recommended for AMR research)

**Performance & reliability:**
- `--retries`: Number of retry attempts for failed downloads (default: 3)
- `--delay`: Delay between requests in seconds (default: 0.5)
- `--log_level`: Logging level (DEBUG, INFO, WARNING, ERROR)

### Examples

**Download bacterial genomes:**
```bash
python ncbi_genome_extractor.py --query "Bacillus subtilis complete genome" --max_results 50
```

**Download viral genomes:**
```bash
python ncbi_genome_extractor.py --query "Influenza A virus genome" --output_dir ./viral_genomes
```

**Download from specific accession:**
```bash
python ncbi_genome_extractor.py --query "NC_045512" --max_results 1
```

**Complex queries with special characters:**
```bash
# For queries with parentheses and operators, use single quotes
python ncbi_genome_extractor.py --query '(Escherichia coli AND marcolide resistance AND Complete genome) AND "Escherichia coli"[porgn:__txid562]' --max_results 10 --output_dir ./ecoli_genomes

# Or escape special characters
python ncbi_genome_extractor.py --query "(Escherichia coli AND marcolide resistance AND Complete genome) AND \"Escherichia coli\"[porgn:__txid562]" --max_results 10 --output_dir ./ecoli_genomes
```

**Note:** When using complex queries with parentheses, operators (AND, OR, NOT), or special characters, wrap the entire query in single quotes to avoid shell interpretation issues.

### Automatic Metadata Extraction

**Metadata is extracted automatically with every genome download!** No extra flags needed.

**Simple command gets you both genomes AND comprehensive metadata:**
```bash
python ncbi_genome_extractor.py --query "Escherichia coli AND amikacin resistance" --max_results 10
```

**This automatically creates:**
- `output/*.fasta` - Genome sequence files
- `output/output_metadata.json` - Comprehensive metadata

**Metadata includes:**
- **BioSample ID**: Unique identifier for biological sample
- **BioProject ID**: Project identifier for related sequences
- **Collection date**: When the sample was collected
- **Organism information**: Genus and species names (auto-parsed)
- **Geographic data**: Country and isolation source
- **MIC data**: Minimum Inhibitory Concentration values for antibiotics
- **Antibiotic resistance**: Specific resistance phenotypes and genotypes

**Choose your preferred format:**
```bash
# JSON format (default)
python ncbi_genome_extractor.py --query "your query" --metadata_format json

# CSV format for Excel
python ncbi_genome_extractor.py --query "your query" --metadata_format csv
```

**Skip metadata only if you really need to:**
```bash
python ncbi_genome_extractor.py --query "your query" --no_metadata
```

**Example metadata JSON structure:**
```json
{
  "genome_id": "12345",
  "accession": "CP123456",
  "organism": "Escherichia coli",
  "genus": "Escherichia",
  "species": "coli",
  "biosample": "SAMN12345678",
  "bioproject": "PRJNA123456",
  "collection_date": "2020-01-15",
  "country": "USA",
  "mic_data": [
    {
      "antibiotic": "amikacin",
      "value": "16",
      "unit": "ug/ml"
    }
  ],
  "antibiotic_resistance": [
    {
      "antibiotic": "gentamicin",
      "resistance": "resistant"
    }
  ]
}
```

### API Usage

You can also use the tool programmatically:

```python
from ncbi_genome_extractor import NCBIExtractor

extractor = NCBIExtractor()

# Basic usage
genome_ids = extractor.search_genomes("your query")
extractor.download_multiple_genomes(genome_ids, "./output")

# With metadata extraction
genome_ids = extractor.search_genomes("Escherichia coli AND amikacin")
extractor.download_multiple_genomes(genome_ids, "./output")

# Extract metadata
metadata = extractor.extract_metadata(genome_ids)

# Save metadata
extractor.save_metadata_to_json(metadata, "genomes_metadata.json")
extractor.save_metadata_to_csv(metadata, "genomes_metadata.csv")
```

## Configuration

The tool uses your NCBI API key for faster downloads. Configure it in the script or environment variables.

## Requirements

- Python 3.7+
- NCBI API key (optional but recommended)

## License

MIT License

## Contributing

Contributions welcome! Please read the contributing guidelines before submitting pull requests.

## Contact

For issues or questions: vihaankulkarni29@gmail.com