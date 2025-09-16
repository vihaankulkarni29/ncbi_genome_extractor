# Configuration for NCBI Genome Extractor

# NCBI Entrez API settings
NCBI_EMAIL = "vihaankulkarni29@gmail.com"
NCBI_API_KEY = "ef7622c2e716fa317fe04d24c42904211107"

# Default settings
DEFAULT_MAX_RESULTS = 100
DEFAULT_RETRIES = 3
DEFAULT_DELAY = 0.5
DEFAULT_OUTPUT_DIR = "./output"

# URLs
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_SEARCH_URL = NCBI_BASE_URL + "esearch.fcgi"
NCBI_FETCH_URL = NCBI_BASE_URL + "efetch.fcgi"