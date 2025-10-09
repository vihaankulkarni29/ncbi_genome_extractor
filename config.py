# Configuration for NCBI Genome Extractor

# NCBI Entrez API settings
NCBI_EMAIL = "vihaankulkarni29@gmail.com"
NCBI_API_KEY = "0a4a31f50b826e949816d1ba53c684546608"

# Default settings
DEFAULT_MAX_RESULTS = 100
DEFAULT_RETRIES = 3
DEFAULT_DELAY = 0.5
DEFAULT_OUTPUT_DIR = "./output"

# URLs
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_SEARCH_URL = NCBI_BASE_URL + "esearch.fcgi"
NCBI_FETCH_URL = NCBI_BASE_URL + "efetch.fcgi"

# BV-BRC API settings
BV_BRC_BASE_URL = "https://www.bv-brc.org/api"
BV_BRC_API_KEY = None  # Set this if you have a BV-BRC API key

# EnteroBase API settings
ENTEROBASE_BASE_URL = "https://enterobase.warwick.ac.uk/api/v2.0"
ENTEROBASE_API_KEY = None  # Set this if you have an EnteroBase API key

# PATRIC API settings
PATRIC_BASE_URL = "https://www.patricbrc.org/api"
PATRIC_API_KEY = None  # Set this if you have a PATRIC API key

# BioCyc API settings
BIOCYC_BASE_URL = "https://websvc.biocyc.org"
BIOCYC_API_KEY = None  # Set this if you have a BioCyc API key