import unittest
from unittest.mock import patch, MagicMock
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from ncbi_genome_extractor import NCBIExtractor

class TestNCBIExtractor(unittest.TestCase):
    def setUp(self):
        self.extractor = NCBIExtractor()

    @patch('requests.get')
    def test_search_query(self, mock_get):
        # Mock response for search
        mock_response = MagicMock()
        mock_response.text = '''<?xml version="1.0"?>
        <eSearchResult>
            <Count>2</Count>
            <IdList>
                <Id>12345</Id>
                <Id>67890</Id>
            </IdList>
        </eSearchResult>'''
        mock_get.return_value = mock_response

        ids = self.extractor.search_genomes("test query")
        self.assertEqual(ids, ['12345', '67890'])

    def test_parse_url(self):
        url = "https://www.ncbi.nlm.nih.gov/nuccore/?term=test+query"
        query = self.extractor.parse_search_url(url)
        self.assertEqual(query, "test query")

if __name__ == '__main__':
    unittest.main()