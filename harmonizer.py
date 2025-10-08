#!/usr/bin/env python3
"""
Data Harmonizer - Standardizes genome metadata from different sources
"""

from typing import List, Dict, Any


def harmonize_data(raw_records: List[Dict[str, Any]], source: str) -> List[Dict[str, Any]]:
    """
    Harmonize raw genome records from different sources into a standardized schema.

    Args:
        raw_records: List of raw genome records from a specific source
        source: Source database ('ncbi' or 'bvbrc')

    Returns:
        List of harmonized genome records with standardized schema
    """
    harmonized_records = []

    for record in raw_records:
        if source == 'ncbi':
            harmonized_record = _harmonize_ncbi_record(record)
        elif source == 'bvbrc':
            harmonized_record = _harmonize_bvbrc_record(record)
        elif source == 'enterobase':
            harmonized_record = _harmonize_enterobase_record(record)
        elif source == 'patric':
            harmonized_record = _harmonize_patric_record(record)
        else:
            raise ValueError(f"Unknown source: {source}")

        # Add database source information
        harmonized_record['database'] = source.upper()

        harmonized_records.append(harmonized_record)

    return harmonized_records


def _harmonize_ncbi_record(record: Dict[str, Any]) -> Dict[str, Any]:
    """Harmonize NCBI record to standard schema"""
    harmonized = {
        'accession': record.get('accession', ''),
        'organism': record.get('organism', ''),
        'strain': _extract_strain_from_title(record.get('title', '')),
        'collection_date': record.get('collection_date', ''),
        'country': record.get('country', ''),
        'host': record.get('host', ''),
        'isolation_source': record.get('isolation_source', ''),
        'amr_phenotypes': _extract_amr_phenotypes(record)
    }

    # Clean up empty strings
    for key, value in harmonized.items():
        if isinstance(value, str) and not value.strip():
            harmonized[key] = None
        elif isinstance(value, list) and not value:
            harmonized[key] = []

    return harmonized


def _harmonize_bvbrc_record(record: Dict[str, Any]) -> Dict[str, Any]:
    """Harmonize BV-BRC record to standard schema"""
    # BV-BRC has different field names, map them to our standard schema
    harmonized = {
        'accession': record.get('accession', record.get('genome_id', '')),
        'organism': record.get('organism', record.get('species', '')),
        'strain': record.get('strain', record.get('isolate', '')),
        'collection_date': record.get('collection_date', record.get('isolation_date', '')),
        'country': record.get('country', record.get('geographic_location', '')),
        'host': record.get('host', record.get('host_name', '')),
        'isolation_source': record.get('isolation_source', record.get('source', '')),
        'amr_phenotypes': _extract_bvbrc_amr_phenotypes(record)
    }

    # Clean up empty strings
    for key, value in harmonized.items():
        if isinstance(value, str) and not value.strip():
            harmonized[key] = None
        elif isinstance(value, list) and not value:
            harmonized[key] = []

    return harmonized


def _extract_strain_from_title(title: str) -> str:
    """Extract strain name from NCBI sequence title"""
    if not title:
        return ''

    # Common patterns for strain extraction
    import re

    # Look for "strain XXX" pattern
    strain_match = re.search(r'strain\s+([^\s,]+)', title, re.IGNORECASE)
    if strain_match:
        return strain_match.group(1).strip()

    # Look for isolate patterns
    isolate_match = re.search(r'isolate\s+([^\s,]+)', title, re.IGNORECASE)
    if isolate_match:
        return isolate_match.group(1).strip()

    return ''


def _extract_amr_phenotypes(record: Dict[str, Any]) -> List[str]:
    """Extract AMR phenotypes from NCBI record"""
    phenotypes = []

    # Check resistance phenotypes
    if record.get('resistance_phenotype'):
        phenotypes.extend(record['resistance_phenotype'])

    # Check antibiotic resistance data
    if record.get('antibiotic_resistance'):
        for resistance in record['antibiotic_resistance']:
            if isinstance(resistance, dict) and resistance.get('antibiotic'):
                phenotypes.append(f"{resistance['antibiotic']}: {resistance.get('resistance', 'resistant')}")

    return phenotypes


def _harmonize_enterobase_record(record: Dict[str, Any]) -> Dict[str, Any]:
    """Harmonize EnteroBase record to standard schema"""
    harmonized = {
        'accession': record.get('accession', ''),
        'organism': record.get('organism', 'Escherichia coli'),
        'strain': record.get('strain', ''),
        'collection_date': record.get('collection_date', ''),
        'country': record.get('country', ''),
        'host': record.get('host', ''),
        'isolation_source': record.get('isolation_source', ''),
        'amr_phenotypes': _extract_enterobase_amr_phenotypes(record)
    }

    # Add EnteroBase-specific fields
    if record.get('serotype'):
        harmonized['serotype'] = record['serotype']
    if record.get('mlst_st'):
        harmonized['mlst_st'] = record['mlst_st']
    if record.get('quality_score'):
        harmonized['quality_score'] = record['quality_score']

    # Clean up empty strings
    for key, value in harmonized.items():
        if isinstance(value, str) and not value.strip():
            harmonized[key] = None
        elif isinstance(value, list) and not value:
            harmonized[key] = []

    return harmonized


def _harmonize_patric_record(record: Dict[str, Any]) -> Dict[str, Any]:
    """Harmonize PATRIC record to standard schema"""
    harmonized = {
        'accession': record.get('accession', record.get('genome_id', '')),
        'organism': record.get('organism', ''),
        'strain': record.get('strain', ''),
        'collection_date': record.get('collection_date', ''),
        'country': record.get('country', ''),
        'host': record.get('host', ''),
        'isolation_source': record.get('isolation_source', ''),
        'amr_phenotypes': _extract_patric_amr_phenotypes(record)
    }

    # Add PATRIC-specific fields
    if record.get('mlst_st'):
        harmonized['mlst_st'] = record['mlst_st']
    if record.get('genome_length'):
        harmonized['genome_length'] = record['genome_length']
    if record.get('quality_score'):
        harmonized['quality_score'] = record['quality_score']

    # Clean up empty strings
    for key, value in harmonized.items():
        if isinstance(value, str) and not value.strip():
            harmonized[key] = None
        elif isinstance(value, list) and not value:
            harmonized[key] = []

    return harmonized


def _extract_enterobase_amr_phenotypes(record: Dict[str, Any]) -> List[str]:
    """Extract AMR phenotypes from EnteroBase record"""
    phenotypes = []

    # EnteroBase has rich AMR data
    if record.get('resistance_phenotype'):
        phenotypes.extend(record['resistance_phenotype'])

    # Add MIC-based phenotypes
    if record.get('mic_data'):
        for mic_entry in record['mic_data']:
            antibiotic = mic_entry.get('antibiotic', '')
            value = mic_entry.get('value', '')
            if antibiotic and value:
                phenotypes.append(f"{antibiotic} MIC: {value}")

    # Add antibiotic resistance entries
    if record.get('antibiotic_resistance'):
        for resistance in record['antibiotic_resistance']:
            if isinstance(resistance, dict):
                antibiotic = resistance.get('antibiotic', '')
                res_type = resistance.get('resistance', '')
                if antibiotic and res_type:
                    phenotypes.append(f"{antibiotic}: {res_type}")
            elif isinstance(resistance, str):
                phenotypes.append(resistance)

    return phenotypes


def _extract_patric_amr_phenotypes(record: Dict[str, Any]) -> List[str]:
    """Extract AMR phenotypes from PATRIC record"""
    phenotypes = []

    # PATRIC AMR data structure
    if record.get('antibiotic_resistance'):
        for resistance in record['antibiotic_resistance']:
            if isinstance(resistance, dict):
                antibiotic = resistance.get('antibiotic', '')
                res_type = resistance.get('resistance', '')
                if antibiotic and res_type:
                    phenotypes.append(f"{antibiotic}: {res_type}")
            elif isinstance(resistance, str):
                phenotypes.append(resistance)

    # Add MIC data if available
    if record.get('mic_data'):
        for mic_entry in record['mic_data']:
            antibiotic = mic_entry.get('antibiotic', '')
            value = mic_entry.get('value', '')
            if antibiotic and value:
                phenotypes.append(f"{antibiotic} MIC: {value}")

    return phenotypes


def _extract_bvbrc_amr_phenotypes(record: Dict[str, Any]) -> List[str]:
    """Extract AMR phenotypes from BV-BRC record"""
    phenotypes = []

    # BV-BRC has different AMR data structure
    # This will need to be customized based on BV-BRC API response format
    amr_data = record.get('amr', record.get('antibiotic_resistance', []))

    if isinstance(amr_data, list):
        for item in amr_data:
            if isinstance(item, dict):
                antibiotic = item.get('antibiotic', item.get('drug', ''))
                resistance = item.get('resistance', item.get('phenotype', ''))
                if antibiotic and resistance:
                    phenotypes.append(f"{antibiotic}: {resistance}")
            elif isinstance(item, str):
                phenotypes.append(item)

    return phenotypes