"""Primer design utilities using Primer3."""

import logging
from typing import Dict

import pandas as pd
import primer3

from phagpcr_package.constants import (
    SEQUENCE_SUBREGION_MARGIN,
    MIN_SEQUENCE_FOR_SUBREGION,
    TARGET_NAME_MAX_LENGTH
)

logger = logging.getLogger(__name__)


def update_primer3(tm: int, primer_size: int, primer_nb: int,
                   kit: str) -> Dict:
    """Configure Primer3 parameters for primer design."""
    primer3_params = {
        # Basic parameters
        'PRIMER_OPT_SIZE': primer_size,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 23,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_OPT_TM': float(tm),
        'PRIMER_MIN_TM': 59.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_PAIR_MAX_DIFF_TM': 100.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        # Product size ranges
        'PRIMER_PRODUCT_SIZE_RANGE': [
            [50, 120], [120, 180], [180, 250], [250, 350], [350, 500]
        ],
        'PRIMER_NUM_RETURN': primer_nb,
        # Thermodynamic parameters
        'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
        'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 0,
        'PRIMER_MAX_LIBRARY_MISPRIMING': 12.00,
        'PRIMER_PAIR_MAX_LIBRARY_MISPRIMING': 20.00,
        'PRIMER_MAX_SELF_ANY_TH': 45.0,
        'PRIMER_MAX_SELF_END_TH': 35.0,
        'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0,
        'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,
        'PRIMER_MAX_HAIRPIN_TH': 24.0,
        'PRIMER_MAX_END_STABILITY': 9.0,
    }

    # Kit-specific adjustments
    if kit == "Invitrogen":
        primer3_params['PRIMER_MIN_TM'] = 58.0
        primer3_params['PRIMER_MAX_TM'] = 61.0

    return primer3_params


def get_primers(header: str, seq: str, primer3_params: Dict) -> Dict:
    """Design primers for given sequence using Primer3."""
    logger.info(f'Designing primers for: {header}')
    # For long sequences, define subregion to avoid terminal regions
    if len(seq) >= MIN_SEQUENCE_FOR_SUBREGION:
        start = SEQUENCE_SUBREGION_MARGIN
        region_length = len(seq) - (start * 2)
        primer3_params['SEQUENCE_INCLUDED_REGION'] = [start,
                                                       region_length]
    else:
        primer3_params['SEQUENCE_INCLUDED_REGION'] = [0, len(seq)]

    # Run Primer3
    try:
        results = primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': header,
                'SEQUENCE_TEMPLATE': seq,
            },
            primer3_params
        )
        return results
    except Exception as e:
        raise RuntimeError(f"Primer3 design failed for {header}: "
                           f"{str(e)}")


def parse_primers(primer_results: Dict, header: str) -> pd.DataFrame:
    """Parse Primer3 results into a DataFrame."""
    # Truncate target name if too long
    target_name = header[:TARGET_NAME_MAX_LENGTH]

    primers = []
    i = 0
    num_returned = primer_results['PRIMER_PAIR_NUM_RETURNED']
    while i < num_returned:
        left_key = f'PRIMER_LEFT_{i}'
        right_key = f'PRIMER_RIGHT_{i}'
        pair_key = f'PRIMER_PAIR_{i}'

        primers.append({
            'header': header,
            'left_primer_name': f"{target_name}_{i}F",
            'left_primer_sequence':
                primer_results[f'{left_key}_SEQUENCE'],
            'right_primer_name': f"{target_name}_{i}R",
            'right_primer_sequence':
                primer_results[f'{right_key}_SEQUENCE'],
            'left position': primer_results[left_key],
            'right position': primer_results[right_key],
            'product size': primer_results[f'{pair_key}_PRODUCT_SIZE'],
            'left_Tm': primer_results[f'{left_key}_TM'],
            'right_Tm': primer_results[f'{right_key}_TM'],
            'left_GC': primer_results[f'{left_key}_GC_PERCENT'],
            'right_GC': primer_results[f'{right_key}_GC_PERCENT'],
            'left_SELF_ANY_TH':
                primer_results[f'{left_key}_SELF_ANY_TH'],
            'right_SELF_ANY_TH':
                primer_results[f'{right_key}_SELF_ANY_TH'],
            'left_SELF_END_TH':
                primer_results[f'{left_key}_SELF_END_TH'],
            'right_SELF_END_TH':
                primer_results[f'{right_key}_SELF_END_TH'],
            'left_HAIRPIN_TH': primer_results[f'{left_key}_HAIRPIN_TH'],
            'right_HAIRPIN_TH':
                primer_results[f'{right_key}_HAIRPIN_TH'],
            'left_END_STABILITY':
                primer_results[f'{left_key}_END_STABILITY'],
            'right_END_STABILITY':
                primer_results[f'{right_key}_END_STABILITY'],
            'left_PENALTY': primer_results[f'{left_key}_PENALTY'],
            'right_PENALTY': primer_results[f'{right_key}_PENALTY'],
            'pair_PENALTY': primer_results[f'{pair_key}_PENALTY'],
            'pair_COMPL_ANY_TH':
                primer_results[f'{pair_key}_COMPL_ANY_TH'],
            'pair_COMPL_END_TH':
                primer_results[f'{pair_key}_COMPL_END_TH'],
        })
        i += 1
    primers_df = pd.DataFrame(primers)
    primers_df = primers_df.round(2)
    return primers_df
