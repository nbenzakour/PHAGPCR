"""
PHAGPCR: qPCR primer design for phage detection in patient serum.

This package provides tools for designing and evaluating qPCR primers
using Primer3 and MFEprimer for specificity screening.
"""

__version__ = "1.0.0"
__author__ = "Nouri L. Ben Zakour"

from phagpcr_package.primer_design import (
    update_primer3,
    get_primers,
    parse_primers,
)
from phagpcr_package.mfe_utils import (
    run_mfe_index,
    run_mfe,
    run_mfe_dimers,
    run_mfe_spec,
    parse_MFEprimers_file,
    parse_MFEprimers_dimer_file,
)
from phagpcr_package.io_utils import (
    mkdir_outdir,
    readfile,
    df_to_fasta,
)

__all__ = [
    # Primer design
    'update_primer3',
    'get_primers',
    'parse_primers',
    # MFE utilities
    'run_mfe_index',
    'run_mfe',
    'run_mfe_dimers',
    'run_mfe_spec',
    'parse_MFEprimers_file',
    'parse_MFEprimers_dimer_file',
    # I/O utilities
    'mkdir_outdir',
    'readfile',
    'df_to_fasta',
]
