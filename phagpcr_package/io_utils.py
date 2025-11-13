"""File I/O utilities for PHAGPCR."""

import logging
import os
from typing import Dict

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

logger = logging.getLogger(__name__)


def mkdir_outdir(output_dir: str) -> None:
    """Create output directory if it doesn't exist."""
    # Use print if logger not yet configured (during initialization)
    if logger.hasHandlers():
        logger.info(f'Creating output directory: {output_dir}')
    if not os.path.exists('{}'.format(output_dir)):
        os.makedirs('{}'.format(output_dir))
    elif logger.hasHandlers():
        logger.warning(f'Output directory already exists: {output_dir}')


def readfile(fasta_file: str) -> Dict[str, str]:
    """Read FASTA file into a dictionary of sequences."""
    logger.info(f'Reading data from: {fasta_file}')
    sequences = {}
    try:
        with open(fasta_file, 'r') as f:
            current_seq = []
            current_header = ''
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_header:
                        sequences[current_header] = ''.join(current_seq)
                    current_header = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_header:
                sequences[current_header] = ''.join(current_seq)
        return sequences
    except IOError as e:
        raise IOError(f"Error reading file {fasta_file}: {str(e)}")


def df_to_fasta(df: pd.DataFrame, output_dir: str,
                output_fasta: str) -> None:
    """Convert DataFrame of primers to FASTA format."""
    left_records = [
        SeqRecord(Seq(row['left_primer_sequence']),
                  id=row['left_primer_name'], description="")
        for index, row in df.iterrows()
    ]
    right_records = [
        SeqRecord(Seq(row['right_primer_sequence']),
                  id=row['right_primer_name'], description="")
        for index, row in df.iterrows()
    ]
    records = left_records + right_records
    with open(f"{output_dir}/{output_fasta}.fna", "w") as handle:
        SeqIO.write(records, handle, 'fasta')
