"""Command-line interface and argument parsing for PHAGPCR."""

import logging
import os
from argparse import ArgumentParser, Namespace
from typing import Optional

from phagpcr_package.constants import DEFAULT_BLAST_DB

logger = logging.getLogger(__name__)


def parse_args() -> Namespace:
    """Parse command line arguments."""
    parser = ArgumentParser(
        description='Design primers using Primer3 and perform '
                    'specificity screening')
    parser.add_argument('-f', '--fasta_file',
                        type=str,
                        help='Path to the input FASTA file')
    parser.add_argument('-o', '--output_dir',
                        type=str,
                        help='Directory where the results will be saved')
    parser.add_argument('-b', '--blast_db',
                        type=str,
                        default=DEFAULT_BLAST_DB,
                        required=False,
                        help='BLAST database for primer screening '
                             '(optional)')
    parser.add_argument('-hg', '--human_genome',
                        action='store_true',
                        help='Screen the primers against the human '
                             'genome - version hg38 (optional)')
    parser.add_argument('-s', '--sequences_dir',
                        type=str, default=None, required=False,
                        help='Directory with sequences for screening '
                             '(optional)')
    parser.add_argument('-r', '--runtype',
                        type=int, default=1, choices=[1, 2, 3],
                        help='Run type: 1=supervised, '
                             '2=multiplex/cocktail, 3=unsupervised')
    parser.add_argument('-nb', '--primer_nb',
                        type=int, default=1, required=False,
                        help='Number of primer pairs per target sequence')
    parser.add_argument('-tm', '--tm_optimal',
                        type=int, default=60, required=False,
                        help='Optimal Tm for primers')
    parser.add_argument('-sz', '--size_optimal',
                        type=int, default=20, required=False,
                        help='Optimal size for primers')
    parser.add_argument('-k', '--kit',
                        type=str, default="",
                        choices=["Quantinova", "Invitrogen"],
                        required=False,
                        help='qPCR kit for Tm specification')
    parser.add_argument('-tm2', '--tm_specificity',
                        type=int, default=40, required=False,
                        help='Minimum Tm for MFEprimer matches')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Enable verbose logging output')
    return parser.parse_args()


def setup_logging(verbose: bool = False,
                  log_file: Optional[str] = None) -> None:
    """
    Configure logging for the application.

    Args:
        verbose: If True, set logging level to DEBUG
        log_file: Optional path to log file
    """
    level = logging.DEBUG if verbose else logging.INFO
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'

    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=level,
        format=log_format,
        datefmt=date_format,
        handlers=handlers,
        force=True
    )


def validate_args(args: Namespace) -> None:
    """
    Validate command line arguments.

    Args:
        args: Parsed command line arguments

    Raises:
        ValueError: If required arguments are missing or invalid
        FileNotFoundError: If input files don't exist
        PermissionError: If files aren't readable
    """
    if not args.fasta_file:
        raise ValueError('Fasta file path (-f/--fasta_file) is required')
    if not args.output_dir:
        raise ValueError('Output directory (-o/--output_dir) is required')

    if not os.path.exists(args.fasta_file):
        raise FileNotFoundError(f'Fasta file not found: '
                                f'{args.fasta_file}')
    if not os.access(args.fasta_file, os.R_OK):
        raise PermissionError(f'Cannot read fasta file: '
                              f'{args.fasta_file}')
