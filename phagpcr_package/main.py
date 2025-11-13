#!/usr/bin/env python3
"""Main entry point for PHAGPCR primer design pipeline."""

import sys

from phagpcr_package.cli import parse_args, setup_logging, validate_args
from phagpcr_package.io_utils import mkdir_outdir, readfile
from phagpcr_package.workflows import run_workflow_1_and_2, run_workflow_3


def main() -> int:
    """
    Main function to run PHAGPCR pipeline.

    Returns:
        Exit code (0 for success, 1 for failure)
    """
    try:
        # Parse arguments
        args = parse_args()

        # Validate arguments
        validate_args(args)

        # Create output directories first
        mkdir_outdir(args.output_dir)
        mkdir_outdir(args.output_dir + "/indiv_results")

        # Setup logging (after directory creation)
        log_file = None
        if args.output_dir:
            log_file = f"{args.output_dir}/phagpcr.log"
        setup_logging(verbose=args.verbose, log_file=log_file)

        # Read input sequences
        sequences = readfile(args.fasta_file)

        # Run appropriate workflow
        if args.runtype in [1, 2]:
            run_workflow_1_and_2(sequences, args)
        elif args.runtype == 3:
            run_workflow_3(sequences, args)

        return 0

    except Exception as e:
        import logging
        logger = logging.getLogger(__name__)
        logger.exception(f"Pipeline failed: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
