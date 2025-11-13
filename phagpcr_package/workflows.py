"""Workflow orchestration for PHAGPCR primer design pipelines."""

import logging
import os
from argparse import Namespace
from typing import Dict

import pandas as pd

from phagpcr_package.primer_design import (
    update_primer3,
    get_primers,
    parse_primers
)
from phagpcr_package.mfe_utils import (
    run_mfe_index,
    run_mfe,
    run_mfe_dimers,
    run_mfe_spec,
    parse_MFEprimers_file,
    parse_MFEprimers_dimer_file
)
from phagpcr_package.io_utils import (
    mkdir_outdir,
    readfile,
    df_to_fasta
)
from phagpcr_package.constants import DEFAULT_HG38_PATH

logger = logging.getLogger(__name__)


def run_workflow_1_and_2(sequences: Dict[str, str], args: Namespace) -> None:
    """
    Run workflow 1 (TARGET) or 2 (COCKTAIL) primer design.

    Args:
        sequences: Dictionary of sequence headers and sequences
        args: Parsed command line arguments
    """
    output_dir = args.output_dir
    blast_db = args.blast_db
    human_genome = args.human_genome
    runtype = args.runtype
    primer_nb = args.primer_nb
    tm = args.tm_optimal
    primer_size = args.size_optimal
    kit = args.kit
    tm_spec = args.tm_specificity

    logger.info("=" * 60)
    logger.info("Step 1: Select primers against target using Primer3")
    logger.info("=" * 60)

    primers = []
    for header, seq in sequences.items():
        primer3_params = update_primer3(tm, primer_size, primer_nb, kit)
        primer_results = get_primers(header, seq, primer3_params)
        prefiltered_primers_df = parse_primers(primer_results, header)

    # Save prefiltered primers
    prefiltered_primers_df.to_csv(
        output_dir + '/prefiltered_primers_file.csv',
        index=False, sep='\t')
    logger.info(f'Prefiltered primers stored in: '
                f'{output_dir}/prefiltered_primers_file.csv')

    # Create multifasta file
    df_to_fasta(prefiltered_primers_df, output_dir,
                'prefiltered_primers')
    logger.info(f'Prefiltered primers in fasta format stored in: '
                f'{output_dir}/prefiltered_primers.fna')

    logger.info("=" * 60)
    logger.info("Step 2: Test specificity using MFEprimer against "
                "custom database")
    logger.info("=" * 60)

    if blast_db is not None:
        run_mfe_index(blast_db)
        logger.info(f"Screening for all features and specificity "
                    f"against Blast database: {blast_db}")

        # Run MFE on separate fasta files
        grouped = prefiltered_primers_df.groupby('left_primer_name')
        df_local_screening = pd.DataFrame()
        for name, group in grouped:
            df_to_fasta(grouped.get_group(name), output_dir, name[:-1])
            run_mfe(output_dir + '/' + str(name[:-1]) + '.fna',
                    blast_db, output_dir + '/indiv_results',
                    str(name[:-1]), tm_spec)
            local_screening = parse_MFEprimers_file(
                output_dir + '/indiv_results/mfe_results_' +
                str(name[:-1]) + '.txt')
            if local_screening is not None and not local_screening.empty:
                df_local_screening = pd.concat([df_local_screening,
                                                 local_screening],
                                                ignore_index=True)
            os.remove(output_dir + '/' + str(name[:-1]) + '.fna')

        print(df_local_screening)
        if not df_local_screening.empty:
            df_local_screening.to_csv(
                output_dir +
                '/mfe_screening_results_against_local_database_summary.csv',
                sep='\t', index=False)
            if 'Fp' in df_local_screening.columns:
                print(df_local_screening.groupby(
                    df_local_screening['Fp'].str[:-1])
                    ['Max_hit_per_pair'].unique())
            else:
                logger.warning("==> No 'Fp' column found in screening "
                               "results")
        else:
            logger.warning("==> No screening results to save. "
                           "DataFrame is empty.")
    else:
        logger.warning("No Blast database provided")

    logger.info("=" * 60)
    logger.info("Step 3: Test specificity using MFEprimer against "
                "human genome")
    logger.info("=" * 60)

    if human_genome:
        run_mfe_spec(output_dir + '/prefiltered_primers.fna',
                     DEFAULT_HG38_PATH, output_dir, 'hg38', tm_spec)
        df_human_spec = parse_MFEprimers_file(
            output_dir + '/mfe_spec_results_hg38.txt')
        if df_human_spec is not None and not df_human_spec.empty:
            df_human_spec = df_human_spec.drop(['Unique_hit',
                                                 'Max_hit_per_pair'],
                                                axis=1)
            df_human_spec.to_csv(
                output_dir + '/mfe_spec_results_hg38_summary.csv',
                sep='\t', index=False)
            print(df_human_spec)
    else:
        logger.warning("No Human Genome screening required")

    if runtype == 2:
        logger.info("=" * 60)
        logger.info("Step 4: Test primer compatibility for multiplex "
                    "and cocktail detection")
        logger.info("=" * 60)

        run_mfe_dimers(output_dir + '/prefiltered_primers.fna',
                       output_dir, 'prefiltered_primers')
        df_cocktail_dimers = parse_MFEprimers_dimer_file(
            output_dir + '/mfe_dimer_results_prefiltered_primers.txt')
        if df_cocktail_dimers is not None and \
           not df_cocktail_dimers.empty:
            df_cocktail_dimers.to_csv(
                output_dir +
                '/mfe_dimer_results_prefiltered_primers_summary.csv',
                sep='\t', index=False)
            print(df_cocktail_dimers)

    logger.info("=" * 60)
    logger.info("Final step: Summarising all results")
    logger.info("=" * 60)


def run_workflow_3(sequences: Dict[str, str], args: Namespace) -> None:
    """
    Run workflow 3 (UNIQUE TARGET) - unique region detection.

    Args:
        sequences: Dictionary of sequence headers and sequences
        args: Parsed command line arguments
    """
    logger.info("Workflow 3: Unique region detection")
    logger.info("This feature is coming soon")
    # TODO: Implement unique region detection
    # Will be implemented in Phase 3
