"""MFEprimer utilities for primer specificity screening."""

import logging
import subprocess
from typing import Optional

import pandas as pd

from phagpcr_package.constants import MFEPRIMER_BIN

logger = logging.getLogger(__name__)


def run_mfe_index(db: str) -> None:
    """Index database for MFEprimer analysis."""
    subprocess.run([MFEPRIMER_BIN, 'index', '-i', db], check=True)


def run_mfe(primers: str, db: str, output_dir: str,
            suffix: str, tm_spec: int) -> None:
    """Run MFEprimer full analysis on primers against database."""
    logger.info(f"Running MFEprimer full analysis for {primers} ...")
    output_file = f"{output_dir}/mfe_results_{suffix}.txt"
    subprocess.run([
        MFEPRIMER_BIN, '-i', primers, '-d', db, '-o', output_file,
        '-S', '500', '-t', str(tm_spec), '--misMatch', '1',
        '--misStart', '2', '--misEnd', '9'
    ], check=True)


def run_mfe_dimers(primers: str, output_dir: str, suffix: str) -> None:
    """Run MFEprimer dimer analysis to detect primer interactions."""
    logger.info("Screening for dimers.")
    logger.info("Running MFEprimer dimers...")
    output_file = f"{output_dir}/mfe_dimer_results_{suffix}.txt"
    subprocess.run([
        MFEPRIMER_BIN, 'dimer', '-i', primers,
        '--dg', '-10', '-o', output_file
    ], check=True)


def run_mfe_spec(primers: str, db: str, output_dir: str,
                 suffix: str, tm_spec: int) -> None:
    """Run MFEprimer specificity check against database."""
    logger.info(f"Screening for specificity against Blast database: "
                f"{db}")
    logger.info("Running MFEprimer specificity...")
    output_file = f"{output_dir}/mfe_spec_results_{suffix}.txt"
    subprocess.run([
        MFEPRIMER_BIN, 'spec', '-i', primers, '-d', db,
        '-o', output_file, '-S', '500', '-t', str(tm_spec)
    ], check=True)


def parse_MFEprimers_file(filename: str) -> Optional[pd.DataFrame]:
    """Parse MFEprimer output file into DataFrame."""
    try:
        with open(filename) as f:
            lines = f.read().splitlines()

        # Find line indicating number of potential amplicons
        start_line_matches = [i for i, line in enumerate(lines)
                              if line.startswith("Descriptions of [")]
        if not start_line_matches:
            logger.warning(f"==> Unexpected file format in {filename}: "
                           f"'Descriptions of [' line not found")
            return pd.DataFrame()

        start_line = start_line_matches[0]

        # Count number of amplicons
        amplicons_count = int(lines[start_line].split(" ")[3])

        if amplicons_count >= 1:
            # Extract details of each amplicon
            amplicon_details = []
            for i in range(amplicons_count):
                amplicon_id = "Amp " + str(i+1)
                amp_line_matches = [j for j, line in enumerate(lines)
                                    if line.startswith(amplicon_id)]
                if not amp_line_matches:
                    logger.warning(f"==> Expected amplicon {amplicon_id} "
                                   f"not found in {filename}")
                    continue
                amp_line = amp_line_matches[0]
                amplicon_fp = lines[amp_line].split(" ")[2].strip()
                amplicon_rp = lines[amp_line].split(" ")[4].strip()
                amplicon_hit = lines[amp_line].split("==>")[1].strip()
                amplicon_size = int(lines[amp_line + 2].split("=")[1]
                                    .split()[0])
                amplicon_gc = round(float(lines[amp_line + 2]
                                          .split("=")[2].strip("%")), 2)
                amplicon_fpTm = round(float(lines[amp_line + 3]
                                            .split("=")[1]
                                            .strip("°C, Delta G ")), 2)
                amplicon_rpTm = round(float(lines[amp_line + 4]
                                            .split("=")[1]
                                            .strip("°C, Delta G ")), 2)
                amplicon_fpDg = round(float(lines[amp_line + 3]
                                            .split("=")[2]
                                            .strip("kcal/mol")), 2)
                amplicon_rpDg = round(float(lines[amp_line + 4]
                                            .split("=")[2]
                                            .strip("kcal/mol")), 2)
                amplicon_binding_sites = (lines[amp_line + 5].split()[2] +
                                          "-" +
                                          lines[amp_line + 5].split()[2])
                amplicon_details.append((amplicon_id, amplicon_fp,
                                         amplicon_rp, amplicon_hit,
                                         amplicon_size, amplicon_gc,
                                         amplicon_fpTm, amplicon_rpTm,
                                         amplicon_fpDg, amplicon_rpDg,
                                         amplicon_binding_sites))
            # Create Pandas DataFrame with amplicon details
            if amplicon_details:
                df = pd.DataFrame(amplicon_details,
                                  columns=["Amplicon ID", "Fp", "Rp",
                                           "Hit ID", "Size", "GC",
                                           "Fp Tm", "Rp Tm", "Fp Dg",
                                           "Rp Dg", "Binding Sites"])
                df['Header'] = df.apply(lambda x: x['Fp'].rsplit('_', 1)[0],
                                        axis=1)
                df['Max_hit_per_pair'] = df.apply(lambda x: amplicons_count,
                                                   axis=1)
                df['Unique_hit'] = df.apply(
                    lambda x: 'No' if amplicons_count > 1 else 'Yes',
                    axis=1)
                return df
            else:
                return pd.DataFrame()
        else:
            logger.warning("==> No primer pairs were found to bind and "
                           "produce amplicons (within the limit specified "
                           "- default 500bp.)")
            return pd.DataFrame()
    except Exception as e:
        logger.error(f"==> Error parsing file {filename}: {str(e)}")
        return pd.DataFrame()


def parse_MFEprimers_dimer_file(filename: str) -> Optional[pd.DataFrame]:
    """Parse MFEprimer dimer output file into DataFrame."""
    try:
        with open(filename) as f:
            lines = f.read().splitlines()

        # Find line indicating number of potential dimers
        start_line_matches = [i for i, line in enumerate(lines)
                              if line.startswith("Dimer List ")]
        if not start_line_matches:
            logger.warning(f"==> Unexpected file format in {filename}: "
                           f"'Dimer List' line not found")
            return pd.DataFrame()

        start_line = start_line_matches[0]

        # Count number of dimers
        dimer_count = int(lines[start_line].split("(")[1].split(")")[0])

        if dimer_count >= 1:
            # Extract details of each dimer
            dimer_details = []
            for i in range(dimer_count):
                dim_id = "Dimer " + str(i+1)
                dim_line_matches = [j for j, line in enumerate(lines)
                                    if line.startswith(dim_id)]
                if not dim_line_matches:
                    logger.warning(f"==> Expected dimer {dim_id} "
                                   f"not found in {filename}")
                    continue
                dim_line = dim_line_matches[0]
                dim_fp = lines[dim_line].split(" ")[2].strip()
                dim_rp = lines[dim_line].split(" ")[4].strip()
                dim_score = lines[dim_line + 2].split(" ")[1].strip(",")
                dim_Dg = round(float(lines[dim_line + 2].split("=")[1]
                                     .strip("kcal/mol")), 2)
                dimer_details.append((dim_id, dim_fp, dim_rp,
                                      dim_score, dim_Dg))

            if dimer_details:
                df = pd.DataFrame(dimer_details,
                                  columns=["Dimer ID", "Fp", "Rp",
                                           "Score", "Dg"])
                df['Problematic'] = df.apply(
                    lambda x: 'Yes' if x['Fp'].rsplit('_', 1)[0] !=
                    x['Rp'].rsplit('_', 1)[0] else 'No', axis=1)
                return df
            else:
                return pd.DataFrame()
        else:
            logger.warning("==> No primer pairs were found to form "
                           "dimers (within the limit specified)")
            return pd.DataFrame()
    except Exception as e:
        logger.error(f"==> Error parsing file {filename}: {str(e)}")
        return pd.DataFrame()
