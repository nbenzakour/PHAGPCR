"""Constants used throughout PHAGPCR package."""

# Sequence analysis parameters
SEQUENCE_SUBREGION_MARGIN = 150  # bp to exclude from each end
MIN_SEQUENCE_FOR_SUBREGION = 500  # bp minimum for subregion design
TARGET_NAME_MAX_LENGTH = 55  # characters

# Default file paths
DEFAULT_HG38_PATH = ('/data/db/blastdb/hg38/'
                     'GCA_000001405.29_GRCh38.p14_genomic.fna')
DEFAULT_BLAST_DB = 'data/blastdb/meta_ILL_ONT_20250521_combined.fna'
MFEPRIMER_BIN = 'bin/mfeprimer'  # Path to mfeprimer executable
