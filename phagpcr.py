#!/usr/bin/env python3
"""
Design and evaluate qPCR primers for phage detection in patient serum.

This pipeline uses Primer3 for primer design and MFEprimer for specificity
screening against custom databases, human genome, and primer dimer detection.
Useful for designing primers compatible with multiplex phage detection in
clinical samples.

Usage:
    python phagpcr.py -f input.fasta -o results/ -b blast_db.fna -hg
    python phagpcr.py -f phage.fasta -o output/ -r 2 -nb 5
"""
import subprocess
import os
import os.path

import pandas as pd
import primer3
from argparse import ArgumentParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from termcolor import colored

# Module constants
SEQUENCE_SUBREGION_MARGIN = 150  # bp to exclude from each end
MIN_SEQUENCE_FOR_SUBREGION = 500  # bp minimum for subregion design
TARGET_NAME_MAX_LENGTH = 55  # characters
DEFAULT_HG38_PATH = '/data/db/blastdb/hg38/\
GCA_000001405.29_GRCh38.p14_genomic.fna'
MFEPRIMER_BIN = 'bin/mfeprimer'  # Path to mfeprimer executable

def parse_args():
    """Parse command line arguments."""
    parser = ArgumentParser(
        description='Design primers using Primer3 and perform \
specificity screening')
    parser.add_argument('-f', '--fasta_file', 
                        type=str, 
                        help='Path to the input FASTA file')
    parser.add_argument('-o', '--output_dir', 
                        type=str, 
                        help='Directory where the results will be saved')
    parser.add_argument('-b', '--blast_db',
                        type=str,
                        default="data/blastdb/\
meta_ILL_ONT_20250521_combined.fna",
                        required=False,
                        help='BLAST database for primer screening \
(optional)')
    parser.add_argument('-hg', '--human_genome',
                        action='store_true',
                        help='Screen the primers against the human genome - \
version hg38 (optional)')
    parser.add_argument('-s', '--sequences_dir',
                        type=str, default=None, required=False,
                        help='Directory with sequences for screening \
(optional)')
    parser.add_argument('-r', '--runtype',
                        type=int, default=1, choices=[1, 2, 3],
                        help='Run type: 1=supervised, \
2=multiplex/cocktail, 3=unsupervised')
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
    return parser.parse_args()

def update_primer3(tm, primer_size, primer_nb, kit):
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
            [100, 150], [150, 175], [175, 200], [200, 225],
            [225, 250], [250, 275], [275, 300]
        ],
        'PRIMER_PRODUCT_OPT_SIZE': 150,
        # Advanced metrics
        'PRIMER_TM_FORMULA': 1,    # SantaLucia 1988
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_SALT_CORRECTIONS': 1,   # SantaLucia 1988
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_POLY_X': 4,
        'PRIMER_INTERNAL_MAX_POLY_X': 4,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        # Complementarity
        'PRIMER_MAX_SELF_ANY': 8.00,
        'PRIMER_MAX_SELF_END': 3.00,
        'PRIMER_MAX_END_STABILITY': 9.0,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        # Output
        'PRIMER_NUM_RETURN': primer_nb,
    }
    return primer3_params

def testing(args):
    """Validate input arguments."""
    # verify fasta file exists
    if not os.path.isfile(args.fasta_file):
        raise FileNotFoundError(f'Fasta file not found: \
{args.fasta_file}')

    # verify fasta file is readable
    if not os.access(args.fasta_file, os.R_OK):
        raise PermissionError(f'Cannot read fasta file: \
{args.fasta_file}')

def mkdir_outdir(output_dir):
    """Create output directory if it doesn't exist."""
    print(colored('1. Creating output directory:', "blue"), output_dir),
    if not os.path.exists('{}'.format(output_dir)):
        os.makedirs('{}'.format(output_dir))
    else:
        print(colored('==> Output directory already exists. \
Carrying on...\n', "yellow"))

def readfile(fasta_file):
    """Read FASTA file into a dictionary of sequences."""
    print(colored("2. Reading data in from ", "blue"), fasta_file)
    sequences = {}
    try:
        with open(fasta_file, 'r') as f:
            current_seq = []
            current_header = ''
            for line in f:
                if line.startswith('>'):
                    if current_header:
                        sequences[current_header] = ''.join(current_seq)
                    current_header = line[1:].strip().replace(" ", "_")
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_header:
                sequences[current_header] = ''.join(current_seq)

        if not sequences:
            raise ValueError(f"No sequences found in {fasta_file}")

        return sequences
    except IOError as e:
        raise IOError(f"Error reading file {fasta_file}: {str(e)}")    

def get_primers(header, seq, primer3_params):
    """Design primers for given sequence using Primer3."""
    print(colored("Designing primers for ", "blue"), header)
    # For long sequences, define subregion to avoid terminal regions
    if len(seq) >= MIN_SEQUENCE_FOR_SUBREGION:
        start = SEQUENCE_SUBREGION_MARGIN
        region_length = len(seq) - (start * 2)
    # Sequence too short for subregion, use full sequence
    else:
        start = 0
        region_length = len(seq)

    primers = []
    # Set the Primer3 sequence arguments
    primer3_input = {
        'SEQUENCE_ID': header,
        'SEQUENCE_TEMPLATE': seq,
        'SEQUENCE_INCLUDED_REGION': [start, region_length],
    }

    # Run Primer3 design - returns dictionary of primer results
    primer_results = primer3.bindings.designPrimers(
        primer3_input, primer3_params)
    return primer_results

def parse_primers(primer_results, header):
    """Parse Primer3 results into a structured DataFrame."""
    target_name = header[:TARGET_NAME_MAX_LENGTH].split(" ")[0].strip()
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
            'left_primer_sequence': primer_results[f'{left_key}_SEQUENCE'],
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
        
def df_to_fasta(df, output_dir, output_fasta):
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

def run_mfe_index(db):
    """Index database for MFEprimer analysis."""
    subprocess.run([MFEPRIMER_BIN, 'index', '-i', db], check=True)

def run_mfe(primers, db, output_dir, suffix, tm_spec):
    """Run MFEprimer full analysis on primers against database."""
    print(colored("Running MFEprimer full analysis for ", "blue") +
          primers + " ...")
    output_file = f"{output_dir}/mfe_results_{suffix}.txt"
    subprocess.run([
        MFEPRIMER_BIN, '-i', primers, '-d', db, '-o', output_file,
        '-S', '500', '-t', str(tm_spec), '--misMatch', '1',
        '--misStart', '2', '--misEnd', '9'
    ], check=True)
    
def run_mfe_dimers(primers, output_dir, suffix):
    """Run MFEprimer dimer analysis to detect primer interactions."""
    print(colored("Screening for dimers.", "blue"))
    print(colored("Running MFEprimer dimers...\n", "blue"))
    output_file = f"{output_dir}/mfe_dimer_results_{suffix}.txt"
    subprocess.run([
        MFEPRIMER_BIN, 'dimer', '-i', primers,
        '--dg', '-10', '-o', output_file
    ], check=True)

def run_mfe_spec(primers, db, output_dir, suffix, tm_spec):
    """Run MFEprimer specificity check against database."""
    print(colored("Screening for specificity against Blast database: ",
                  "blue") + str(db))
    print(colored("Running MFEprimer specificity...\n", "blue"))
    output_file = f"{output_dir}/mfe_spec_results_{suffix}.txt"
    subprocess.run([
        MFEPRIMER_BIN, 'spec', '-i', primers, '-d', db,
        '-o', output_file, '-S', '500', '-t', str(tm_spec)
    ], check=True)   

def parse_MFEprimers_file(filename):
    # read the file and split it into lines
    try:
        with open(filename) as f:
            lines = f.read().splitlines()
        
        # find the line that indicates the number of potential amplicons
        start_line_matches = [i for i, line in enumerate(lines) if line.startswith("Descriptions of [")]
        if not start_line_matches:
            print(colored(f"==> Unexpected file format in {filename}: 'Descriptions of [' line not found", "yellow"))
            return pd.DataFrame()
        
        start_line = start_line_matches[0]

        # count the number of amplicons
        amplicons_count = int(lines[start_line].split(" ")[3])
        
        if amplicons_count >= 1:
            # extract the details of each amplicon
            amplicon_details = []
            for i in range(amplicons_count):
                amplicon_id = "Amp " + str(i+1)
                amp_line_matches = [j for j, line in enumerate(lines) if line.startswith(amplicon_id)]
                if not amp_line_matches:
                    print(colored(f"==> Expected amplicon {amplicon_id} not found in {filename}", "yellow"))
                    continue
                amp_line = amp_line_matches[0]
                amplicon_fp = lines[amp_line].split(" ")[2].strip()
                amplicon_rp = lines[amp_line].split(" ")[4].strip()
                amplicon_hit = lines[amp_line].split("==>")[1].strip()
                amplicon_size = int(lines[amp_line + 2].split("=")[1].split()[0])
                amplicon_gc = round(float(lines[amp_line + 2].split("=")[2].strip("%")), 2)
                amplicon_fpTm = round(float(lines[amp_line + 3].split("=")[1].strip("°C, Delta G ")), 2)
                amplicon_rpTm = round(float(lines[amp_line + 4].split("=")[1].strip("°C, Delta G ")), 2)
                amplicon_fpDg = round(float(lines[amp_line + 3].split("=")[2].strip("kcal/mol")), 2)
                amplicon_rpDg = round(float(lines[amp_line + 4].split("=")[2].strip("kcal/mol")), 2)
                amplicon_binding_sites = lines[amp_line + 5].split()[2] + "-" + lines[amp_line + 5].split()[2]
                amplicon_details.append((amplicon_id, amplicon_fp, amplicon_rp, amplicon_hit, amplicon_size, amplicon_gc, amplicon_fpTm, amplicon_rpTm, amplicon_fpDg, amplicon_rpDg, amplicon_binding_sites))    
            # create a Pandas DataFrame with the amplicon details
            if amplicon_details:
                df = pd.DataFrame(amplicon_details, columns=["Amplicon ID", "Fp", "Rp", "Hit ID", "Size", "GC", "Fp Tm", "Rp Tm", "Fp Dg", "Rp Dg", "Binding Sites"])
                df['Header'] = df.apply(lambda x: x['Fp'].rsplit('_', 1)[0], axis=1)
                df['Max_hit_per_pair'] = df.apply(lambda x: amplicons_count, axis=1)
                df['Unique_hit'] = df.apply(lambda x: 'No' if amplicons_count >1 else 'Yes', axis=1)
                #print(df, "\n")
                return df
            else:
                return pd.DataFrame()
        else:
            print(colored("==> No primer pairs were found to bind and produce amplicons (within the limit specified - default 500bp.)", "yellow"))
            return pd.DataFrame()
    except Exception as e:
        print(colored(f"==> Error parsing file {filename}: {str(e)}", "red"))
        return pd.DataFrame()

def parse_MFEprimers_dimer_file(filename):
    """Parse MFEprimer dimer output file into DataFrame."""
    try:
        with open(filename) as f:
            lines = f.read().splitlines()

        # find the line that indicates the number of potential dimers
        start_line_matches = [i for i, line in enumerate(lines)
                              if line.startswith("Dimer List ")]
        if not start_line_matches:
            print(colored(f"==> Unexpected file format in {filename}: \
'Dimer List' line not found", "yellow"))
            return pd.DataFrame()

        start_line = start_line_matches[0]

        # count the number of dimers
        dimer_count = int(lines[start_line].split("(")[1].split(")")[0])

        if dimer_count >= 1:
            # extract the details of each dimer
            dimer_details = []
            for i in range(dimer_count):
                dim_id = "Dimer " + str(i+1)
                dim_line_matches = [j for j, line in enumerate(lines)
                                    if line.startswith(dim_id)]
                if not dim_line_matches:
                    print(colored(f"==> Expected dimer {dim_id} \
not found in {filename}", "yellow"))
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
            print(colored("==> No primer pairs were found to form \
dimers (within the limit specified)", "yellow"))
            return pd.DataFrame()
    except Exception as e:
        print(colored(f"==> Error parsing file {filename}: {str(e)}",
                      "red"))
        return pd.DataFrame()

if __name__ == '__main__':
    
    args = parse_args()
    
    fasta_file = args.fasta_file
    output_dir = args.output_dir
    blast_db = args.blast_db
    sequences_dir = args.sequences_dir
    human_genome = args.human_genome
    runtype = args.runtype
    primer_nb = args.primer_nb
    tm = args.tm_optimal
    primer_size = args.size_optimal
    kit = args.kit
    tm_spec = args.tm_specificity

    # Testing arguments
    testing(args)
    
    # Creare output directory  
    mkdir_outdir(output_dir)
    mkdir_outdir(output_dir + "/indiv_results")
    
    # Read fasta_file in
    sequences = readfile(fasta_file)
    print("Number of sequences submitted: ",len(sequences), "\n")
    
    # Primer selection and testing for runtype 1 and 2
    if runtype == 1 or runtype == 2:
        print("--------------------------------------------------------")
        print("Step 1 - select all primers against target using Primer3")
        print("--------------------------------------------------------")
        primers=[]
        for header, seq in sequences.items():
            primer3_params = update_primer3(tm, primer_size, primer_nb, kit)
            primer_results = get_primers(header, seq, primer3_params)         
            prefiltered_primers_df = parse_primers(primer_results, header)
        
        # create prefiltered primers file
        prefiltered_primers_df.to_csv(output_dir + '/prefiltered_primers_file.csv', index=False, sep ='\t')
        print(colored('\nPrefiltered primers stored in: ', "yellow"), output_dir + '/prefiltered_primers_file.csv')
        
        # create multifasta file with all predicted priners
        df_to_fasta(prefiltered_primers_df, output_dir, 'prefiltered_primers')
        print(colored('Prefiltered primers in fasta format stored in: ', "yellow"), output_dir + '/prefiltered_primers.fna')
 
        print("\n--------------------------------------------------------")               
        print("Step 2 - test specificity of primers using MFEprimer against custom database")
        print("--------------------------------------------------------")
        if blast_db != None:
            ## by default, mfeprimer_index will check if the indexing has already been performed before proceeding
            run_mfe_index(blast_db)
            print(colored("Screening for all features and specificity against Blast database: ", "blue") + str(blast_db) + "\n")
            # run_mfe(output_dir + '/prefiltered_primers.fna',blast_db,output_dir,'blastdb', tm_spec)
            # testing running MFE on separate fasta files         
            grouped = prefiltered_primers_df.groupby('left_primer_name')
            df_local_screening = pd.DataFrame()
            for name, group in grouped:
                df_to_fasta(grouped.get_group(name), output_dir, name[:-1])
                run_mfe(output_dir + '/' + str(name[:-1]) + '.fna', blast_db, output_dir + '/indiv_results', str(name[:-1]), tm_spec)
                local_screening = parse_MFEprimers_file(output_dir + '/indiv_results/mfe_results_' + str(name[:-1]) + '.txt')
                if local_screening is not None and not local_screening.empty:
                    df_local_screening = pd.concat([df_local_screening, local_screening], ignore_index=True)
                os.remove(output_dir + '/' + str(name[:-1]) + '.fna')
            print(df_local_screening)
            if not df_local_screening.empty:
                df_local_screening.to_csv(output_dir + '/mfe_screening_results_against_local_database_summary.csv', sep='\t', index=False)
                if 'Fp' in df_local_screening.columns:
                    print(df_local_screening.groupby(df_local_screening['Fp'].str[:-1])['Max_hit_per_pair'].unique())
                else:
                    print(colored("==> No 'Fp' column found in screening results", "yellow"))
            else:
                print(colored("==> No screening results to save. DataFrame is empty.", "yellow"))
        else:
            print(colored("No Blast database provided", "yellow"))
            
        print("\n--------------------------------------------------------")               
        print("Step 3 - test specificity of primers using MFEprimer against human genome")
        print("--------------------------------------------------------")    
        if human_genome:
            # Run MFEprimer specificity against human genome
            run_mfe_spec(output_dir + '/prefiltered_primers.fna',
                         DEFAULT_HG38_PATH, output_dir, 'hg38', tm_spec)
            df_human_spec = parse_MFEprimers_file(output_dir + '/mfe_spec_results_hg38.txt')
            if df_human_spec is not None and not df_human_spec.empty:
                df_human_spec = df_human_spec.drop(['Unique_hit', 'Max_hit_per_pair'], axis=1)
                df_human_spec.to_csv(output_dir + '/mfe_spec_results_hg38_summary.csv', sep='\t', index=False)
                print(df_human_spec)
        else:
            print(colored("No Human Genome screening required", "yellow"))
 
        if runtype == 2:
            print("\n--------------------------------------------------------")               
            print("Step 4 - test primer compatibility for multiplex and cocktail detection")
            print("--------------------------------------------------------")               
            # run MFEprimer dimers on combined file 
            run_mfe_dimers(output_dir + '/prefiltered_primers.fna', output_dir, 'prefiltered_primers')
            df_cocktail_dimers = parse_MFEprimers_dimer_file(output_dir + '/mfe_dimer_results_prefiltered_primers.txt')
            if df_cocktail_dimers is not None and not df_cocktail_dimers.empty:
                df_cocktail_dimers.to_csv(output_dir + '/mfe_dimer_results_prefiltered_primers_summary.csv', sep='\t', index=False)
                print(df_cocktail_dimers)
           
        print("\n--------------------------------------------------------")
        print("Final step - summarising all results")
        print("--------------------------------------------------------")       
    
    # Type 3: unsupervised design providing phage genome
    elif runtype == 3:
        # for header, seq in sequences.items():
            # search for primase gene
            # if present
            #   primase = ...
            #   get_primers(header_primase, seq_primase, primer3_params)
            # search for polymerase gene
            # if present
            #   polymerase = ...
            #   get_primers(header_polymerase, seq_polymerase, primer3_params)           
        print("new feature to come")
        
