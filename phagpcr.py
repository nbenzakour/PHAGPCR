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

def parse_args():
    parser = ArgumentParser(description='Design primers using Primer3 and perform BLAST search on the designed primers')
    parser.add_argument('-f', '--fasta_file', 
                        type=str, 
                        help='Path to the input FASTA file')
    parser.add_argument('-o', '--output_dir', 
                        type=str, 
                        help='Directory where the results will be saved')
    parser.add_argument('-b', '--blast_db', 
                        type=str, default="~/REPOS/PHAGPCR/test/phiau_all_phages_20230418_nr.fna", required=False, 
                        help='BLAST database for searching the primers (optional)')
    parser.add_argument('-hg', '--human_genome',
                        action='store_true',
                        help='Screen the primers against the human genome - \
version hg38 (optional)')
    parser.add_argument('-s', '--sequences_dir', 
                        type=str, default=None, required=False, 
                        help='Directory containing sequences for screening against the primers (optional)')
    parser.add_argument('-r', '--runtype', 
                        type=int, default=1, choices=[1, 2, 3], 
                        help='Run type for primer design: 1 = supervised design providing target, 2 = cocktail detection with multiplex primers, 3 = unsupervised design providing phage genome')
    parser.add_argument('-nb', '--primer_nb', 
                        type=int, default=1, required=False, 
                        help='Nb of primer pairs to design per target sequence')
    parser.add_argument('-tm', '--tm_optimal', 
                        type=int, default=60, required=False, 
                        help='Optimal Tm for primers')
    parser.add_argument('-sz', '--size_optimal', 
                        type=int, default=20, required=False, 
                        help='Optimal size for primers')
    parser.add_argument('-k', '--kit', 
                        type=str, default="", choices=["Quantinova", "Invitrogen"], required=False, 
                        help='qPCR kit used to inform Tm specification')
    parser.add_argument('-tm2', '--tm_specificity', 
                        type=int, default=40, required=False, 
                        help='Minimum Tm to report primers-template matches screened with MFEprimer')
    return parser.parse_args()

def update_primer3(tm, primer_size, primer_nb, kit):
    primer3_params = {
        # Basic parameters
        'PRIMER_OPT_SIZE' : primer_size, 
        'PRIMER_MIN_SIZE' : 18,
        'PRIMER_MAX_SIZE' : 23,
        'PRIMER_PICK_INTERNAL_OLIGO' : 0,
        'PRIMER_OPT_TM' : float(tm),
        'PRIMER_MIN_TM' : 59.0,
        'PRIMER_MAX_TM' : 63.0,
        'PRIMER_PAIR_MAX_DIFF_TM' : 100.0,
        'PRIMER_MIN_GC' : 40.0, 
        'PRIMER_MAX_GC' : 60.0,
        # Product
        'PRIMER_PRODUCT_SIZE_RANGE' : [[100,150],[150,175],[175,200],[200,225],[225,250],[250,275],[275,300]],
        'PRIMER_PRODUCT_OPT_SIZE' : 150,
        # Advanced metrics
        'PRIMER_TM_FORMULA' : 1,    # SantaLucia 1988
        'PRIMER_SALT_MONOVALENT' : 50.0,
        'PRIMER_SALT_CORRECTIONS' : 1,   # SantaLucia 1988
        'PRIMER_DNA_CONC' : 50.0,
        'PRIMER_MAX_POLY_X' : 4,
        'PRIMER_INTERNAL_MAX_POLY_X' : 4,
        'PRIMER_MAX_NS_ACCEPTED' : 0,
        # Complementarity
        'PRIMER_MAX_SELF_ANY' : 8.00,
        'PRIMER_MAX_SELF_END' : 3.00,
        'PRIMER_MAX_END_STABILITY' : 9.0,
        'PRIMER_PAIR_MAX_COMPL_ANY' : 12,
        'PRIMER_PAIR_MAX_COMPL_END' : 8,
        # Output
        'PRIMER_NUM_RETURN' : primer_nb,
    }
    return primer3_params

def testing(args):
    # verify fasta file exists
    if not os.path.isfile(args.fasta_file):
        raise Exception('Fasta file was not found')

def mkdir_outdir(output_dir):
    # Creating the output directory if not present
    print(colored('1. Creating output directory:', "blue"), output_dir), 
    if not os.path.exists('{}'.format(output_dir)):
        os.makedirs('{}'.format(output_dir))
    else:
        print(colored('==> Output directory already exists. Carrying on...\n',"yellow"))

def readfile(fasta_file):
    # Read the FASTA file into a dictionary
    print(colored("2. Reading data in from ", "blue"), fasta_file)
    sequences = {}
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
        sequences[current_header] = ''.join(current_seq)
    return sequences    

def get_primers(header, seq, primer3_params):
    print(colored("Designing primers for ", "blue"), header)
    # if sequence is equal or greater to 500 bp, define a subregion to design primers (=minimum 200bp) 
    if len(seq) >= 500:
        start = 150
        region_length = len(seq) - (start*2)
    # sequence is too short to define a subregion, use full sequence instead
    else:
        start = 0
        region_length = len(seq)
        
    primers = []
    # Set the Primer3 sequence arguments
    primer3_input = {
        'SEQUENCE_ID': header,
        'SEQUENCE_TEMPLATE': seq,
        'SEQUENCE_INCLUDED_REGION' : [start,region_length],
    }
    
    # Run the Primer3 design process - return a dictionary of primer results
    primer_results = primer3.bindings.designPrimers(primer3_input, primer3_params)
    return primer_results

def parse_primers(primer_results, header):
    target_name = header[:55].split(" ")[0].strip()
    # Parse the primers
    primers = []
    i = 0
    while i < primer_results['PRIMER_PAIR_NUM_RETURNED']:       
        primers.append({'header': header, 
                        'left_primer_name': target_name + "_" + str(i) + "F",
                        'left_primer_sequence': primer_results['PRIMER_LEFT_' + str(i) + '_SEQUENCE'], 
                        'right_primer_name': target_name + "_" + str(i) + "R", 
                        'right_primer_sequence': primer_results['PRIMER_RIGHT_' + str(i) + '_SEQUENCE'],
                        'left position': primer_results['PRIMER_LEFT_' + str(i)], 
                        'right position': primer_results['PRIMER_RIGHT_' + str(i)], 
                        'product size': primer_results['PRIMER_PAIR_' + str(i) + '_PRODUCT_SIZE'],
                        'left_Tm' : primer_results['PRIMER_LEFT_' + str(i) + '_TM'],
                        'right_Tm' : primer_results['PRIMER_RIGHT_' + str(i) + '_TM'],
                        'left_GC' : primer_results['PRIMER_LEFT_' + str(i) + '_GC_PERCENT'],
                        'right_GC' : primer_results['PRIMER_RIGHT_' + str(i) + '_GC_PERCENT'],
                        'left_SELF_ANY_TH' : primer_results['PRIMER_LEFT_' + str(i) + '_SELF_ANY_TH'],
                        'right_SELF_ANY_TH' : primer_results['PRIMER_RIGHT_' + str(i) + '_SELF_ANY_TH'],
                        'left_SELF_END_TH' : primer_results['PRIMER_LEFT_' + str(i) + '_SELF_END_TH'],
                        'right_SELF_END_TH' : primer_results['PRIMER_RIGHT_' + str(i) + '_SELF_END_TH'],
                        'left_HAIRPIN_TH' : primer_results['PRIMER_LEFT_' + str(i) + '_HAIRPIN_TH'],
                        'right_HAIRPIN_TH' : primer_results['PRIMER_RIGHT_' + str(i) + '_HAIRPIN_TH'],
                        'left_END_STABILITY' : primer_results['PRIMER_LEFT_' + str(i) + '_END_STABILITY'],
                        'right_END_STABILITY' : primer_results['PRIMER_RIGHT_' + str(i) + '_END_STABILITY'],
                        'left_PENALTY' : primer_results['PRIMER_LEFT_' + str(i) + '_PENALTY'],
                        'right_PENALTY' : primer_results['PRIMER_RIGHT_' + str(i) + '_PENALTY'],
                        'pair_PENALTY' : primer_results['PRIMER_PAIR_' + str(i) + '_PENALTY'],
                        'pair_COMPL_ANY_TH' : primer_results['PRIMER_PAIR_' + str(i) + '_COMPL_ANY_TH'],
                        'pair_COMPL_END_TH' : primer_results['PRIMER_PAIR_' + str(i) + '_COMPL_END_TH'],})
        i = i+1
    primers_df = pd.DataFrame(primers)
    primers_df = primers_df.round(2)
    return(primers_df)
        
def df_to_fasta(df, output_dir, output_fasta):
    records = [SeqRecord(Seq(row['left_primer_sequence']), id=row['left_primer_name'], description="") for index, row in df.iterrows()] + \
              [SeqRecord(Seq(row['right_primer_sequence']), id=row['right_primer_name'], description="") for index, row in df.iterrows()]
    with open(f"{output_dir}/{output_fasta}.fna", "w") as handle:
        SeqIO.write(records, handle, 'fasta')

def run_mfe_index(db):
    """Index database for MFEprimer analysis."""
    subprocess.run(['mfeprimer', 'index', '-i', db], check=True)

def run_mfe(primers, db, output_dir, suffix, tm_spec):
    """Run MFEprimer full analysis on primers against database."""
    print(colored("Running MFEprimer full analysis for ", "blue") +
          primers + " ...")
    output_file = f"{output_dir}/mfe_results_{suffix}.txt"
    subprocess.run([
        'mfeprimer', '-i', primers, '-d', db, '-o', output_file,
        '-S', '500', '-t', str(tm_spec), '--misMatch', '1',
        '--misStart', '2', '--misEnd', '9'
    ], check=True)
    
def run_mfe_dimers(primers, output_dir, suffix):
    """Run MFEprimer dimer analysis to detect primer interactions."""
    print(colored("Screening for dimers.", "blue"))
    print(colored("Running MFEprimer dimers...\n", "blue"))
    output_file = f"{output_dir}/mfe_dimer_results_{suffix}.txt"
    subprocess.run([
        'mfeprimer', 'dimer', '-i', primers,
        '--dg', '-10', '-o', output_file
    ], check=True)

def run_mfe_spec(primers, db, output_dir, suffix, tm_spec):
    """Run MFEprimer specificity check against database."""
    print(colored("Screening for specificity against Blast database: ",
                  "blue") + str(db))
    print(colored("Running MFEprimer specificity...\n", "blue"))
    output_file = f"{output_dir}/mfe_spec_results_{suffix}.txt"
    subprocess.run([
        'mfeprimer', 'spec', '-i', primers, '-d', db,
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
    # read the file and split it into lines
    with open(filename) as f:
        lines = f.read().splitlines()
    
    # find the line that indicates the number of potential dimers
    start_line = [i for i, line in enumerate(lines) if line.startswith("Dimer List ")][0]

    # count the number of dimers
    dimer_count = int(lines[start_line].split("(")[1].split(")")[0])
    
    if dimer_count >= 1:
        # extract the details of each dimer
        dimer_details = []
        for i in range(dimer_count):
            dim_id = "Dimer " + str(i+1)
            dim_line = [j for j, line in enumerate(lines) if line.startswith(dim_id)][0]
            dim_fp = lines[dim_line].split(" ")[2].strip()
            dim_rp = lines[dim_line].split(" ")[4].strip()
            dim_score = lines[dim_line + 2].split(" ")[1].strip(",")
            dim_Dg = round(float(lines[dim_line + 2].split("=")[1].strip("kcal/mol")), 2)
            dimer_details.append((dim_id, dim_fp, dim_rp, dim_score, dim_Dg))    
        # create a Pandas DataFrame with the amplicon details
        df = pd.DataFrame(dimer_details, columns=["Dimer ID", "Fp", "Rp", "Score", "Dg"])
        df['Problematic'] = df.apply(lambda x: 'Yes' if x['Fp'].rsplit('_', 1)[0] != x['Rp'].rsplit('_', 1)[0] else 'No', axis=1)
        #print(df, "\n")
        return df
    else:
        print(colored("==> No primer pairs were found to form dimers (within the limit specified)", "yellow"))

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
            # run MFEprimer specificity on combined file against human genome
            run_mfe_spec(output_dir + '/prefiltered_primers.fna','/data/db/blastdb/hg38/GCA_000001405.29_GRCh38.p14_genomic.fna',output_dir,'hg38', tm_spec)
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
        
