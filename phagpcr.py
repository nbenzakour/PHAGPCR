import sys
import os
import os.path
import pandas as pd
import primer3
from argparse import ArgumentParser

primer3_params = {
    # Basic parameters
    'PRIMER_OPT_SIZE' : 20, 
    'PRIMER_MIN_SIZE' : 19,
    'PRIMER_MAX_SIZE' : 23,
    'PRIMER_PICK_INTERNAL_OLIGO' : 0,
    'PRIMER_OPT_TM' : 60.0,
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
    'PRIMER_NUM_RETURN' : 4,
}

def parse_args():
    parser = ArgumentParser(description='Design primers using Primer3 and perform BLAST search on the designed primers')
    parser.add_argument('-f', '--fasta_file', 
                        type=str, 
                        help='Path to the input FASTA file')
    parser.add_argument('-o', '--output_dir', 
                        type=str, 
                        help='Directory where the results will be saved')
    parser.add_argument('-b', '--blast_db', 
                        type=str, default=None, required=False, 
                        help='BLAST database for searching the primers (optional)')
    parser.add_argument('-H', '--human_genome', 
                        type=bool, default=False, required=False, 
                        help='Screen the primers against the human genome - verison hg38 (True or False)(optional)')
    parser.add_argument('-s', '--sequences_dir', 
                        type=str, default=None, required=False, 
                        help='Directory containing sequences for screening against the primers (optional)')
    parser.add_argument('-r', '--runtype', 
                        type=int, default=1, choices=[1, 2, 3], 
                        help='Run type for primer design: 1 = supervised design providing target, 2 = unsupervised design providing phage genome, 3 = cocktail detection with multiplex primers')
    return parser.parse_args()

def testing(args):
    # verify fasta file exists
    if not os.path.isfile(args.fasta_file):
        raise Exception('Fasta file was not found')

def mkdir_outdir(output_dir):
    # Creating the output directory if not present
    print('1. Creating output directory:', output_dir)
    if not os.path.exists('{}'.format(output_dir)):
        os.makedirs('{}'.format(output_dir))
    else:
        print('==> Output directory already exists. Carrying on...')

def readfile(fasta_file):
    # Read the FASTA file into a dictionary
    print("2. Reading data in from ", fasta_file)
    sequences = {}
    with open(fasta_file, 'r') as f:
        current_seq = ''
        current_header = ''
        for line in f:
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = current_seq
                current_header = line[1:].strip()
                current_seq = ''
            else:
                current_seq += line.strip()
        sequences[current_header] = current_seq
    return sequences     

def get_primers(header, seq, primer3_params):
    print("3. Designing primers for ", header)
    primers = []
    # Set the Primer3 sequence arguments
    primer3_input = {
        'SEQUENCE_ID': header,
        'SEQUENCE_TEMPLATE': seq,
    }
    
    # Run the Primer3 design process - return a dictionary of primer results
    primer_results = primer3.bindings.designPrimers(primer3_input, primer3_params)
    return primer_results

def parse_primers(primer_results, header):
    target_name = header[:30].split(" ")[0].strip()
    # Parse the primers
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
    return(primers)

#def run_mfe(primer_fasta_file, blast_db, human_genome)

if __name__ == '__main__':
    
    args = parse_args()
    
    fasta_file = args.fasta_file
    output_dir = args.output_dir
    blast_db = args.blast_db
    sequences_dir = args.sequences_dir
    human_genome = args.human_genome
    runtype = args.runtype

    # Testing arguments
    testing(args)
    
    # Creare output directory  
    mkdir_outdir(output_dir)
    
    # Read fasta_file in
    sequences = readfile(fasta_file)
    print("Number of sequences submitted: ",len(sequences))
    
    # Type 1: supervised design providing target
    if runtype == 1:
        print("Run type for primer design: 1 = supervised design providing target")
        print("--------------------------------------------------------")
        print("Step 1 - select all primers against target using Primer3")
        print("--------------------------------------------------------")
        primers=[]
        for header, seq in sequences.items():
            primer_results = get_primers(header, seq, primer3_params)         
            primer_summary = parse_primers(primer_results, header)
        all_primers = pd.DataFrame(primer_summary)
        # create prefiltered primers file
        all_primers.to_csv(output_dir + '/prefiltered_primers_file.csv', index=False)
        print('Prefiltered primers stored in:', output_dir + '/prefiltered_primers_file.csv')
 
        print("--------------------------------------------------------")               
        print("Step 2 - test specificity of primers using MFEprimer")
        print("--------------------------------------------------------")
        if blast_db == True:
            print("Screening against Blast database:" + blast_db)
        else:
            print("No Blast database provided")    
        if human_genome == "yes":
            print("Screening against Human Genome version Hg38")
        else:
            print("No Human Genome screening required")
        if (blast_db == True) | (human_genome == "yes"):
            ## run mfe
            print("run MFE")
        else:
            print("Specificity of primers using MFEprimer not required")
            
        print("--------------------------------------------------------")
        print("Step 3 - summarising all results")
        print("--------------------------------------------------------")       
        
    
    # Type 2: unsupervised design providing phage genome
    elif runtype == 2:
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
        
    # Type 3: cocktail detection with multiplex primers')
    else:
        print("other new feature to come")