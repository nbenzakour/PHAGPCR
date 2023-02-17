import subprocess
import sys
import os
import os.path
import pandas as pd
import primer3
from argparse import ArgumentParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

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
    parser.add_argument('-hg', '--human_genome', 
                        type=bool, default=False, required=False, 
                        help='Screen the primers against the human genome - version hg38 (True or False)(optional)')
    parser.add_argument('-s', '--sequences_dir', 
                        type=str, default=None, required=False, 
                        help='Directory containing sequences for screening against the primers (optional)')
    parser.add_argument('-r', '--runtype', 
                        type=int, default=1, choices=[1, 2, 3], 
                        help='Run type for primer design: 1 = supervised design providing target, 2 = cocktail detection with multiplex primers, 3 = unsupervised design providing phage genome')
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

def update_primer3(tm, primer_size, kit):
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
        'PRIMER_NUM_RETURN' : 2,
    }
    return primer3_params

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
        print('==> Output directory already exists. Carrying on...\n')

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
                current_header = line[1:].strip().replace(" ", "_")
                current_seq = ''
            else:
                current_seq += line.strip()
        sequences[current_header] = current_seq
    return sequences     

def get_primers(header, seq, primer3_params):
    print("Designing primers for ", header)
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
    target_name = header[:50].split(" ")[0].strip()
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
    primers_df = pd.DataFrame(primers)
    return(primers_df)

def df_to_fasta(df, output_dir, output_fasta):
    records = []
    for index, row in df.iterrows():
        seq_record = SeqRecord(Seq(row['left_primer_sequence']),
                               id=row['left_primer_name'],
                               description="")
        records.append(seq_record)
        seq_record = SeqRecord(Seq(row['right_primer_sequence']),
                               id=row['right_primer_name'],
                               description="")
        records.append(seq_record)
    with open(output_dir + "/" + output_fasta + ".fna", "w") as handle:
        SeqIO.write(records, handle, 'fasta')

def run_mfe_index(db):
    subprocess.call('mfeprimer index -i {}'.format(db), shell=True)

def run_mfe(primers, db, output_dir, suffix, tm_spec):
    print("Running MFEprimer full analysis for " + primers + " ...")
    subprocess.call('mfeprimer -i {} -d {} -o {}/mfe_results_{}.txt -S 500 -t {}'.format(primers, db, output_dir, suffix, tm_spec), shell=True)
    
def run_mfe_dimers(primers, output_dir, suffix):
    print("Screening for dimers.")
    print("Running MFEprimer dimers...\n")
    subprocess.call('mfeprimer dimer -i {} --dg -10 -o {}/mfe_dimer_results_{}.txt'.format(primers, output_dir, suffix), shell=True)

def run_mfe_spec(primers, db, output_dir, suffix, tm_spec):
    print("Screening for specificity against Blast database: " + str(db))
    print("Running MFEprimer specificity...\n")
    subprocess.call('mfeprimer spec -i {} -d {} -o {}/mfe_spec_results_{}.txt -S 500 -t {}'.format(primers, db, output_dir, suffix, tm_spec), shell=True)

if __name__ == '__main__':
    
    args = parse_args()
    
    fasta_file = args.fasta_file
    output_dir = args.output_dir
    blast_db = args.blast_db
    sequences_dir = args.sequences_dir
    human_genome = args.human_genome
    runtype = args.runtype
    tm = args.tm_optimal
    primer_size = args.size_optimal
    kit = args.kit
    tm_spec = args.tm_specificity

    # Testing arguments
    testing(args)
    
    # Creare output directory  
    mkdir_outdir(output_dir)
    
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
            primer3_params = update_primer3(tm, primer_size, kit)
            primer_results = get_primers(header, seq, primer3_params)         
            all_primers_df = parse_primers(primer_results, header)
        
        # create prefiltered primers file
        all_primers_df.to_csv(output_dir + '/prefiltered_primers_file.csv', index=False, sep ='\t')
        print('\nPrefiltered primers stored in: ', output_dir + '/prefiltered_primers_file.csv')
        
        # create multifasta file with all predicted priners
        df_to_fasta(all_primers_df, output_dir, 'all_primers')
        print('Prefiltered primers in fasta format stored in: ', output_dir + '/all_primers.fna')
 
        print("\n--------------------------------------------------------")               
        print("Step 2 - test specificity of primers using MFEprimer against custom database")
        print("--------------------------------------------------------")
        if blast_db != None:
            ## by default, mfeprimer_index will check if the indexing has already been performed before proceeding
            run_mfe_index(blast_db)
            print("Screening for all features and specificity against Blast database: " + str(blast_db) + "\n")
            # run_mfe(output_dir + '/all_primers.fna',blast_db,output_dir,'blastdb', tm_spec)
            # testing running MFE on separate fasta files         
            grouped = all_primers_df.groupby('left_primer_name')
            for name, group in grouped:
                df_to_fasta(grouped.get_group(name), output_dir, name[:-1])
                run_mfe(output_dir + '/' + str(name[:-1]) + '.fna', blast_db, output_dir, str(name[:-1]), tm_spec)
                os.remove(output_dir + '/' + str(name[:-1]) + '.fna')
        else:
            print("No Blast database provided")
            
        print("\n--------------------------------------------------------")               
        print("Step 3 - test specificity of primers using MFEprimer against human genome")
        print("--------------------------------------------------------")    
        if human_genome == True:
            # run MFEprimer specificity on combined file against human genome
            run_mfe_spec(output_dir + '/all_primers.fna','/data/db/blastdb/hg38/GCA_000001405.29_GRCh38.p14_genomic.fna',output_dir,'hg38', tm_spec)
        else:
            print("No Human Genome screening required")
 
        if runtype == 2:
            print("\n--------------------------------------------------------")               
            print("Step 4 - test primer compatibility for multiplex and cocktail detection")
            print("--------------------------------------------------------")               
            # run MFEprimer dimers on combined file 
            run_mfe_dimers(output_dir + '/all_primers.fna', output_dir, 'all_primers')


           
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
        
