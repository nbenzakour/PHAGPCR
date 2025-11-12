#!/usr/bin/python3

#title           :bz_gene_extractor.py
#description     :This script will extract a gene fasta entry from a genbank file, based on a keyword
#author          :Nouri L. Ben Zakour
#date            :20240216
#version         :0.1
#=======================================================================================================

import argparse
import os
from Bio import SeqIO
import os

def extract_gene_sequence(input_file, word, outdir, skip):
    # Create the output directory if it does not exist
    os.makedirs(outdir, exist_ok=True)
    records = SeqIO.parse(input_file, "genbank")
    basename = os.path.splitext(os.path.basename(input_file))[0]
    word_found = False  # Initialize word_found as False for each file
    for record in records:
        for feature in record.features:
            if "product" in feature.qualifiers and word in feature.qualifiers["product"][0]:
                sequence = feature.location.extract(record).seq
                locus_tag = feature.qualifiers["locus_tag"][0]
                if len(sequence) < skip:
                    print(f"{locus_tag} \t-> matches '{word}' \tlength: {len(sequence)} bp but below cutoff {skip} bp")
                    continue
                output_file = os.path.join(outdir, f"{locus_tag}_{word.replace(' ', '_')}_{len(sequence)}bp.fna")
                with open(output_file, "w") as f:
                    f.write(f">{locus_tag}_{word.replace(' ', '_')}_{len(sequence)}bp\n{sequence}\n")
                print(f"{locus_tag} \t-> matches '{word}' \tlength: {len(sequence)} bp \t saved to {output_file}")
                word_found = True  # Set word_found to True when the word is found
    return word_found

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-i', '--input', help='Input file')
group.add_argument('-d', '--directory', help='Input directory')
parser.add_argument('-w', '--word', required=True, help='Word to search for in the product qualifiers')
parser.add_argument('-w2', '--word2', help='Second word to search for if the first word is not found')
parser.add_argument('-o', '--outdir', required=True, help='Output directory')
parser.add_argument("-sk", "--skip", default=0, type=int, help="Minimum sequence length to consider")

args = parser.parse_args()

# Then in your function call, you can check which one is not None and use it
if args.input is not None:
    word_found = extract_gene_sequence(args.input, args.word, args.outdir, args.skip)
    if not word_found and args.word2:
        word_found = extract_gene_sequence(args.input, args.word2, args.outdir, args.skip)
    if not word_found:
        print(f"{args.input} \t-> no CDS product matching {args.word} or {args.word2}. This may be an annotation issue or a true absence of that function.")
elif args.directory is not None:
    for file in os.listdir(args.directory):
        if file.endswith((".gb",".gbk")):
            word_found = extract_gene_sequence(os.path.join(args.directory, file), args.word, args.outdir, args.skip)
            if not word_found and args.word2:
                word_found = extract_gene_sequence(os.path.join(args.directory, file), args.word2, args.outdir, args.skip)
            if not word_found:
                print(f"{file} \t-> no CDS product matching {args.word} or {args.word2}. This may be an annotation issue or a true absence of that function.")
        print('\n')