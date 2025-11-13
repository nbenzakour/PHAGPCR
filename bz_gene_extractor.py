#!/usr/bin/env python3
"""
Extract gene sequences from GenBank files based on product qualifier keywords.

This script searches GenBank annotation files for genes matching specific
product keywords and extracts their sequences to FASTA files. Useful for
bulk extraction of specific genes (e.g., primase, polymerase) from phage
genome annotations.

Usage:
    python bz_gene_extractor.py -i phage.gb -w "primase" -o output/
    python bz_gene_extractor.py -d genbank_dir/ -w "DNA polymerase" \
-w2 "polymerase" -o output/ -sk 200

Author: Nouri L. Ben Zakour
Version: 0.2
Date: 2024-02-16
"""
import argparse
import os

from Bio import SeqIO

def extract_gene_sequence(input_file, word, outdir, skip):
    """
    Extract gene sequences from GenBank file based on product qualifier.

    Args:
        input_file: Path to GenBank file
        word: Keyword to search in product qualifiers
        outdir: Output directory for extracted sequences
        skip: Minimum sequence length threshold (bp)

    Returns:
        bool: True if keyword found, False otherwise
    """
    # Validate input file
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")

    # Create the output directory if it does not exist
    os.makedirs(outdir, exist_ok=True)

    try:
        records = SeqIO.parse(input_file, "genbank")
        basename = os.path.splitext(os.path.basename(input_file))[0]
        word_found = False

        for record in records:
            for feature in record.features:
                if "product" in feature.qualifiers and \
                   word in feature.qualifiers["product"][0]:
                    sequence = feature.location.extract(record).seq
                    locus_tag = feature.qualifiers["locus_tag"][0]

                    if len(sequence) < skip:
                        print(f"{locus_tag} \t-> matches '{word}' \t\
length: {len(sequence)} bp but below cutoff {skip} bp")
                        continue

                    output_file = os.path.join(
                        outdir,
                        f"{locus_tag}_{word.replace(' ', '_')}_\
{len(sequence)}bp.fna")
                    with open(output_file, "w") as f:
                        f.write(f">{locus_tag}_{word.replace(' ', '_')}_\
{len(sequence)}bp\n{sequence}\n")
                    print(f"{locus_tag} \t-> matches '{word}' \tlength: \
{len(sequence)} bp \t saved to {output_file}")
                    word_found = True

        return word_found
    except Exception as e:
        print(f"Error processing {input_file}: {str(e)}")
        return False

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-i', '--input', help='Input GenBank file')
group.add_argument('-d', '--directory',
                   help='Input directory containing GenBank files')
parser.add_argument('-w', '--word', required=True,
                    help='Word to search for in product qualifiers')
parser.add_argument('-w2', '--word2',
                    help='Second word to search if first not found')
parser.add_argument('-o', '--outdir', required=True,
                    help='Output directory for extracted sequences')
parser.add_argument("-sk", "--skip", default=0, type=int,
                    help="Minimum sequence length to consider (bp)")

args = parser.parse_args()

# Process input file or directory
if args.input is not None:
    try:
        word_found = extract_gene_sequence(args.input, args.word,
                                           args.outdir, args.skip)
        if not word_found and args.word2:
            word_found = extract_gene_sequence(args.input, args.word2,
                                               args.outdir, args.skip)
        if not word_found:
            print(f"{args.input} \t-> no CDS product matching \
{args.word} or {args.word2}. Annotation issue or true absence.")
    except FileNotFoundError as e:
        print(f"Error: {e}")
elif args.directory is not None:
    if not os.path.isdir(args.directory):
        print(f"Error: Directory not found: {args.directory}")
    else:
        for file in os.listdir(args.directory):
            if file.endswith((".gb", ".gbk")):
                try:
                    file_path = os.path.join(args.directory, file)
                    word_found = extract_gene_sequence(file_path,
                                                       args.word,
                                                       args.outdir,
                                                       args.skip)
                    if not word_found and args.word2:
                        word_found = extract_gene_sequence(file_path,
                                                           args.word2,
                                                           args.outdir,
                                                           args.skip)
                    if not word_found:
                        print(f"{file} \t-> no CDS product matching \
{args.word} or {args.word2}. Annotation issue or true absence.")
                except Exception as e:
                    print(f"Error processing {file}: {str(e)}")