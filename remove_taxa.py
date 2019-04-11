#!/usr/bin/env python3
"""
  18-09-14, Max Schön, <max-emil.schon@icm.uu.se>
  example usage:
   python3 parse_barrnap.py -l barrnap.log -r contigs.fasta -o rRNAs.fasta -f filter -t 1

"""
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", required=True, type=str, nargs='+', help="fasta file(s) to be searched")
parser.add_argument("-t", "--taxa", required=True, type=str, nargs='+', help="taxa to be removed")
args = parser.parse_args()


def remove_taxa_from_fasta(infile, taxa):
    selected_taxa = []
    for rec in SeqIO.parse(infile, 'fasta'):
        if not rec.id in taxa:
            selected_taxa.append(rec)
    with open(infile, 'w') as out:
        SeqIO.write(selected_taxa, out, 'fasta')

def main(args):
    for f in args.fasta:
        remove_taxa_from_fasta(f, args.taxa)


if __name__ == '__main__':
    main(args)
