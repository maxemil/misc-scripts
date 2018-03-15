#! /usr/bin/env python3
from Bio import SeqIO
import argparse

parser.add_argument("-i", "--input", required=True, help="fasta file")
parser.add_argument("-o", "--output", required=True, help="output sorted fasta file")
args = parser.parse_args()

if __name__ == "__main__":
    seq_iterator = SeqIO.parse(args.input, 'fasta')
    seqs = [f for f in sorted(seq_iterator, key=lambda x: x.id)]
    SeqIO.write(seqs, open(args.output, 'w'), 'fasta')
