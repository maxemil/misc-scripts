#! /usr/bin/python3

import argparse
import sys
from Bio import SearchIO, SeqIO
from collections import defaultdict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True,
                    help="input blast xml")
parser.add_argument("-o", "--output", required=True,
                    help="output fasta")
parser.add_argument("-f", "--fasta", required=False,
                    help="input fasta with original sequences (database)")
parser.add_argument("-m", "--blast_method", required=True,
                    help="choose from [tblastn, blastp]")
parser.add_argument("-t", "--threshold", required=False, default=0, type=int,
                    help="minumm length for sequences to be written")
args = parser.parse_args()


def get_unique_blastp_hits(infile, fasta):
    hits = set()
    for aln in SearchIO.parse(infile, 'blast-xml'):
        for hsp in aln.hsps:
            hits.add(hsp.hit_id)
    seqs = {rec.id:rec for rec in SeqIO.parse(fasta, 'fasta')}
    return [seqs[hit] for hit in hits]


def get_unique_tblastn_hits(infile):
    hits = defaultdict(lambda: (float("inf"),0,0, ""))
    for aln in SearchIO.parse(infile, 'blast-xml'):
        for hsp in aln.hsps:
            if hsp.hit_start < hits[hsp.hit_id][0] and hsp.hit_end > hits[hsp.hit_id][1]:
                hits[hsp.hit_id] = (hsp.hit_start, hsp.hit_end, hsp.aln_span, hsp.hit)
            elif hsp.aln_span > hits[hsp.hit_id][2]:
                hits[hsp.hit_id] = (hsp.hit_start, hsp.hit_end, hsp.aln_span, hsp.hit)
                print("hsp was longer that the other")
            elif hsp.hit_start < hits[hsp.hit_id][0]:
                print("{} < {} but {} not > {}".format(hsp.hit_start, hits[hsp.hit_id][0], hsp.hit_end, hits[hsp.hit_id][1]))
            elif hsp.hit_end > hits[hsp.hit_id][1]:
                print("{} not < {} but {} > {}".format(hsp.hit_start, hits[hsp.hit_id][0], hsp.hit_end, hits[hsp.hit_id][1]))
    return [hit[3] for hit in hits.values()]


def main(args):
    if args.blast_method == 'tblastn':
        hits = get_unique_tblastn_hits(args.input)
    elif args.blast_method == 'blastp' and args.fasta:
        hits = get_unique_blastp_hits(args.input, args.fasta)
    else:
        print('error, check input files')
        sys.exit()
    with open(args.output, 'w') as outhandle:
        for hit in hits:
            hit.seq = hit.seq.ungap('-').ungap('*')
            hit.description = ''
            if len(hit.seq) > args.threshold:
                SeqIO.write(hit, outhandle, 'fasta')

if __name__ == '__main__':
    main(args)
