#!/usr/bin/env python3
"""
  17-02-16, Max Schön, <max-emil.schon@icm.uu.se>
  example usage:
   python3 parse_barrnap.py -l barrnap.log -r contigs.fasta -o rRNAs.fasta -f filter -t 1
   python3 parse_barrnap.py -l <(barrnap contigs.fasta) -r contigs.fasta -o rRNAs.fasta -f filter -t 1

  output:
   a fasta file containing the sequences of all predicted rRNA genes (or filtered)
"""
from Bio import SeqIO
import pandas as pd
import multiprocessing
import numpy as np
import time
import argparse
import textwrap

parser = argparse.ArgumentParser(
    description=textwrap.dedent(
        "parse a barrnap file and report rRNA sequences"
    ),
    epilog=textwrap.dedent(__doc__),
    formatter_class=argparse.RawDescriptionHelpFormatter
)

parser.add_argument("-l", "--logfile", required=True, type=str, help="barrnap log file or pipep-in output")
parser.add_argument("-r", "--reference", required=True, type=str, help="original sequences in fasta format")
parser.add_argument("-o", "--output", required=False, type=str, default='rRNA.fasta',
                    help="output file name, fasta format, default=rRNA.fasta")
parser.add_argument("-f", "--filter", required=False, default='rRNA',
                    choices=['16S', '5S', '23S', 'rRNA', '18S', '28S', 'SSU', 'LSU'],
                    help="filter term, default=rRNA (no filtering)", nargs='+')
parser.add_argument("-t", "--threads", required=False, type=int, default=1, help="number of threads, default=1")

args = parser.parse_args()

REF_DICT = {}

if 'SSU' in args.filter:
    args.filter += ['18S', '16S']
elif 'LSU' in args.filter:
    args.filter += ['28S', '23S']

def get_putative_rrna(contig, start, end, tag, local_list):
    if any([f  in  tag for f in args.filter]):
        local_list.append(SeqIO.SeqRecord(REF_DICT[contig].seq[start:end], contig + ":%s-%s" % (start, end),
                                          description=tag))


def parallel_apply(df):
    local_list = []
    df.apply(lambda x: get_putative_rrna(x[0], x[3], x[4], x[8], local_list), axis=1)
    return local_list


def load_references(ref_fasta):
    t0 = time.time()
    ref_dict = {r.id: r for r in SeqIO.parse(ref_fasta, 'fasta')}
    t1 = time.time()
    print("Loading sequences from the contigs took %s seconds" % (t1-t0))
    return ref_dict


def load_positions(barrnap_log):
    t0 = time.time()
    rRNAs = pd.read_csv(barrnap_log, header=None, skiprows=1, sep='\t')
    t1 = time.time()
    print("Loading rRNA position from the barrnap logfile took %s seconds" % (t1-t0))
    return rRNAs


def cut_sequences(rRNAs):
    t0 = time.time()
    # rRNAs.apply(lambda x: get_putativeRRNA(x[0], x[3], x[4], x[8]), axis=1)
    pool = multiprocessing.Pool(processes=args.threads)
    pool_output = pool.map(parallel_apply, np.array_split(rRNAs, args.threads))
    pool.close()
    pool.join()
    t1 = time.time()
    print("cutting rRNA sequences from contigs took %s seconds" % (t1-t0))
    return pool_output


def write_sequences(outfile, pool_output):
    t0 = time.time()
    if any([l for l in pool_output]):
        with open(outfile, 'w') as f:
            for local_list in pool_output:
                for rec in local_list:
                    SeqIO.write(rec, f, 'fasta')
    else:
        print("No output file written, no sequences extracted")
    t1 = time.time()
    print("writing rRNA sequences to file took %s seconds" % (t1-t0))


if __name__ == "__main__":
    REF_DICT = load_references(args.reference)
    rRNAs = load_positions(args.logfile)
    sequences = cut_sequences(rRNAs)
    write_sequences(args.output, sequences)
