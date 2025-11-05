#! /usr/bin/env python3

from Bio import AlignIO, SeqIO, Seq, SeqRecord
from collections import OrderedDict
import Bio
import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, nargs='+',
    help='fasta files to be concatenated')
parser.add_argument('-out', type=str, default='concat.fasta',
        help='if no output file is given, print alignment to concat.fasta')
parser.add_argument('-s', action='store_true',
        help='get basic stats for alignments')
parser.add_argument('-t', '--threshold', type=int, default=1,
        help='minimum number of genes that must be present in an organism')
parser.add_argument('-ps', '--partition_size', type=int, default=0,
        help='size of the partitions')
parser.add_argument('-p', '--partitions', action='store_true',
        help='print a partitions file in RAxML style')
parser.add_argument('-sep', '--seperator', default='_',
        help='seperator between species and gene id')
args = parser.parse_args()

def main(args):
    df = pd.DataFrame(columns = args.fasta)

    alignments = {}
    lengths = {}

    for f in args.fasta:
        try:
            aln = AlignIO.read(f, 'fasta')
            if aln.get_alignment_length() == 0:
                print('alignment file %s contains no sequences' % f)
                continue
            alignments[f] = aln
            lengths[f] = '-' * aln.get_alignment_length()
        except:
            print('alignment file %s not readable' % f)

        for rec in SeqIO.parse(f, 'fasta'):
            if rec.id.split(args.seperator)[0] in df.index:
                df.loc[rec.id.split(args.seperator)[0], f] = str(rec.seq)
            else:
                df.loc[rec.id.split(args.seperator)[0]] = pd.Series(dtype=str)
                df.loc[rec.id.split(args.seperator)[0], f] = str(rec.seq)
    if args.s:
        for col in df:
            try:
                print(col)
                print('# of unique sequences:', df[col].dropna().size)
                print('# of positions:', len(df[col].dropna()[0]))
            except:
                print('no sequences found in %s' % col)
    df = df.loc[df.count(axis=1) >= args.threshold]
    print("{} ({:.2f}%) missing genes".format(df.size - sum(df.count(0)),
                        ((df.size - sum(df.count(0))) / df.size) * 100),
                        file=sys.stderr)
    df.fillna(value=lengths, axis='rows', inplace=True)

    with open(args.out, 'w') as out:
        for row in df.itertuples():
            rec = SeqRecord.SeqRecord(id=row[0],
                                      description='',
                                      seq=Seq.Seq(''.join(row[1:])))
            SeqIO.write(rec, out, 'fasta')

    if args.partitions:
        partitions = OrderedDict()
        start = 1
        for i, gene in enumerate(df.iteritems()):
            if args.partition_size > 0:
                p_count = 1
                gene_end = start + len(gene[1][0]) - 1
                while start + args.partition_size < gene_end:
                    end = start + args.partition_size
                    partitions[gene[0] + str(p_count)] = (start, end)
                    start = end + 1
                    p_count += 1
                if gene_end - start < 30:
                    start, end = partitions[gene[0] + str(p_count - 1)]
                    partitions[gene[0] + str(p_count - 1)] = (start, gene_end)
                else:
                    partitions[gene[0] + str(p_count)] = (start, gene_end)
                start = gene_end + 1
            else:
                end = start + len(gene[1][0]) - 1
                partitions[gene[0] + str(p_count)] = (start, end)
                start = end + 1
        with open("{}.partitions".format(args.out), 'w') as pout:
            print("#nexus", file=pout)
            print("begin sets;", file=pout)
            for k, v in partitions.items():
                print("charset {} = {}-{};".format(k, v[0], v[1]), file=pout)
            print("end;", file=pout)

if __name__ == '__main__':
    main(args)
