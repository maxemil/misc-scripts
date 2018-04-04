#! /usr/bin/env python3

from Bio import AlignIO, SeqIO, Seq, SeqRecord
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
args = parser.parse_args()

def main(args):
    df = pd.DataFrame(columns = args.fasta)

    alignments = {}
    lengths = {}

    for f in args.fasta:
        try:
            aln = AlignIO.read(f, 'fasta')
            alignments[f] = aln
            lengths[f] = '-' * aln.get_alignment_length()
        except:
            print('alignment file %s not readable' % f)

        for rec in SeqIO.parse(f, 'fasta'):
            if rec.id in df.index:
                df[f].loc[rec.id] = str(rec.seq)
            else:
                df.loc[rec.id] = pd.Series()
                df[f].loc[rec.id] = str(rec.seq)
    if args.s:
        for col in df:
            print(col)
            print('# of unique sequences:', df[col].dropna().size)
            print('# of positions:', len(df[col].dropna()[0]))

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


if __name__ == '__main__':
    main(args)
