#!/usr/bin/env python3
"""
  17-02-16, Max Schön, <max-emil.schon@icm.uu.se>
  example usage:
"""
from Bio import SeqIO
import argparse
from itertools import groupby
import copy

parser = argparse.ArgumentParser()

parser.add_argument("-b", "--blast", required=True, type=str, help="blast result file")
parser.add_argument("-f", "--fasta", required=True, type=str, help="fasta input file")
parser.add_argument("-n", "--header", required=True, type=str, nargs='+', help="blast tab header")
parser.add_argument("-o", "--outfile", required=True, type=str, help="output with with longest alignment translated")
parser.add_argument("-t", "--translate",  action='store_true', help="output with with longest alignment translated")

args = parser.parse_args()

class Hsp:
    """
    Store information about single HSP in an alignment hit.
    """

    def __init__(self,entry, header):
        bt_fields = entry.split()
        self.members = {k:v for k,v in zip(header, bt_fields)}

    def __str__(self):
        try:
            f = "{}\t{}".format(self.members['qseqid'], self.members['sseqid'])
        except KeyError:
            f = "No idea what the name of query and subject are"
        return f


class BlastRecord:
    """
    Object representing a Blast Record.
    """

    def __init__(self, qseqid=None, hits=None):
        """Initialize Blast Record instance"""
        self.qseqid = qseqid
        self.hits = hits

    def __str__(self):
        """Return output string of BLAST record"""
        return self.qseqid

def parse_blast(handle, header):
    for qseqid, blasts in groupby(handle, lambda l: l.split()[0]):
        hits = []
        for line in blasts:
            hsp = Hsp(line.strip(), header)
            hits.append(hsp)
        yield BlastRecord(qseqid=qseqid, hits=hits)


def get_cds(seq, start, end, frame):
    cds_start = frame - 1
    while cds_start < start:
        cds_start += 3
    cds_end = frame -1
    while cds_end < end:
        cds_end += 3
        if cds_end + 3 > end:
            break
    return seq[cds_start:cds_end]


def get_seqs_fasta(fastafile):
    return {rec.id: rec for rec in SeqIO.parse(fastafile, 'fasta')}


def get_max_bitscore(blastfile, pos):
    max_bit = 0.0
    for line in open(blastfile):
        if float(line.split()[pos-1]) > max_bit:
            max_bit = float(line.split()[pos-1])
    return max_bit

def main(args):
    seqs = get_seqs_fasta(args.fasta)
    longest_cds = {}
    max_bit = get_max_bitscore(args.blast, 13)
    for rec in parse_blast(open(args.blast), args.header):
        for hsp in rec.hits:
            if float(hsp.members['bitscore']) > 0.9 * max_bit:
                sseqid = hsp.members['sseqid']
                sstart = int(hsp.members['sstart'])
                send = int(hsp.members['send'])
                sframe = int(hsp.members['sframe'])
                cds = get_cds(seqs[sseqid], sstart, send, sframe)
                try:
                    if longest_cds[sseqid][1] < len(cds):
                        longest_cds[sseqid] = (cds, len(cds))
                except:
                    longest_cds[sseqid] = (cds, len(cds))
    with open(args.outfile, 'w') as outfile:
        for k,v in longest_cds.items():
            seq_prot = copy.deepcopy(longest_cds[k][0])
            if args.translate:
                seq_prot.seq = longest_cds[k][0].seq.translate()
            else:
                seq_prot.seq = longest_cds[k][0].seq
            SeqIO.write(seq_prot, outfile, 'fasta')


if __name__ == '__main__':
    main(args)
