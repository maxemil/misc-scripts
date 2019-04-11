#! /usr/bin/env python3
import re
from Bio import SeqIO, Seq
from ete3 import ncbi_taxonomy
ncbi = ncbi_taxonomy.NCBITaxa()


def main(db_file, query_file):

    seqdict = {rec.id: rec for rec in SeqIO.parse(db_file, 'fasta')}
    pattern = re.compile('OX=[0-9]*')
    seqs = [line.strip() for line in open(query_file)]

    for seq in seqs:
        rec = seqdict[seq]
        seqid = rec.id.split('|')[1]
        seqtax = pattern.search(rec.description).group(0).split('=')[1]
        rec.id = '{}.{}'.format(seqtax, seqid)
        rec.description = ''
        print(rec.id)

def is_taxid_in_clade(taxid, clade):
    lineage = ncbi.get_lineage(taxid)
    for l in ncbi.get_taxid_translator(lineage).values():
        if l == clade:
            return True
    return False

if __name__ == '__main__':
    main("/local/two/databases/uniprot_all.fasta", "silix_unlab_60_90_169.ids")
