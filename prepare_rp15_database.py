from Bio import SeqIO
import glob
import re

def replace(string, patterns):
    for p in patterns:
        string = string.replace(p[0], p[1])
    return string


patternfile = '/local/one/people/MaxEmil/Environmental_Genomes_RP15/Alphaproteobacteria_BM001/ref_names_updated.txt'
patterns = {}
for line in open(patternfile):
    line = line.strip().split('\t')
    patterns[line[0]] = line[1]

infiles = glob.glob('*.faa')
for f in infiles:
    with open(f.replace('.faa', '.clean.faa'), 'w') as outfile:
        for rec in SeqIO.parse(f, 'fasta'):
            rec.id = re.sub("_gi\|.*", "", rec.id)
            rec.id = re.sub("=", "_", rec.id)
            rec.id = re.sub("[50S|30S]_ribosomal_protein_.+", "", rec.id)
            rec.id = re.sub("_intID[0-9]+", "", rec.id)
            rec.id = re.sub("_[0-9]+_[0-9]+", "", rec.id)
            rec.id = patterns[rec.id]
            SeqIO.write(rec, outfile, 'fasta')
