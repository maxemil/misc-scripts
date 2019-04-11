#! /usr/bin/python3

import pandas as pd
from collections import defaultdict
import sys


def main(fnodes, species_map_file):
    species_map = {line.split()[0]: line.split()[1] for line in open(species_map_file)}

    species = set([species_map[line.split('\t')[1].split('&')[0]] for line in open(fnodes)])
    families = set([line.split('\t')[0] for line in open(fnodes)])
    members = pd.DataFrame(columns=species, index=families, dtype=int)
    members = members.fillna(0)

    for line in open(fnodes):
        family = line.split('\t')[0]
        spec = species_map[line.split('\t')[1].split('&')[0]]
        members[spec][family] += 1

    members.to_csv('families_count.tab', sep='\t', header=True, index_label='family')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
