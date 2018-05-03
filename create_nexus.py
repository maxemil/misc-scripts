from collections import defaultdict
import ete3
from ETE3_Utils import *
import argparse


def defaultdict_defaultdict():
    return defaultdict(defaultdict_string)


def get_clade_dict(tax_map):
    ncbi = ete3.ncbi_taxonomy.NCBITaxa()
    taxon_clade = defaultdict(defaultdict_defaultdict)
    for k, v in tax_map.items():
        lineage = ncbi.get_lineage(v)
        for l in lineage:
            rank = ncbi.get_rank([l])[l]
            if rank != 'no rank':
                taxon_clade[k][rank] = ncbi.get_taxid_translator([l])[l]
    return taxon_clade


def collapse_groups_at_rank(tree, taxon_clade, annotation_dict, rank):
    groups = set([d[rank] for d in taxon_clade.values()])
    max_dist = max([tree.get_distance(l) for l in tree.get_leaves()])
    for g in groups:
        species_set = [k for k,v in taxon_clade.items() if v[rank] == g]
        ancestor = tree.get_common_ancestor(species_set)
        if set(species_set) == set([l for l in ancestor.get_leaves()]):
            dist = max([tree.get_distance(l) for l in species_set])
            ancestor.name = 'internal_%s' % g.replace(' ', '_')
            annotation_dict[ancestor.name] = "[&label=%i,!collapse={\"collapsed\",%f},!name=\"%s\"]" % (ancestor.support, max_dist-dist, g)
        else:
            # print('found a nonmonophyletic group: %s' % g)
            pass

def collapse_groups(tree, taxon_clade):
    annotation_dict = {}
    for rank in ['genus', 'family', 'order', 'class', 'phylum', 'superkingdom']:
        collapse_groups_at_rank(tree, taxon_clade, annotation_dict, rank)

    tree_string = tree.write(format=1)
    for k,v in annotation_dict.items():
        tree_string = tree_string.replace(k,v)
    return tree_string


def get_taxa_block_nexus(tree):
    taxa_block = "begin taxa;\n"
    taxa_block += "\tdimensions ntax=%s;\n" % len(tree.get_leaves())
    taxa_block += "\ttaxlabels\n"
    for t in tree.get_leaves():
        if t.name.startswith("'lcl"):
            taxa_block += "\t%s[&!color=#ff0000]\n" % t.name
        else:
            taxa_block += "\t%s\n" % t.name
    taxa_block += ";\n"
    taxa_block += "end;\n"
    return taxa_block

def get_figtree_block_nexus():
    return """
begin figtree;
	set nodeLabels.displayAttribute="label";
	set nodeLabels.isShown=true;
	set nodeLabels.fontSize=8;
	set tipLabels.fontSize=10;
	set trees.orderType="increasing";
end;
"""


def write_nexus(outfile, tree_string, tree):
    with open(outfile, 'w') as out:
        print('#NEXUS', file=out)
        print(get_taxa_block_nexus(tree), file=out)
        print('begin trees;', file=out)
        print('\ttree tree_1 = [&R]', file=out)
        print("\t%s" % tree_string, file=out)
        print('end;', file=out)
        print(get_figtree_block_nexus(), file=out)


def prepare_node_leave_names(tree):
    for l in tree.get_leaves():
        l.name = "'%s'" % l.name
    for n in tree.traverse():
        if not n.name:
            n.name = 'INTERNAL_%i_SUPPORT' % n.support


def get_tax_ids(tree):
    tax_map = defaultdict(defaultdict_string)
    for l in tree.get_leaves():
        try:
            s = l.name.split('.')
            tax_map[l] = int(s[1])
        except:
            print(l.name)
    return tax_map

def main(treefile, ingroup, outfile):
    tree = ete3.PhyloTree(treefile, format=2)

    tax_map = get_tax_ids(tree)
    taxon_clade = get_clade_dict(tax_map)


    root_tree(tree, ingroup)
    tree.ladderize(direction=1)

    prepare_node_leave_names(tree)

    tree_string = collapse_groups(tree, taxon_clade)
    tree_string = tree_string.replace('INTERNAL_', "[&label=").replace('_SUPPORT',']')
    write_nexus(outfile, tree_string, tree)

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True,
                    help="input tree file, newick format")
parser.add_argument("-o", "--output", required=False, default="tree.nex",
                    help="output file in nexus format, prepared for FigTree")

args = parser.parse_args()

if __name__ == '__main__':
    ingroup = ['Magnetococcus_marinus_MC_1.156889',
               'Brucella_abortus_S19.430066',
               'Zymomonas_mobilis_subsp__mobilis_ZM4___ATCC_31821.264203',
               'Rhodospirillum_rubrum_ATCC_11170.269796',
               'Holospora_obtusa.49893',
               'Rickettsia_prowazekii_str__Madrid_E.272947',
               'Neorickettsia_sennetsu_str__Miyayama.222891']
    main(args.input ,ingroup, args.output)
