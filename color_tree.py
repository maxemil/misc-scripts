import argparse
import ete3

def parse_colors(tree, colors):
    tax2color = {}
    for line in open(colors):
        line = line.strip().split()
        tax2color[line[0]] = line[1]
    return tax2color
    
def get_color_taxon(name, tax2color):
    if name in tax2color:
        return tax2color[name]
    elif any([t in name for t in tax2color.keys()]):
        for t in tax2color.keys():
            if t in name:
                return tax2color[t]
    else:
        return "#000000"

def get_taxa_block_nexus(tree, tax2color):
    taxa_block = "begin taxa;\n"
    taxa_block += "\tdimensions ntax={};\n".format(len(tree.get_leaves()))
    taxa_block += "\ttaxlabels\n"
    for t in tree.get_leaves():
        taxa_block += "\t{}[&!color={}]\n".format(t.name, get_color_taxon(t.name, tax2color))
    taxa_block += ";\n"
    taxa_block += "end;\n"
    return taxa_block

def get_figtree_block_nexus():
    return """
begin figtree;
	set branchLabels.displayAttribute="label";
	set branchLabels.isShown=true;
	set branchLabels.fontSize=10;
	set tipLabels.fontSize=12;
	set trees.orderType="increasing";
end;
"""

def write_nexus(outfile, tree, tax2color):
    with open(outfile, 'w') as out:
        print('#NEXUS', file=out)
        print(get_taxa_block_nexus(tree, tax2color), file=out)
        print('begin trees;', file=out)
        print('\ttree tree_1 = [&R]', file=out)
        print("\t{}".format(tree.write(format=2)), file=out)
        print('end;', file=out)
        print(get_figtree_block_nexus(), file=out)

def main(treefile, outfile, colors):
    tree = ete3.PhyloTree(treefile, format=2)
    tree.ladderize(direction=1)
    tax2color = parse_colors(tree, colors)
    write_nexus(outfile, tree, tax2color)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True,
                    help="input tree file, newick format")
parser.add_argument("-c", "--colors", required=True,
                    help="input file with color codes for taxa, as a tab separated list")
parser.add_argument("-o", "--output", required=False, default="tree.nex",
                    help="output file in nexus format, prepared for FigTree")
args = parser.parse_args()

if __name__ == '__main__':
    main(args.input, args.output, args.colors)