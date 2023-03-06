import ete3
import pandas as pd
import seaborn as sns
import numpy as np
import pylab as pl
from matplotlib import colors

cmap = sns.color_palette("magma_r", as_cmap=True)
norm = colors.Normalize(vmin=0, vmax=60)

def node_style_basic():
    """
    It generates a NodeStyle
    :return: ETE3 NodeStyle Object
    """
    style = ete3.NodeStyle()
    style["fgcolor"] = "#000000"
    style["size"] = 0
    # style["vt_line_color"] = "#000000"
    # style["hz_line_color"] = "#000000"
    style["vt_line_width"] = 2
    style["hz_line_width"] = 2
    style["vt_line_type"] = 0  # 0 solid, 1 dashed, 2 dotted
    style["hz_line_type"] = 0
    return style


def get_ancestor(nodes):
    ancestor = nodes[0]
    for n in nodes:
        ancestor = ancestor.get_common_ancestor(n)
    return ancestor


def plot_colormap():
    a = np.array([[0,60]])
    pl.figure(figsize=(9, 1.5))
    img = pl.imshow(a, cmap=cmap)
    pl.gca().set_visible(False)
    cax = pl.axes([0.1, 0.2, 0.8, 0.6])
    pl.colorbar(orientation='horizontal', cax=cax)
    pl.savefig("colorbar.pdf")

# def set_node_style(tree, style, leaves=False, condition=None):
#     for n in tree.traverse():
#         if leaves:
#             if n.is_leaf() is False:
#                 continue
#         if condition is None:
#             n.img_style = style()
#         else:
#             operator, node_feature, value = condition
#             if node_feature not in n.features:
#                 continue
#             if operator(getattr(n, node_feature), value):
#                 n.img_style = style()    

def simple_tree_style():
    ts = ete3.TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    # ts.scale = 10
    ts.mode = 'r'
    ts.allow_face_overlap = True
    ts.min_leaf_separation = 0
    ts.branch_vertical_margin = 0
    ts.legend.add_face(ete3.CircleFace(label='1000 copies',radius=scale(1000)/2, color="#65C18C"), column=0)
    ts.legend.add_face(ete3.CircleFace(label='5000 copies',radius=scale(5000)/2, color="#65C18C"), column=1)
    ts.legend.add_face(ete3.CircleFace(label='40 originations',radius=scale(40), color="#B667F1"), column=0)
    ts.legend.add_face(ete3.CircleFace(label='200 originations',radius=scale(200), color="#B667F1"), column=1)
    ts.legend.add_face( ete3.RectFace(label='100 duplications', width=scale(100)*4, height=scale(100)*2, fgcolor="#86C6F4", bgcolor="#86C6F4"), column=0)
    ts.legend.add_face( ete3.RectFace(label='500 duplications', width=scale(500)*4, height=scale(500)*2, fgcolor="#86C6F4", bgcolor="#86C6F4"), column=1)
    ts.legend.add_face( ete3.RectFace(label='100 transfers', width=scale(100)*4, height=scale(100)*2, fgcolor="#D3DEDC", bgcolor="#D3DEDC"), column=0)
    ts.legend.add_face( ete3.RectFace(label='500 transfers', width=scale(500)*4, height=scale(500)*2, fgcolor="#D3DEDC", bgcolor="#D3DEDC"), column=1)
    return ts


def layout(node):
    style = node.img_style
    # if node.is_leaf():
    N = ete3.AttrFace("name", fsize=40, fgcolor="black")
    ete3.faces.add_face_to_node(N, node, 1, position="branch-right")
    if all([feat in node.features for feat in ['copies', 'transfers', 'duplications', 'percL', 'originations']]):
    # if 'copies' in node.features:
        #C = CircleFace(radius=log(node.weight), color="Red", style="circle")
        copies = ete3.CircleFace(radius=node.copies/2, color="#65C18C")
        originations = ete3.CircleFace(radius=node.originations, color="#B667F1")
        duplications = ete3.RectFace(width=node.duplications*4, height=node.duplications*2, fgcolor="#86C6F4", bgcolor="#86C6F4")
        transfers = ete3.RectFace(width=node.transfers*4, height=node.transfers*2, fgcolor="#D3DEDC", bgcolor="#D3DEDC")
        # losses = ete3.RectFace(width=2, height=node.losses, fgcolor="Orange", bgcolor="Orange")

        style["hz_line_color"] = colors.to_hex(cmap(norm(node.percL)))
        style["vt_line_color"] = colors.to_hex(cmap(norm(node.percL)))
        style["size"] = 0
        style["vt_line_width"] = 10
        style["hz_line_width"] = 10

        node.add_face(copies, column=0, position="branch-right")
        node.add_face(originations, column=0, position="branch-right")
        # ete3.faces.add_face_to_node(losses, node, 2, position="float")
        ete3.faces.add_face_to_node(duplications, node, 0, position="branch-top")
        ete3.faces.add_face_to_node(transfers, node, 0, position="branch-bottom")
    # node.set_style(style)

def scale(x):
    # return int((x / 200)) * 2
    return np.sqrt((x / 200))  * 10

def add_genome_size(treeO, treeC, species_map, events, filename):
    events.to_csv(filename.replace('pdf', 'csv'), sep='\t')
    for n in treeO.traverse():
        name = ''
        if n.is_leaf():
            name = species_map[n.name]
            n.name = n.name.replace('_', ' ')
        else:
            name = get_ancestor(
                [treeC.get_leaves_by_name(species_map[l.name])[0] for l in n.get_leaves()]).name
            n.name = name
        n.add_features(copies=scale(events.loc[name]['copies']))
        n.add_features(transfers=scale(events.loc[name]['transfers']))
        n.add_features(duplications=scale(events.loc[name]['duplications']))
        n.add_features(losses=scale(events.loc[name]['losses']))
        n.add_features(percL=events.loc[name]['percL'])
        n.add_features(originations=scale(events.loc[name]['originations']))

    ts = simple_tree_style()
    # set_node_style(treeO, node_style_basic)
    for node in treeO.traverse():
        node.dist = 0.1
    treeO.ladderize(direction=1)
    treeO.render(filename, tree_style=ts)




df = pd.read_csv('output_new/events_summary.tsv', sep='\t')
df.index = df['node']    
df['percL'] = df['losses']*100/(df['copies']+df['losses']-df['transfers']-df['originations']-df['duplications'])

species_map = {line.split()[0]: line.split()[1] for line in open("input/map_species.txt")}
troot = ete3.PhyloTree('species_trees/NM-nDE.fasta.LGC60G4F.treeForALE_root.tree', format=0)
# tclean = ete3.PhyloTree('species_trees/NM-nDE.fasta.LGC60G4F.treeForALE_OrigNames.tree')
tclean = ete3.PhyloTree('species_trees/NM-nDE.fasta.LGC60G4F.treeForALE_clean.tree', format=1)

add_genome_size(troot, tclean, species_map, df, "Asgard_ALE_ete.pdf")
# plot_colormap()
