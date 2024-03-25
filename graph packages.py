import networkx as nx
from ete3 import Tree
import matplotlib.pyplot as plt

def visualise_nx(tree): #visualises networkx graphs as hierarchical trees with root at the top
    pos = nx.nx_agraph.graphviz_layout(tree, prog="dot")
    nx.draw(tree, pos, with_labels=True)
    plt.savefig(str(tree)+'.png')
    plt.show()
    return


def nx_to_nw(root, graph): #converts a networkx graph to a newick string, a format for describing trees that is used by many packages designed for evolutionary/phylogenetic studies
    newick = ''
    a, b = graph.neighbors(root)
    if graph.out_degree(a) == 0:
        newick += '('+str(a)+','
    else:
        newick += '('+nx_to_nw(a, graph)+','
    if graph.out_degree(b) == 0:
        newick += str(b)+')'+str(root)
    else:
        newick += nx_to_nw(b, graph)+')'+str(root)
    return newick
def nw_to_ete3(root, graph): #uses the previously defined function to convert from networkx to ete3 graphs
    newickString = nx_to_nw(root, graph)
    newickString+=';'
    t = Tree(newickString, format = 8)
    return(t)
def ete3_to_nx(eteTree): #converts from ete3 to networkx graphs
    nxTree = nx.DiGraph()
    root = eteTree.get_tree_root()
    nxTree.add_node(root.name)
    for node in eteTree.traverse("levelorder"):
        if len(node.children) == 0:
            continue
        else:
            for j in node.children:
                nxTree.add_edge(node.name, j.name)
    return nxTree


def sanger_to_nx(eteTree): #converts from ete3 to networkx in a way that does not require all nodes to be labelled
    nxTree = nx.DiGraph()
    root = eteTree.get_tree_root()
    dummyName = '0'
    if root.name == '':
        nxTree.add_node(dummyName)
        dummyName = dummyName + '0'
    else:
        nxTree.add_node(root.name)
    for node in eteTree.traverse("levelorder"):
        if node.name == '':
            node.name = dummyName
            dummyName = dummyName + '0'
        if len(node.children) == 0:
            continue
        else:
            for j in node.children:
                if j.name =='' :
                    j.name = dummyName
                    dummyName = dummyName + '0'
                nxTree.add_edge(node.name, j.name)
    return nxTree