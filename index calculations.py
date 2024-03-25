import numpy as np
import networkx as nx
from networkx.algorithms.traversal.depth_first_search import dfs_tree

#Colless index
def colless(tree):
    colless = 0
    for i in tree.nodes:
        if tree.out_degree(i) == 0: #skip leaves
            continue
        else:
            nl = []
            for j in tree.neighbors(i):
                subtree = dfs_tree(tree, j) # create subtree rooted at j
                num=0
                for k in subtree.nodes:
                    if subtree.out_degree(k) == 0:
                        num+=1
                nl.append(num)
            diff = abs(nl[0] - nl[1])
        colless+= diff
    return(colless)

#J1 (Lemant et al, 2022)
def S_i(tree, i):
    subtree = dfs_tree(tree, i)
    fv = 0
    for v in subtree.nodes:
        fv+=1
    return(fv)
def S_star_i (tree, i):
    Si = S_i(tree, i)
    return (Si-1)

def W_i(tree, i):
    Si = S_star_i(tree, i)
    Wij = 0
    if tree.out_degree(i)==0:
        return(0)
    for j in tree.neighbors(i):
        pij = S_i (tree, j)/Si
        Wij+= -1*pij*np.log2(pij)
    return(Wij)
def J1(tree):
    num = 0
    denom = 0
    for k in tree.nodes:
        if tree.out_degree(k)==0: #skip le
            continue
        else:
            Sk = S_star_i(tree, k)
            Wi = W_i(tree, k)
            num+= Sk*Wi
            denom += Sk
    if denom == 0:
        return(np.nan)
    else:
        return (num/denom)