import argparse
from numpy.random import random
from networkx import nx

#setting default parameter values
p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument('-v', '--verbose', action='store_true', default = False)
p.add_argument('--pGS', type = float, default = 0.01, help = 'PGC specification probability')
p.add_argument('--pdeath', type = float, default = 0.1, help = 'cell death probability')
p.add_argument('--GScycles', default = '2, 4', help = 'PGC specification cycles: start,end')
p.add_argument('--tc', type = float, default = 0.5, help = 'cycle time (d)')
p.add_argument('--start_day', type = float, default = 0.0, help = 'initial time (d)')
p.add_argument('--start_reps', type = int, default = 0, help = 'initial reps (from zygote)')
p.add_argument('--max_reps_follow', type = int, default = 6, help = 'max reps to follow')
p.add_argument('-n', '--nclones', type = int, default = 500, help = 'number of clones to simulate')
args = p.parse_args(args=[])


(GS_start, GS_end) = [int(x) for x in args.GScycles.split(',')]

#model of germline specification
def cbranch(clone, ctype, time, reps, branch, graph):
    if random() < args.pdeath:  #simulation of clone ends if cell dies
        return (graph)
    time += args.tc
    reps += 1
    if ctype == 'S' and reps >= GS_start and reps <= GS_end: # test for germline specification
        if random() < args.pGS:
            ctype = 'G'
    if ctype == 'S' and reps >= GS_end: # don't follow somatic lineages after germline specification
        return (graph)
    if reps >= args.max_reps_follow:
        return (graph)
    graph.add_edge(branch, branch+'1')
    graph.add_edge(branch, branch+'0')
    graph = cbranch(clone, ctype, time, reps, branch + '1', graph)
    graph = cbranch(clone, ctype, time, reps, branch + '0', graph)

    return(graph)

#simpler cell proliferation model
def simpleModel(ctype, time, reps, branch, graph):
    if random() < args.pdeath:  #simulation of clone ends if cell dies
        return (graph)
    time += args.tc
    reps += 1
    if reps >= args.max_reps_follow:
        return (graph)
    graph.add_edge(branch, branch+'1')
    graph.add_edge(branch, branch+'0')
    graph = simpleModel(ctype, time, reps, branch + '1', graph)
    graph = simpleModel(ctype, time, reps, branch + '0', graph)
    return(graph)
