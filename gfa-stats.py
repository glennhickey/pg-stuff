#!/usr/bin/env python3

"""Super quick-and-dirty graph stats getter. 

number of nodes
number of components
number of edges
Length (bp)
Total path length (bp)
number of paths
Number of path steps

Usage  : zcat graph.gz.fz | gfa-stats.py

"""

import os,sys

node_to_length = {}
node_count = 0
step_count = 0
edge_count = 0
total_length = 0
total_path_length = 0
path_count = 0
step_count = 0

for line in sys.stdin:
    if line.startswith('S'):
        toks = line.split('\t')
        name, length = toks[1], len(toks[2].strip())
        node_to_length[name] = length
        node_count += 1
        total_length += length;
    elif line.startswith('L'):
        edge_count += 1
    elif line.startswith('P'):
        toks = line.split('\t')
        step_toks = toks[2].split(',')
        step_count += len(step_toks)
        for step in step_toks:
            total_path_length += node_to_length[step[:-1]]
        path_count += 1
    elif line.startswith('W'):
        toks = line.split('\t')
        stpe_count += toks[6].count('>') + toks[6].count('<')
        total_path_length += int(toks[5]) - int(toks[4])
        path_count += 1

print("number of nodes: {}".format(node_count))
print("number of edges: {}".format(edge_count))
print("Length (bp): {}".format(total_length))
print("Total path length (bp): {}".format(total_path_length))
print("number of paths: {}".format(path_count))
print("number of path steps: {}".format(step_count))

        
        
        
        
    
