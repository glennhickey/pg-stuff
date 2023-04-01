#!/usr/bin/env python3

# parse the minigraph-split.log to get a bed file of filtered contigs

import os, sys
from collections import defaultdict

assert len(sys.argv) == 2

contig_name = None
prev_line_type = None
uf_fail = False
cur_table = {}
with open(sys.argv[1], 'r') as split_log:
    for line in split_log:
        if line.startswith('Query contig is ambiguous'):
            toks = line.replace('uf= infinity', 'uf=9999').strip().split()
            contig_name = toks[4]
            contig_length = int(toks[5][4:])
            query_coverage = float(toks[6][4:])
            coverage_threshold = float(toks[8][:-1])
            uf = float(toks[9][3:])
            uf_threshold = float(toks[11][:-1])
            pipe_pos = contig_name.find('|')
            dot_pos = contig_name.find('.')
            if dot_pos > 0 and pipe_pos > 0 and dot_pos < pipe_pos:
                contig_name = contig_name.replace('.', '#', 1)
            contig_name = contig_name.replace('id=', '').replace('|', '#')
            print('{}\t{}\t{}\t{}'.format(contig_name, 0, contig_length, '{};{};{};{}'.format(query_coverage, coverage_threshold, uf, uf_threshold)))
