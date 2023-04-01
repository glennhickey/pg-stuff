#!/usr/bin/env python3

# parse the minigraph-split.log to get stats for the top-two contigs by converage for each ambiguous mapping

import os, sys
from collections import defaultdict

assert len(sys.argv) == 2

single_stats_count = defaultdict(int)
single_stats_bp = defaultdict(int)
pair_stats_count = defaultdict(int)
pair_stats_bp = defaultdict(int)
single_stats_uf_fail_count = defaultdict(int)
single_stats_uf_fail_bp = defaultdict(int)
pair_stats_uf_fail_count = defaultdict(int)
pair_stats_uf_fail_bp = defaultdict(int)

contig_name = None
prev_line_type = None
uf_fail = False
cur_table = {}
with open(sys.argv[1], 'r') as split_log:
    for line in split_log:
        toks = line.replace('uf= infinity', 'uf=9999').strip().split()
        # classify the line
        line_type = None
        if line.startswith('Assigned contig'):
            line_type = 'A'
            contig_name = toks[4]
        elif line.startswith('Query contig is ambiguous'):
            line_type = 'Q'
            contig_name = toks[4]
            query_coverage = float(toks[6][4:])
            coverage_threshold = float(toks[8][:-1])
            uf = float(toks[9][3:])
            uf_threshold = float(toks[11][:-1])
            assert query_coverage < coverage_threshold or uf < uf_threshold
            uf_fail = uf < uf_threshold
        elif line.startswith(' Reference contig mappings:'):
            line_type = 'R'
            assert contig_name is not None
        elif len(toks) == 2 and prev_line_type in ['R', 'C']:
            line_type = 'C'
            ref_contig_name = toks[0][:-1] # strip the :
            coverage = int(toks[1])

        # update the coverage table
        if line_type == 'C':
            # add a contig
            cur_table[ref_contig_name] = coverage
        else:            
            if len(cur_table) > 1:
                # add the table to our stats
                vals = sorted(cur_table.values(), reverse=True)[:2]
                # top 2 by coverage (python fail)
                items = []
                for contig, cov in cur_table.items():
                    if cov in vals and len(items) < 2:
                        items.append((contig, cov))
                assert len(items) == 2
                if items[0][0] > items[1][0]:
                    items = (items[1], items[0])

                # single coverage
                for contig, cov in items:
                    single_stats_count[contig] += 1                        
                    single_stats_bp[contig] += cov
                    if uf_fail:
                        single_stats_uf_fail_count[contig] += 1                        
                        single_stats_uf_fail_bp[contig] += cov

                #pair coverage
                contig_pair = items[0][0], items[1][0]
                pair_stats_count[contig_pair] += 1
                pair_stats_bp[contig_pair] += items[0][1] + items[1][1]
                if uf_fail:
                    pair_stats_uf_fail_count[contig_pair] += 1
                    pair_stats_uf_fail_bp[contig_pair] += items[0][1] + items[1][1]
            # new cur table
            cur_table = {}
            
        prev_line_type = line_type

# header
print("ref-contig1\t\tref-contig-2\tbp\tcontigs\tbp-uf\tcontigs-uf")

for ref_contig in sorted(single_stats_bp.keys()):
    print("{}\t{}\t{}\t{}\t{}\t{}".format(ref_contig, '.',
                                          single_stats_bp[ref_contig],
                                          single_stats_count[ref_contig],
                                          single_stats_uf_fail_bp[ref_contig],
                                          single_stats_uf_fail_count[ref_contig]))

for ref_contig1, ref_contig2 in sorted(pair_stats_bp.keys()):
    print("{}\t{}\t{}\t{}\t{}\t{}".format(ref_contig1, ref_contig2,
                                          pair_stats_bp[(ref_contig1, ref_contig2)],
                                          pair_stats_count[(ref_contig1, ref_contig2)],
                                          pair_stats_uf_fail_bp[(ref_contig1, ref_contig2)],
                                          pair_stats_uf_fail_count[(ref_contig1, ref_contig2)]))
    
        
