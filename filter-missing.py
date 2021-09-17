#!/usr/bin/env python
"""
Filter out sites with > x% missing alleles
"""

import os, sys

if len(sys.argv) != 3:
    sys.stderr.write('usage: {} <vcf> <filter>\n'.format(sys.argv[0]))
    sys.exit(1)

vcf_path = sys.argv[1]
threshold = float(sys.argv[2])
assert threshold >= 0 and threshold <= 1
sys.stderr.write("[filter-missing.py] Input: {} Threshold: {}\n".format(vcf_path, threshold))

f_count = 0
t_count = 0
with open(vcf_path, 'r') as vcf_file:
    for line in vcf_file:
        filter = False
        if len(line) > 5 and not line.startswith('#'):
            gts = line.split('\t')[10:]
            total_count = 0
            missing_count = 0
            for gt in gts:
                for allele in gt.split('|'):
                    total_count += 1
                    if allele == '.':
                        missing_count += 1
            filter = float(missing_count) / float(total_count) > threshold
            t_count += 1
            if filter:
                f_count += 1
        if not filter:
            sys.stdout.write(line)                
sys.stderr.write("[filter-missing.py] Filtered {}/{} records\n".format(f_count, t_count))

            
        
