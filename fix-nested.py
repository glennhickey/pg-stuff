#!/usr/bin/env python3
"""
Correct LV and PS tags to account for missing sites in the VCF. 

Since snarl information in VCF is incomplete (why we are running this in the first place)
we use an interval tree to get the information
"""

import os, sys, gzip
from intervaltree import Interval, IntervalTree

if len(sys.argv) not in [2]:
    sys.stderr.write('fix-nested.py: correct LV and PS tags to account for missing sites in the VCF\n\n')
    sys.stderr.write('usage: {} <vcf> | bgzip > fixed.vcf.gz \n'.format(sys.argv[0]))
    sys.exit(1)

vcf_path = sys.argv[1]

if not vcf_path.endswith('vcf.gz'):
    sys.stderr.write('only .vcf.gz input supported\n')
    sys.exit(1)

def get_parent(toks):
    info = toks[7]
    itoks = info.split(';')
    parent = None
    for itok in itoks:
        if itok.startswith('PS='):
            return itok[3:]
    return None

def set_parent_lv(toks, parent, lv):
    info = toks[7]
    itoks = info.split(';')
    for i, itok in enumerate(itoks):
        if itok.startswith('PS='):
            if parent and lv > 0:
                itoks[i] = 'PS={}'.format(parent)
            else:
                itoks[i] = None
        if itok.startswith('LV='):
            itoks[i] = 'LV={}'.format(lv)
    toks[7] = ';'.join([itok for itok in itoks if itok])

# pass one: compute the nesting tree from reference coordinates
chrom_count = 0
interval_count = 0
interval_trees_by_chrom = {}
with gzip.open(vcf_path, 'r') as vcf_file:
    for line in vcf_file:
        line = line.decode('utf-8')
        if len(line) > 5 and not line.startswith('#'):
            toks = line.split('\t')
            chrom = toks[0]
            pos = int(toks[1])
            name = toks[2]
            ref_len = len(toks[3])
            if chrom not in interval_trees_by_chrom:
                interval_trees_by_chrom[chrom] = IntervalTree()
                chrom_count += 1
            interval_trees_by_chrom[chrom].addi(pos, pos + ref_len, name)
            interval_count += 1

sys.stderr.write("[fix-nested.py]: Computed {} interval trees containing a total of {} intervals\n".format(chrom_count, interval_count))

# pass two: use the interval tree to sort out LV and PS
with gzip.open(vcf_path, 'r') as vcf_file:
    for line in vcf_file:
        line = line.decode('utf-8')
        if len(line) > 5 and not line.startswith('#'):
            toks = line.split('\t')
            chrom = toks[0]
            pos = int(toks[1])
            name = toks[2]
            ref_len = len(toks[3])
            interval = Interval(pos, pos + ref_len, name)
            # find all overlapping intervals            
            overlaps = interval_trees_by_chrom[chrom].overlap(interval)
            # count containing intervals (and find smallest one)
            parent_interval = None
            num_parents = 0
            found_self = False
            for overlap in overlaps:
                if overlap != interval:
                    if overlap.begin <= interval.begin and overlap.end >= interval.end:
                        if overlap.begin < interval.begin or overlap.end > interval.end:
                            num_parents += 1
                            if parent_interval is None or (parent_interval.begin <= overlap.begin and parent_interval.end >= overlap.end):
                                parent_interval = overlap
                        else:
                            sys.stderr.write("[fix-nested.py]: Warning found duplicate intervals: {} and {}\n".format(interval, overlap))
                else:
                    found_self = True

            # sanity check
            assert(found_self)

            lv = num_parents
            ps = parent_interval.data if parent_interval else None
            assert (ps == None) == (lv == 0)
                    
            set_parent_lv(toks, ps, lv)
            sys.stdout.write('\t'.join(toks))
        else:
            sys.stdout.write(line)

            
        
