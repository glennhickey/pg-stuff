#!/usr/bin/env python
"""
Filter out nested records from vg deconstruct output.  A record is considered nested if it has a PS
INFO field, and the variant from that field is present in the VCF
"""

import os, sys

if len(sys.argv) != 3:
    sys.stderr.write('usage: {} <vcf> <max_allele_len> \n'.format(sys.argv[0]))
    sys.exit(1)

vcf_path = sys.argv[1]
max_len = int(sys.argv[2])

def get_parent(toks):
    info = toks[7]
    itoks = info.split(';')
    parent = None
    for itok in itoks:
        if itok.startswith('PS='):
            return itok[3:]
    return None

# pass one: find big alleles and keep track of parents
too_big = set()
parents = {}
with open(vcf_path, 'r') as vcf_file:
    for line in vcf_file:
        if len(line) > 5 and not line.startswith('#'):
            toks = line.split('\t')
            ref_len = len(toks[3])
            if ref_len > max_len or max([len(x) for x in toks[4].split(',')]) > max_len:
                too_big.add(toks[2])
            parent = get_parent(toks)
            if parent:
                parents[toks[2]] = parent

# make sure if an allele is too big, its parents get flagged for deletion
big_parents = set()
for big_site in too_big:
    parent = parents[big_site] if big_site in parents else None
    while parent:
        big_parents.add(parent)
        parent = parents[parent] if parent in parents else None
too_big = too_big.union(big_parents)

# pass two: keep only stuff that is either top level or whose immediate parent has been deleted
sys.stderr.write("[strip-nested.py] Found {} sites with an allele > {} (or that are parents of such sites)\n".format(len(too_big), max_len))

# pass two: filter out records whose parents are in the file
f_count = 0
t_count = 0
with open(vcf_path, 'r') as vcf_file:
    for line in vcf_file:
        filter = False
        if len(line) > 5 and not line.startswith('#'):
            toks = line.split('\t')
            name = toks[2]
            if name in too_big:
                filter = True
            else:
                parent = get_parent(toks)
                if parent is not None and parent not in too_big:
                    filter = True
            if filter:
                f_count += 1
            t_count += 1
        if not filter:
            sys.stdout.write(line)                
sys.stderr.write("[strip-nested.py] Filtered {}/{} records\n".format(f_count, t_count))

            
        
