#!/usr/bin/env python3
"""
Correct LV and PS tags to account for missing sites in the VCF. 

If parent in VCF:
   LV = parent LV + 1
else
   LV = 0
"""

import os, sys, gzip

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

# pass one: compute the nesting tree from PS tags
parents_by_chrom = {}
with gzip.open(vcf_path, 'r') as vcf_file:
    for line in vcf_file:
        line = line.decode('utf-8')
        if len(line) > 5 and not line.startswith('#'):
            toks = line.split('\t')
            chrom = toks[0]
            if chrom not in parents_by_chrom:
                parents_by_chrom[chrom] = {}
            name = toks[2]
            parent = get_parent(toks)
            parents_by_chrom[chrom][name] = parent

# walk up parents (from dict) counting levels
# (important: but only walk up if parent was actually seen in vcf)
def get_lv(chrom, name):
    lv = 0
    while name:
        parent = parents_by_chrom[chrom][name]
        if parent not in parents_by_chrom[chrom]:
            parent = None
        if parent:
            lv += 1
        name = parent
    return lv

# pass two: correct the LV and PS tags
with gzip.open(vcf_path, 'r') as vcf_file:
    for line in vcf_file:
        line = line.decode('utf-8')
        filter = False
        if len(line) > 5 and not line.startswith('#'):
            toks = line.split('\t')
            chrom = toks[0]
            name = toks[2]
            parent = get_parent(toks)
            lv = get_lv(chrom, name)
            if lv > 0:
                assert parent is not None
            set_parent_lv(toks, parent, lv)
            sys.stdout.write('\t'.join(toks))
        else:
            sys.stdout.write(line)

            
        
