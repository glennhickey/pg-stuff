#!/usr/bin/env python3
"""
Filter out nested records from vg deconstruct output.  A record is considered nested if it has a PS
INFO field, and the variant from that field is present in the VCF.


"""

import os, sys, gzip

if len(sys.argv) not in [2, 3, 4]:
    sys.stderr.write('strip-nested.py: keep only highest-level sites (whose ref-alleles are <= max_allele_len if specified)\n\n')
    sys.stderr.write('usage: {} <vcf> [max_ref_allele_len] [max_alt_allele_len] | bgzip > top-level.vcf.gz \n'.format(sys.argv[0]))
    sys.exit(1)

vcf_path = sys.argv[1]
max_ref_len = int(sys.argv[2]) if len(sys.argv) > 2 else sys.maxsize
max_alt_len = int(sys.argv[3]) if len(sys.argv) > 3 else sys.maxsize

if not vcf_path.endswith('vcf.gz'):
    sys.stderr.write('only .vcf.gz input supported\n')
    sys.exit(1)

def get_parent(toks):
    info = toks[7]
    itoks = info.split(b';')
    parent = None
    for itok in itoks:
        if itok.startswith(b'PS='):
            return itok[3:]
    return None

# pass one: find big alleles and keep track of parents
too_big = set()
parents = {}
ids = set()
with gzip.open(vcf_path, 'rb') as vcf_file:
    for line in vcf_file:
        if len(line) > 5 and not line.startswith(b'#'):
            toks = line.split(b'\t')
            name = toks[2]
            ref_len = len(toks[3])
            alt_len = max([len(alt) for alt in toks[3].split(',')]) if max_alt_len != sys.maxsize else -1
            parent = get_parent(toks)
            if ref_len > max_ref_len or alt_len > max_alt_len:
                too_big.add(toks[2])
            if parent:
                parents[toks[2]] = parent
            ids.add(name)

# make sure if an allele is too big, its parents get flagged for deletion
big_parents = set()
for big_site in too_big:
    parent = parents[big_site] if big_site in parents else None
    while parent:
        big_parents.add(parent)
        parent = parents[parent] if parent in parents else None
too_big = too_big.union(big_parents)
        
# pass two: keep only stuff that is either top level or whose immediate parent has been deleted
if max_ref_len != sys.maxsize:
    sys.stderr.write("[strip-nested.py] Found {} sites with ref allele length > {}".format(len(too_big), max_ref_len))
    if max_alt_len != sys.maxsize:
        sys.sderr.write(" or alt length > {}".format(max_alt_len))
    sys.stderr.write(" (or that are parents of such sites)\n")

# pass two: filter out records whose parents are in the file
f_count = 0
t_count = 0
warning_set = set()
with gzip.open(vcf_path, 'rb') as vcf_file:
    for line in vcf_file:
        filter = False
        if len(line) > 5 and not line.startswith(b'#'):
            toks = line.split(b'\t')
            name = toks[2]
            parent = get_parent(toks)
            if parent is not None and parent not in ids and parent not in warning_set:
                sys.stderr.write("PS {} not in vcf\n".format(parent.decode("utf-8")))
                warning_set.add(parent)
            if name in too_big:
                filter = True
            elif parent is not None and parent in ids and parent not in too_big:
                filter = True
        if not filter:
            sys.stdout.write(line.decode("utf-8"))
        else:
            f_count += 1
        t_count += 1
sys.stderr.write("[strip-nested.py] Found {} PS ids not present in VCF\n".format(len(warning_set)))
sys.stderr.write("[strip-nested.py] Filtered {}/{} records\n".format(f_count, t_count))

            
        
