#!/usr/bin/env python
"""
Filter out nested records from vg deconstruct output.  A record is considered nested if it has a PS
INFO field, and the variant from that field is present in the VCF
"""

import os, sys

if len(sys.argv) != 2:
    sys.stderr.write('usage: {} <vcf>\n'.format(sys.argv[0]))
    sys.exit(1)

vcf_path = sys.argv[1]

# pass one: store every id
vcf_ids = set()
with open(vcf_path, 'r') as vcf_file:
    for line in vcf_file:
        if len(line) > 5 and not line.startswith('#'):
            vcf_ids.add(line.split('\t')[2])

# pass two: filter out records whose parents are in the file
f_count = 0
with open(vcf_path, 'r') as vcf_file:
    for line in vcf_file:
        filter = False
        if len(line) > 5 and not line.startswith('#'):
            info = line.split('\t')[7]
            itoks = info.split(';')
            parent = None
            for itok in itoks:
                if itok.startswith('PS='):
                    parent = itok[3:]
                    break
            if parent and parent in vcf_ids:
                filter = True
                f_count += 1
        if not filter:
            sys.stdout.write(line)                
sys.stderr.write("Filtered {}/{} records\n".format(f_count, len(vcf_ids)))

            
        
