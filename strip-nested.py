#!/usr/bin/env python3
"""
Filter out nested records from vg deconstruct output.  A record is considered nested if it has a PS
INFO field, and the variant from that field is present in the VCF
"""

import vcf
import os, sys

if len(sys.argv) != 2:
    sys.stderr.write('usage: {} <vcf>\n'.format(sys.argv[0]))
    sys.exit(1)

vcf_path = sys.argv[1]

# pass one: store every id
vcf_ids = set()
vcf_reader = vcf.Reader(filename=vcf_path)
for record in vcf_reader:
    vcf_ids.add(record.ID)

# pass two: filter out records whose parents are in the file
f_count = 0
vcf_reader = vcf.Reader(filename=vcf_path)
vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
for record in vcf_reader:
    if "PS" not in record.INFO or record.INFO["PS"] not in vcf_ids:
        vcf_writer.write_record(record)
    else:
        f_count += 1
sys.stderr.write(f"Filtered {f_count}/{len(vcf_ids)} records\n")

            
        
