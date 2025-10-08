#!/usr/bin/env python3
"""
Get some basic VCF Stats (since bcftools stats doesn't, ex, count alleles)
"""

import os, sys
import subprocess
import pysam
import argparse
import shutil
from collections import defaultdict

    
def main(command_line=None):                     
    parser = argparse.ArgumentParser('VCF stats: print some summary info to stdout')
    parser.add_argument('--vcf', required=True,
                        help='VCF file')
    parser.add_argument('--ref', required=True,
                        help='Reference prefix')
    parser.add_argument('--tsv',
                        help='Output tsv file')
    parser.add_argument('--sv-len', type=int, default=50,
                        help='Minimum length (delta between ref/alt) to be SV [50]')
    
    args = parser.parse_args(command_line)

    num_sv_sites = 0
    num_sv_alleles = 0
    num_small_sites = 0
    num_small_alleles = 0
    num_smnp_sites = 0
    num_smnp_alleles = 0    

    vcf_file = pysam.VariantFile(args.vcf, 'rb')
    tsv_file = open(args.tsv, 'w') if args.tsv else None    
    for var in vcf_file.fetch():
        ref = var.contig.startswith(args.ref)
        lens = [len(a) for a in var.alleles]
        if 'AF' in var.info:
            af = round(max(var.info['AF']), 5)
        else:
            af = 0
        ref_len = lens[0]
        if var.alleles[0] == '.' and 'NR' in var.info:
            ref_len = len(var.info['NR'])
        alt_lens = [0]
        if len(var.alleles) > 1:
            alt_lens = lens[1:]       
        delta = max(abs(alt_len - ref_len) for alt_len in alt_lens)
        if tsv_file:
            tsv_file.write(f'{var.contig}\t{int(ref)}\t{delta}\t{len(lens)}\t{af}\n')
        if delta >= args.sv_len:
            num_sv_sites += 1
            num_sv_alleles += len(lens)
        else:
            if max(lens) < args.sv_len:
                num_small_sites += 1
                num_small_alleles += len(lens)
            else:
                num_smnp_sites += 1
                num_smnp_alleles += len(lens)

    print(f'{num_small_sites}\t{num_small_alleles}\t{num_sv_sites}\t{num_sv_alleles}\t{num_smnp_sites}\t{num_smnp_alleles}\n')
        
    vcf_file.close()
    if tsv_file:
        tsv_file.close()
            
if __name__ == '__main__':
    main()
