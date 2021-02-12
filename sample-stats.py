#!/bin/python

import sys

tot_sam_len = {}
chrom_sam_lens = {}
samples = set()
for i, arg in enumerate(sys.argv[1:]):
         chrom_sam_len = {}
         with open(arg, 'r') as paths_file:
                  for line in paths_file:
                           toks = line.strip().split('\t')
                           sample = toks[0]
                           path_len = int(toks[1])
                           st = sample.split('.')
                           name = '.'.join(st[:-1])
                           if name not in tot_sam_len:
                                    tot_sam_len[name] = 0
                           tot_sam_len[name] += path_len
                           if name not in chrom_sam_len:
                                    chrom_sam_len[name] = 0
                           chrom_sam_len[name] += path_len
                           samples.add(name)
         chrom_sam_lens[arg] = chrom_sam_len

for sample in samples:
         sys.stdout.write(sample)
         for chrom in chrom_sam_lens.keys():
                  l = chrom_sam_lens[chrom][sample] if sample in chrom_sam_lens[chrom] else 0
                  sys.stdout.write("\t{}".format(l))
         sys.stdout.write("\t{}\n".format(tot_sam_len[sample]))
                                   
