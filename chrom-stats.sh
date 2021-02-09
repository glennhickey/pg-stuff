#!/bin/bash

# must end with /
CHROM_BUCKET=$1

NPROC=4

# download the vg files
for i in `aws s3 ls ${CHROM_BUCKET} | grep vg | awk '{print $4}'`; do echo ${CHROM_BUCKET}${i} ; done | parallel -j $NPROC "aws s3 cp {} ."

# make the python script (don't trust awk)
echo "import sys
tot = 0
for line in sys.stdin:
    toks = line.strip().split()
    if toks:
       tot += int(toks[-1])
print(tot)" > tot.py

# get path length for each one
ls *.vg | parallel -j $NPROC "vg paths -Ev {} > {}.paths.tsv"

# get the basic stats
ls *.vg | parallel -j $NPROC "vg stats -lz {} > {}.stats.tsv"

# append some path stats
for i in *.vg ; do printf "paths\t" >> ${i}.stats.tsv ; cat ${i}.paths.tsv | wc -l >> ${i}.stats.tsv ; done

for i in *.vg ; do printf "contigs\t" >> ${i}.stats.tsv ; cat ${i}.paths.tsv | awk '{print $1}' | while read x ; do echo ${x%:*-*}; done | sort | uniq | wc -l >> ${i}.stats.tsv ; done

for i in *.vg ; do printf "samples\t" >> ${i}.stats.tsv ; cat ${i}.paths.tsv | awk '{print $1}' | while read x ; do echo ${x%.*}; done | sort | uniq | wc -l >> ${i}.stats.tsv ; done

for i in *.vg ; do printf "path-length\t" >> ${i}.stats.tsv ; cat ${i}.paths.tsv | python tot.py >> ${i}.stats.tsv ; done

for i in *.vg ; do printf "hg38-length\t" >> ${i}.stats.tsv ; cat ${i}.paths.tsv | grep hg38 | awk '{print $2}' >> ${i}.stats.tsv ; done
