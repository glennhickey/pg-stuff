#!/bin/bash

# must end with /
CHROM_BUCKET=$1

NPROC=8

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

echo "import sys
totals = {}
for line in sys.stdin:
    toks = line.strip().split()
    if toks:
       name = '.'.join(toks[0].split('.')[:-1])
       unique = int(toks[1]) - int(toks[2])
       if name not in totals:
       	  totals[name] = 0
       totals[name] += unique
for key in totals.keys():
    print('{}\t{}'.format(key, totals[key]))" > unique_cov.py
       
# get path length for each one
ls *.vg | parallel -j $NPROC "vg paths -Ev {} > {}.paths.tsv"

# get the basic stats
ls *.vg | parallel -j $NPROC "vg stats -lz {} > {}.stats.tsv"

# get the path coverage
ls *.vg | parallel -j $NPROC "vg paths -cv {} > {}.coverage.tsv"

# get the binned path coverage
#ls *.vg | parallel -j $NPROC "vg depth {} -b 1000 > {}.coverage.bed"

# append some path stats
for i in *.vg ; do printf "paths\t" >> ${i}.stats.tsv ; cat ${i}.paths.tsv | wc -l >> ${i}.stats.tsv ; done

for i in *.vg ; do printf "contigs\t" >> ${i}.stats.tsv ; cat ${i}.paths.tsv | awk '{print $1}' | while read x ; do echo ${x%:*-*}; done | sort | uniq | wc -l >> ${i}.stats.tsv ; done

for i in *.vg ; do printf "samples\t" >> ${i}.stats.tsv ; cat ${i}.paths.tsv | awk '{print $1}' | while read x ; do echo ${x%.*}; done | sort | uniq | wc -l >> ${i}.stats.tsv ; done

for i in *.vg ; do printf "path-length\t" >> ${i}.stats.tsv ; cat ${i}.paths.tsv | python tot.py >> ${i}.stats.tsv ; done

for i in *.vg ; do printf "hg38-length\t" >> ${i}.stats.tsv ; cat ${i}.paths.tsv | grep hg38 | awk '{print $2}' >> ${i}.stats.tsv ; done

for i in *.vg ; do cat ${i}.coverage.tsv | python unique_cov.py | grep -v hg38 > ${i}.unique.tsv; done
for i in *.vg ; do python barchart.py ${i}.unique.tsv --save ${i}.unique.png --x_sideways --width 16 --height 16  --title ${i}.unique; done

