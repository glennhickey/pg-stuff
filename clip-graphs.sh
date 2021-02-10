#!/bin/bash

NPROC=4

# download the vg files
for i in `aws s3 ls ${CHROM_BUCKET} | grep vg | awk '{print $4}'`; do echo ${CHROM_BUCKET}${i} ; done | parallel -j $NPROC "aws s3 cp {} ."

# clip them
ls *.vg | parallel -j $NPROC "clip-vg {} $MASK_BED -r 'mat>2' -r 'pat>1' -r 'CHM13>CHM13.0' > {}.clip"

# for a CHM13 reference:
#ls *.vg | parallel -j $NPROC "clip-vg {} $MASK_BED -r 'mat>2' -r 'pat>1' > {}.clip"
