#!/bin/bash

# must end with /
BED_BUCKET=$1

# download the bed files
for i in `aws s3 ls ${BED_BUCKET} | grep mask.bed | awk '{print $4}'`; do aws s3 cp ${BED_BUCKET}${i} . ; done

# total up the masked bases for each one
rm mask-stats.tsv ; for i in *.mask.bed ; do printf "${i::-12}\t" >> mask-stats.tsv; cat $i | awk '{sum += ($3 - $2)} END {print sum}' >> mask-stats.tsv ; done


