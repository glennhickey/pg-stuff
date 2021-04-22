#!/bin/bash

set -xe

#example ./clean-vcf.sh cactus.vcf.gz GRCh38 1000 hg38.fa
VCF=$1
PREFIX=$2
LEN=$3
REF=$4
THREADS=8

# make the sdf
rm -rf ${REF}.sdf
rtg format $REF -o ${REF}.sdf

# then decompose the small variants
RES1=${VCF::-7}.decomposed.max${LEN}.not.normalized.vcf.gz
rm -rf ${RES1} ${RES1}.tbi
bcftools annotate -x "INFO/AF,INFO/AC" $VCF -e "STRLEN(REF)>${LEN} || STRLEN(ALT)>${LEN}" | sed -e "s/${PREFIX}.//g" | rtg vcfdecompose -i - -t ${REF}.sdf --break-indels --break-mnps -o ${RES1}

# then we normalize and sort
RES2=${VCF::-7}.decomposed.max${LEN}.vcf.gz
bcftools norm --rm-dup all -f $REF ${RES1} | bcftools sort - | bgzip --threads $THREADS > $RES2
tabix -fp vcf $RES2


