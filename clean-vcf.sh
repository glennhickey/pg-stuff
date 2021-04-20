#!/bin/bash

set -x

#example ./clean-vcf.sh cactus.vcf GRCh38 1000 hg38.sdf
VCF=$1
PREFIX=$2
LEN=$3
# comes from rtg format
REF=$4
THREADS=20

# first we split out the big variants
bcftools view $VCF -i "STRLEN(REF)>${LEN} || STRLEN(ALT)>${LEN}" | sed -e "s/${PREFIX}.//g" | bgzip --threads ${THREADS} > big.vcf.gz

# then we decompose the small variants
bcftools view $VCF -e "STRLEN(REF)>${LEN} || STRLEN(ALT)>${LEN}" | sed -e "s/${PREFIX}.//g" | rtg vcfdecompose -i - -t $REF --break-indels --break-mnps -o small.dc.vcf.gz

# then we make a merged vcf
bcftools view small.dc.vcf.gz > temp.vcf
bcftools view -H big.vcf.gz >> temp.vcf

# then sort it
bcftools sort temp.vcf -T . | bgzip --threads ${THREADS} > merged.dc.vcf.gz
tabix -f -p vcf merged.dc.vcf.gz
rm temp.vcf
