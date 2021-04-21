#!/bin/bash

set -x

#example ./clean-vcf.sh cactus.vcf GRCh38 1000 hg38.sdf
VCF=$1
PREFIX=$2
LEN=$3
REF=$4
THREADS=20

# make the sdf
rtg format $REF -o ${REF}.sdf

# first we split out the big variants
bcftools annotate -x "INFO/AF,INFO/AC" $VCF -i "STRLEN(REF)>${LEN} || STRLEN(ALT)>${LEN}" | sed -e "s/${PREFIX}.//g" | bgzip --threads ${THREADS} > big.vcf.gz

# then we decompose the small variants
bcftools annotate -x "INFO/AF,INFO/AC" $VCF -e "STRLEN(REF)>${LEN} || STRLEN(ALT)>${LEN}" | sed -e "s/${PREFIX}.//g" | rtg vcfdecompose -i - -t ${REF}.sdf --break-indels --break-mnps -o small.dc.vcf.gz

# then we make a merged vcf
bcftools view small.dc.vcf.gz > temp.vcf
bcftools view -H big.vcf.gz >> temp.vcf

# then sort it
bcftools sort temp.vcf  | bgzip --threads ${THREADS} > merged.dc.vcf.gz
tabix -f -p vcf merged.dc.vcf.gz
rm temp.vcf

# then normalize it
bcftools norm merged.dc.vcf.gz --rm-dup all -f ${REF} --threads ${THREADS} -Oz > merged.norm.vcf.gz
tabix -f -p vcf merged.norm.vcf.gz

