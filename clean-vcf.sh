#!/bin/bash

set -xe

#example ./clean-vcf.sh cactus.vcf.gz 1000 hg38.fa
VCF=$1
LEN=$2
REF=$3
THREADS=8

# make the sdf out of the fasta reference
rm -rf ${REF}.sdf
rtg format $REF -o ${REF}.sdf

# filter out big alleles.  also remove info fields that will be corrupted by decompose, and strip chrom prefixes
RES1=${VCF::-7}.decomposed.max${LEN}.temp1.vcf.gz
bcftools annotate -x "INFO/AF,INFO/AC" $VCF -e "STRLEN(REF)>${LEN} || STRLEN(ALT)>${LEN}" | sed -e "s/GRCh38.//g" -e "s/CHM13.//g" | bgzip --threads $THREADS > $RES1
tabix -fp vcf $RES1

# strip out nested variants
RES2=${VCF::-7}.decomposed.max${LEN}.temp2.vcf.gz
./strip-nested.py $RES1 | bgzip --threads $THREADS > $RES2
tabix -fp vcf $RES2

# then decompose the remaining variants
RES3=${VCF::-7}.decomposed.max${LEN}.temp3.vcf.gz
rm -rf ${RES3} ${RES3}.tbi
rtg vcfdecompose -i $RES2 -t ${REF}.sdf --break-indels --break-mnps -o ${RES3}

# then we normalize and sort
RES4=${VCF::-7}.decomposed.max${LEN}.vcf.gz
bcftools norm --rm-dup all -f $REF ${RES3} | bcftools sort - | bgzip --threads $THREADS > $RES4
tabix -fp vcf $RES4


