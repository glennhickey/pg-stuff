#!/bin/bash

#
# Run deconstruct in parallel on slurm, using nesting options
# Note that this script is a bit brittle in that it takes input using
# wildcard to select all vg files in a directory (and relies on chr1.vg existing to get sample set from)
#

set -ex

VG_DIR=$1
REF=$2
L=$3
OUTPUT_DIR=$4
OUTPUT_NAME=$5
OUTPUT_FA=${OUTPUT_NAME::-7}.fa

mkdir -p $OUTPUT_DIR
rm -rf ${OUTPUT_DIR}/*.vcf*

for VG in ${VG_DIR}/*.vg; do
    BASE=$(basename $VG)
    if [[ $BASE != "chrEBV.vg" ]]; then
	BASE=${BASE::-3}
	VCF=${OUTPUT_DIR}/${BASE}.${REF}.${L}.vcf.gz
	FASTA=${OUTPUT_DIR}/${BASE}.${REF}.${L}.fa
	sbatch -W ./slurm-deconstruct-one.sh $VG "-R -n -L ${L} -P ${REF} -f ${FASTA}" $VCF &
	#./slurm-deconstruct-one.sh $VG "-n -L ${L} -P ${REF} -f ${FASTA}" $VCF &
    fi
done

wait

cat ${OUTPUT_DIR}/*.fa > ${OUTPUT_FA}
cat ${OUTPUT_DIR}/*.fa.nesting.tsv > ${OUTPUT_FA}.nesting.tsv

for VG in ${VG_DIR}/*.vg; do
    BASE=$(basename $VG)
    if [[ $BASE == "chrX.vg" || $BASE == "chrY.vg" || $BASE == "chrM.vg" || $BASE == "chrOther.vg" || $BASE == "chrEBV.vg" ]]; then
	BASE=${BASE::-3}
	VCF=${OUTPUT_DIR}/${BASE}.${REF}.${L}.vcf.gz
	CHR1_VCF=${OUTPUT_DIR}/chr1.${REF}.${L}.vcf.gz
	# get the samples from chr1
	bcftools query -l ${CHR1_VCF} | sort > ${OUTPUT_NAME}.all-samples
	# get the samples from the vcf
	bcftools query -l ${VCF}  | sort > ${OUTPUT_NAME}.chrom.samples
	if [[ $(diff ${OUTPUT_NAME}.all-samples ${OUTPUT_NAME}.chrom.samples) != 0 ]]; then
	    # samples that are missing from this chromosome
	    comm -32 ${OUTPUT_NAME}.all-samples ${OUTPUT_NAME}.chrom.samples | awk '{print $1}' > ${OUTPUT_NAME}.missing-samples
	    bcftools view -h ${CHR1_VCF} -S ${OUTPUT_NAME}.missing-samples | bgzip > ${OUTPUT_NAME}.missing-header.vcf.gz
	    tabix -fp vcf  ${OUTPUT_NAME}.missing-header.vcf.gz
	    # add these samples to the header
	    bcftools merge ${VCF} ${OUTPUT_NAME}.missing-header.vcf.gz -Oz > ${VCF}.merge.vcf.gz
	    tabix -fp vcf  ${VCF}.merge.vcf.gz
	    # finally, sort these samples to be same as chr1
	    bcftools view ${VCF}.merge.vcf.gz -S ${OUTPUT_NAME}.all-samples -Oz > ${OUTPUT_DIR}/${BASE}.${REF}.${L}.vcf.gz
	    tabix -fp vcf ${OUTPUT_DIR}/${BASE}.${REF}.${L}.vcf.gz
	    # remove temp sutff since we're using wildcards below
	    rm -f  ${VCF}.merge.vcf.gz ${OUTPUT_NAME}.missing-header.vcf.gz
	fi	
    fi
done

bcftools concat ${OUTPUT_DIR}/*.vcf.gz | bgzip > ${OUTPUT_NAME}
tabix -fp vcf ${OUTPUT_NAME}
rm -rf ${OUTPUT_DIR}
    
 
