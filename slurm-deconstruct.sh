#!/bin/bash
set -ex

VG_DIR=$1
REF=$2
L=$3
OUTPUT_DIR=$4
OUTPUT_NAME=$5

mkdir -p $OUTPUT_DIR
rm -rf ${OUTPUT_DIR}/*.vcf*

for VG in ${VG_DIR}/*.vg; do
    if [[ $VG != "chrY.vg" ]]; then
	BASE=$(basename $VG)
	BASE=${BASE::-3}
	VCF=${OUTPUT_DIR}/${BASE}.${REF}.${L}.vcf.gz
	#$sbatch -W ./slurm-deconstruct-one.sh $VG "-n -L ${L} -P ${REF}" $VCF &
	./slurm-deconstruct-one.sh $VG "-n -L ${L} -P ${REF}" $VCF &
    fi
done

wait

bcftools concat ${OUTPUT_DIR}/*.vcf.gz | bgzip > ${OUTPUT_NAME}
tabix -fp vcf ${OUTPUT_NAME}.vcf.gz
rm -rf ${OUTPUT_DIR}
    
    
