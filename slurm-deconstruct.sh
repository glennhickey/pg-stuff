#!/bin/bash
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
    if [[ $BASE != "chrY.vg" && $BASE != "chrM.vg" && $BASE != "chrOther.vg" && $BASE != "chrEBV.vg" ]]; then
	BASE=${BASE::-3}
	VCF=${OUTPUT_DIR}/${BASE}.${REF}.${L}.vcf.gz
	FASTA=${OUTPUT_DIR}/${BASE}.${REF}.${L}.fa
	sbatch -W ./slurm-deconstruct-one.sh $VG "-n -L ${L} -P ${REF} -f ${FASTA}" $VCF &
	#./slurm-deconstruct-one.sh $VG "-n -L ${L} -P ${REF} -f ${FASTA}" $VCF &
    fi
done

wait

cat ${OUTPUT_DIR}/*.fa > ${OUTPUT_FA}
cat ${OUTPUT_DIR}/*.fa.nesting.tsv > ${OUTPUT_FA}.nesting.tsv

bcftools concat ${OUTPUT_DIR}/*.vcf.gz | bgzip > ${OUTPUT_NAME}
tabix -fp vcf ${OUTPUT_NAME}
rm -rf ${OUTPUT_DIR}
    
 
