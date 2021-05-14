# Construct a pangenome graph, splitting on chromosomes

#!/bin/bash

# I/O options
TRUtH=""
CALLS=""
OUTDIR=""
TRUTH_SAMPLE=""
CALLS_SAMPLE=""
REFERENCE="s3://vg-k8s/vgamb/wg/cactus/GRCh38-f1-90/fasta-gap-masked/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] -t <TRUTH> -c <CALLS> -o <OUTDIR> -a <TRUTH-SAMPLE> -b <CALLS-SAMPLE>\n";
	 printf "Options:\n"a
	 printf "   -t TRUTH        truth vcf\n"
	 printf "   -c CALLS        calls vcf\n"
	 printf "   -o OUTDIR       all output here\n"
	 printf "   -a TRUTH-SAMPLE sample for truth\n"
	 printf "   -b CALLS-SAMPLE sample for calls\n"
    exit 1
}

while getopts "t:c:o:a:b:" o; do
    case "${o}" in
        t)
            TRUTH=${OPTARG}
            ;;
        c)
            CALLS=${OPTARG}
            ;;
        o)
            OUTDIR=${OPTARG}
            ;;		  
		  a)
				TRUTH_SAMPLE=${OPTARG}
				;;
		  b)
				CALLS_SAMPLE=${OPTARG}
				;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

# Check required options
if [[ $TRUTH == "" ]]; then
	 printf "truth must be specified with -t\n"
	 usage
elif [[ $CALLS == "" ]]; then
	 printf "Seqfile must be specified with -s\n"
       usage
elif [[ $OUTDIR == "" ]]; then
	 printf "Outdir must be specified with -o\n"
       usage	 
elif [[ $TRUTH_SAMPLE == "" ]]; then
	 printf "Truth sample must be specified with -a\n"
	 usage
elif [[ $CALLS_SAMPLE == "" ]]; then
	 printf "Calls sample must be specified with -b\n"
	 usage
fi

set -ex

# make the vcfeval template
LOCAL_REFERENCE="$(basename $REFERENCE)"
TEMPLATE="${LOCAL_REFERENCE}.sdf"

if [ ! -e "$TEMPLATE" ]; then
	 aws s3 cp $REFERENCE $LOCAL_REFERENCE
	 rtg format $LOCAL_REFERENCE -o $TEMPLATE
fi

# download the beds
HARD_BED="hard.bed"
EASY_BED="easy.bed"

if [ ! -f "$HARD_BED" ]; then
	 wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/union/GRCh38_alldifficultregions.bed.gz
	 gzip -dc GRCh38_alldifficultregions.bed.gz > $HARD_BED
fi

if [ ! -f "$EASY_BED" ]; then
	 wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/union/GRCh38_notinalldifficultregions.bed.gz
	 gzip -dc GRCh38_notinalldifficultregions.bed.gz > $EASY_BED
fi

# download the truth
LOCAL_TRUTH="truth.vcf.gz"
LOCAL_TRUTH_SAMPLE="truth_sample.vcf.gz"

if [ ! -f "$LOCAL_TRUTH_SAMPLE" ]; then
	 aws s3 cp $TRUTH $LOCAL_TRUTH
	 if [[ $TRUTH_SAMPLE != "" ]]; then
		  bcftools view $LOCAL_TRUTH -s $TRUTH_SAMPLE -Oz > $LOCAL_TRUTH_SAMPLE
	 else
		  mv $LOCAL_TRUTH $LOCAL_TRUTH_SAMPLE
	 fi
	 tabix -fp vcf  $LOCAL_TRUTH_SAMPLE
fi

mkdir -p $OUTDIR
# download the calls
LOCAL_CALLS="${OUTDIR}/calls.vcf.gz"
LOCAL_CALLS_SAMPLE="${OUTDIR}/calls_sample.vcf.gz"

if [ ! -f "$LOCAL_CALLS_SAMPLE" ]; then
	 aws s3 cp $CALLS $LOCAL_CALLS
	 aws s3 cp ${CALLS}.tbi ${LOCAL_CALLS}.tbi
	 bcftools view $LOCAL_CALLS -a -s $CALLS_SAMPLE | bcftools view -e 'STRLEN(REF)>2000 || STRLEN(ALT)>2000 || GT~"\." || GT="ref"' > ${OUTDIR}/temp.vcf
	 strip-nested.py ${OUTDIR}/temp.vcf | sed -e "s/GRCh38.//g" | bgzip --threads 8 > ${LOCAL_CALLS_SAMPLE}
	 tabix -fp vcf ${LOCAL_CALLS_SAMPLE}
fi

# do the whole genome eval
rm -rf ${OUTDIR}/eval-easy
rtg vcfeval -t $TEMPLATE -b $LOCAL_TRUTH_SAMPLE -c $LOCAL_CALLS_SAMPLE -o ${OUTDIR}/eval-easy -e $EASY_BED &
rm -rf ${OUTDIR}/eval-hard
rtg vcfeval -t $TEMPLATE -b $LOCAL_TRUTH_SAMPLE -c $LOCAL_CALLS_SAMPLE -o ${OUTDIR}/eval-hard -e $HARD_BED &

# do the chromosome eval
for i in `seq 22`; do
	 CHROM="chr${i}"
	 bcftools view -r $CHROM $LOCAL_TRUTH_SAMPLE -Oz > ${OUTDIR}/truth_${CHROM}.vcf.gz
	 tabix -fp vcf ${OUTDIR}/truth_${CHROM}.vcf.gz
	 bcftools view -r $CHROM $LOCAL_CALLS_SAMPLE -Oz > ${OUTDIR}/calls_${CHROM}.vcf.gz
	 tabix -fp vcf ${OUTDIR}/calls_${CHROM}.vcf.gz

	 rm -rf ${OUTDIR}/eval-easy-${CHROM}
	 rtg vcfeval -t $TEMPLATE -b ${OUTDIR}/truth_${CHROM}.vcf.gz -c ${OUTDIR}/calls_${CHROM}.vcf.gz -o ${OUTDIR}/eval-easy-${CHROM} -e $EASY_BED
	 rm -rf ${OUTDIR}/eval-hard-${CHROM}
	 rtg vcfeval -t $TEMPLATE -b ${OUTDIR}/truth_${CHROM}.vcf.gz -c ${OUTDIR}/calls_${CHROM}.vcf.gz -o ${OUTDIR}/eval-hard-${CHROM} -e $HARD_BED
done

set +x

# print the results
echo "Whole genome results"
printf "Easy\t"
cat ${OUTDIR}/eval-easy/summary.txt
printf "Hard\t"
tail -1 ${OUTDIR}/eval-hard/summary.txt
echo "Easy (by chrom)"
for i in `seq 22`; do
	 printf "chr${i}\t"
	 tail -1 ${OUTDIR}/eval-easy-chr{$i}/summary.txt
done
echo "Hard (by chrom)"
for i in `seq 22`; do
	 printf "chr${i}\t"
	 tail -1 ${OUTDIR}/eval-hard-chr${i}/summary.txt
done
