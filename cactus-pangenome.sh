# Construct a pangenome graph, splitting on chromosomes

#!/bin/bash

# I/O options
JOBSTORE=""
SEQFILE=""
MINIGRAPH=""
OUTPUT_BUCKET=""
OUTPUT_NAME=""
ALIGN_NAME=""
REFERENCE=""
VCF_REFERENCE=""
DECOY=""
CONFIG=""
GAP_MASK=""
CHM13_Y=""
NORMALIZE_ITERATIONS="0"

# Workflow options
PHASE=""

# Parameters
# Leave runs of softamsked sequences unaligned if they are at least this long
# They will also be excluded from the .vg output
MASK_LEN=100000

# general toil options
TOIL_OPTS="--batchSystem mesos --provisioner aws --defaultPreemptable --betaInertia 0 --targetTime 1 --realTimeLogging"
# cactus jobs get run on trusty old r3 clusters 
TOIL_R3_OPTS="--nodeTypes r3.8xlarge:0.63 --maxNodes 25"
TOIL_R4_OPTS="--nodeTypes r5.8xlarge:1.26 --maxNodes 25 --nodeStorage 1000"
# except join, which needs a little more RAM for the whole-genome indexing
TOIL_JOIN_OPTS="--nodeTypes r5.16xlarge --maxNodes 1 --nodeStorage 2000"

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [OPTIONS] -j <JOBSTORE> -s <SEQFILE> -m <MINIGRAPH> -o <OUTPUT-BUCKET> -n <OUTPUT-NAME> -r <REFERENCE> \n"
	 printf "I/O Options:\n"
	 printf "   -j JOBSTORE       Use the given jobstore.  ex: aws:us-west-2:my-job-store\n"
	 printf "   -s SEQFILE        Cactus input seqfile.  Ideally, preprocessed.  Must be local (not S3) \n"
	 printf "   -m MINIGRAPH      Use this minigraph. ex: ftp://ftp.dfci.harvard.edu/pub/hli/minigraph/HPRC-f1/GRCh38-f1-90.gfa.gz \n"	 
	 printf "   -o OUTPUT         Output bucket.  ex: s3://cactus-output/GRCh38-pangenome\n"
	 printf "   -n NAME           Output name.  All output files will be prefixed with this name\n"
	 printf "   -a ALIGN-NAME     Output name for cactus-align and everything after (by default, same as -n)\n"
	 printf "   -r REFERENCE      Reference genome name.  This must be present in the SEQFILE.  ex: GRCh38\n"
	 printf "   -v VCF_REFERENCE  Reference genome name for VCF export (is REFERENCE by default)\n"
	 printf "   -d DECOY          Path to graph of decoy sequences\n"
	 printf "   -c CONFIG         Cactus configuration file (applied to all commands)\n"
	 printf "Workflow Options:\n"
	 printf "   -p PHASE          Resume workflow starting with given phase {map, mask, split, align, join}\n"
	 printf "   -M MASK           Don't align softmasked sequence stretches greater than MASK. 0 to disable [default = 100000]\n"
	 printf "   -g                Run gap-masking step to prevent bar for handling large minimizer gaps (clumsy but improves precision)\n"
	 printf "   -y                Assume CHM13 has chrY\n"
	 printf "   -N ITERATIONS     Normalize N interations with vg\n"
    exit 1
}

while getopts "j:s:m:o:n:a:r:v:d:c:p:M:gyN:" o; do
    case "${o}" in
        j)
            JOBSTORE=${OPTARG}
            ;;
        s)
            SEQFILE=${OPTARG}
            ;;
        m)
            MINIGRAPH=${OPTARG}
            ;;		  
		  o)
				OUTPUT_BUCKET=${OPTARG}
				;;
		  n)
				OUTPUT_NAME=${OPTARG}
				;;
		  a)
				ALIGN_NAME=${OPTARG}
				;;		  
		  r)
				REFERENCE=${OPTARG}
				;;
		  v)
				VCF_REFERENCE=${OPTARG}
				;;
		  d)
				DECOY=${OPTARG}
				;;
		  c)
				CONFIG=${OPTARG}
				;;		  
		  p)
				PHASE=${OPTARG}
				;;
		  M)
				MASK_LEN=${OPTARG}
				;;
		  g)
				GAP_MASK="1"
				;;
		  y)
				CHM13_Y="1"
				;;
		  N)
				NORMALIZE_ITERATIONS=${OPTARG}
				;;
		
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

# Check required options
if [[ $JOBSTORE == "" ]]; then
	 printf "Jobstore must be specified with -j\n"
	 usage
elif [[ $SEQFILE == "" ]]; then
	 printf "Seqfile must be specified with -s\n"
	 usage
elif [[ $MINIGRAPH == "" ]]; then
	 printf "Minigraph must be specified with -m\n"
	 usage	 
elif [[ $OUTPUT_BUCKET == "" ]]; then
	 printf "Output bucket must be specified with -o\n"
	 usage
elif [[ $OUTPUT_NAME == "" ]]; then
	 printf "Output name must be specified with -n\n"
	 usage
elif [[ $REFERENCE == "" ]]; then
	 printf "Reference must be specified with -r\n"
	 usage
fi

if [[ $REFERENCE == "CHM13" ]]; then
	 REFCONTIGS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrM"
	 if [[ $CHM13_Y == "1" ]]; then
	 	REFCONTIGS="${REFCONTIGS} chrY"
	 fi
else
	 REFCONTIGS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
fi

if [[ $CONFIG != "" ]]; then
	 TOIL_OPTS="${TOIL_OPTS} --configFile $CONFIG"
fi

if [[ $MASK_LEN == "0" ]]; then
	 MASK_LEN=4000000000
fi

date

set -ex

# phase 1: map contigs to minigraph
if [[ $PHASE == "" || $PHASE == "map" ]]; then
	 cactus-graphmap $JOBSTORE $SEQFILE $MINIGRAPH ${OUTPUT_BUCKET}/${OUTPUT_NAME}.paf --outputFasta ${OUTPUT_BUCKET}/${OUTPUT_NAME}.gfa.fa  --logFile ${OUTPUT_NAME}.graphmap.log ${TOIL_OPTS} ${TOIL_R3_OPTS} --maskFilter ${MASK_LEN}
	 aws s3 cp  ${OUTPUT_NAME}.graphmap.log ${OUTPUT_BUCKET}/logs-${OUTPUT_NAME}/
	 aws s3 cp $SEQFILE ${OUTPUT_BUCKET}/
fi

# phase 2: mask coverage gaps (so bar doesn't try to realign them)
if [[ $GAP_MASK == "1" ]]; then
	 MASK_SEQFILE=${SEQFILE}.${OUTPUT_NAME}.mask
	 if [[ $PHASE == "" || $PHASE == "mask" || $PHASE == "map" ]]; then
		  cat $SEQFILE | tail -n +2 | awk -F"\t |/" '{print $1, $NF}' | sed -e 's/s3://g' -e 's/https://g' -e 's/http://g' | awk -v obucket=${OUTPUT_BUCKET} -v oname=${OUTPUT_NAME} '{print $1 "\t" obucket "/fa-masked-" oname "/" $2}' > $MASK_SEQFILE
		  cactus-preprocess $JOBSTORE $SEQFILE $MASK_SEQFILE  --realTimeLogging --logFile ${OUTPUT_NAME}.gapmask.log ${TOIL_OPTS} ${TOIL_R3_OPTS} --maskFile ${OUTPUT_BUCKET}/${OUTPUT_NAME}.paf --minLength ${MASK_LEN}
	 fi
	 SEQFILE=${MASK_SEQFILE}
fi

# phase 3: divide fasta and PAF into chromosomes
if [[ $PHASE == "" || $PHASE == "mask" || $PHASE == "map" || $PHASE == "split" ]]; then
	 cactus-graphmap-split $JOBSTORE $SEQFILE $MINIGRAPH ${OUTPUT_BUCKET}/${OUTPUT_NAME}.paf  --refContigs ${REFCONTIGS} --otherContig chrOther --reference $REFERENCE --outDir ${OUTPUT_BUCKET}/chroms-${OUTPUT_NAME} --logFile ${OUTPUT_NAME}.graphmap-split.log ${TOIL_OPTS} ${TOIL_R3_OPTS} --maskFilter ${MASK_LEN}
	 aws s3 cp  ${OUTPUT_NAME}.graphmap-split.log ${OUTPUT_BUCKET}/logs-${OUTPUT_NAME}/
fi

if [[ $REFERENCE == "GRCh38" ]]; then
    REFCONTIGS="${REFCONTIGS} chrOther"
fi

if [[ $ALIGN_NAME == "" ]]; then
	 ALIGN_NAME=${OUTPUT_NAME}
fi

# phase 4: align each chromosome with Cactus, producing output in both HAL and vg
if [[ $PHASE == "" || $PHASE == "mask" || $PHASE == "map" || $PHASE == "split" || $PHASE == "align" ]]; then
	 aws s3 cp ${OUTPUT_BUCKET}/chroms-${OUTPUT_NAME}/chromfile.txt ./chromfile-${ALIGN_NAME}.txt
	 aws s3 sync ${OUTPUT_BUCKET}/chroms-${OUTPUT_NAME}/seqfiles ./seqfiles-${ALIGN_NAME} 
	 sed -i -e "s/seqfiles/seqfiles-${ALIGN_NAME}/g" ./chromfile-${ALIGN_NAME}.txt
	 if [[ $REFERENCE == "CHM13" ]]; then
	     grep -v ^chrOther ./chromfile-${ALIGN_NAME}.txt > ./chromfile-${ALIGN_NAME}.txt.temp
	     mv ./chromfile-${ALIGN_NAME}.txt.temp ./chromfile-${ALIGN_NAME}.txt
	 fi
	 cactus-align-batch $JOBSTORE ./chromfile-${ALIGN_NAME}.txt ${OUTPUT_BUCKET}/align-batch-${ALIGN_NAME} --alignCores 32 --alignOptions "--pafInput --pangenome --outVG --realTimeLogging --barMaskFilter ${MASK_LEN} --reference ${REFERENCE} --retryCount 0" --logFile ${ALIGN_NAME}.align.log ${TOIL_OPTS} ${TOIL_R4_OPTS}
	 aws s3 cp  ${ALIGN_NAME}.align.log ${OUTPUT_BUCKET}/logs-${ALIGN_NAME}/
fi

JOIN_OPTS="--clipLength ${MASK_LEN} --wlineSep . --indexCores 63 --normalizeIterations ${NORMALIZE_ITERATIONS}"
if [[ $DECOY != "" ]]; then
	 JOIN_OPTS="--decoyGraph ${DECOY} ${JOIN_OPTS}"
fi
if [[ $REFERENCE == "CHM13" ]]; then
	 JOIN_OPTS="--rename GRCh38>GRCh38.0 ${JOIN_OPTS}"
else
	 JOIN_OPTS="--rename CHM13>CHM13.0 ${JOIN_OPTS}"
fi

if [[ $VCF_REFERENCE != "" ]]; then
    JOIN_OPTS="--vcfReference $VCF_REFERENCE ${JOIN_OPTS}"
fi

set +x
VGFILES=""
for i in $REFCONTIGS
do VGFILES="${VGFILES} ${OUTPUT_BUCKET}/align-batch-${ALIGN_NAME}/${i}.vg"
done
HALFILES=""
for i in $REFCONTIGS
do HALFILES="${HALFILES} ${OUTPUT_BUCKET}/align-batch-${ALIGN_NAME}/${i}.hal"
done
set -x

# phase 5: merge the chromosome output into whole genome HAL, GFA, VCF, XG, SNARLS and GBWT
cactus-graphmap-join $JOBSTORE --outDir $OUTPUT_BUCKET --outName $ALIGN_NAME --reference $REFERENCE $JOIN_OPTS --vg $VGFILES --hal $HALFILES --logFile ${ALIGN_NAME}.join.log ${TOIL_OPTS} ${TOIL_JOIN_OPTS}

aws s3 cp  ${ALIGN_NAME}.join.log ${OUTPUT_BUCKET}/logs-${ALIGN_NAME}/

date
