# Construct a pangenome graph, splitting on chromosomes

#!/bin/bash

# I/O options
JOBSTORE=""
SEQFILE=""
MINIGRAPH=""
OUTPUT_BUCKET=""
OUTPUT_NAME=""
REFERENCE=""
VCF_REFERENCE=""
DECOY=""
CONFIG=""

# Workflow options
PHASE=""

# Parameters
# Leave runs of softamsked sequences unaligned if they are at least this long
# They will also be excluded from the .vg output
MASK_LEN=100000

# general toil options
TOIL_OPTS="--batchSystem mesos --provisioner aws --defaultPreemptable --betaInertia 0 --targetTime 1 --realTimeLogging"
# jobs get run on r3 clusters
TOIL_R3_OPTS="--nodeTypes r3.8xlarge:0.7 --maxNodes 25"
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
	 printf "   -r REFERENCE      Reference genome name.  This must be present in the SEQFILE.  ex: GRCh38\n"
	 printf "   -v VCF_REFERENCE  Reference genome name for VCF export (is REFERENCE by default)\n"
	 printf "   -d DECOY          Path to graph of decoy sequences\n"
	 printf "   -c CONFIG         Cactus configuration file (applied to all commands)\n"
	 printf "Workflow Options:\n"
	 printf "   -p PHASE          Resume workflow starting with given phase {map, split, align, join}\n"
	 printf "   -M MASK           Don't align softmasked sequence stretches greater than MASK. 0 to disable [default = 100000]\n"	 
    exit 1
}

while getopts "j:s:m:o:n:r:v:d:c:p:M:" o; do
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

# phase 2: divide fasta and PAF into chromosomes
if [[ $PHASE == "" || $PHASE == "map" || $PHASE == "split" ]]; then
	 cactus-graphmap-split $JOBSTORE $SEQFILE $MINIGRAPH ${OUTPUT_BUCKET}/${OUTPUT_NAME}.paf  --refContigs ${REFCONTIGS} --otherContig chrOther --reference $REFERENCE --outDir ${OUTPUT_BUCKET}/chroms-${OUTPUT_NAME} --logFile ${OUTPUT_NAME}.graphmap-split.log ${TOIL_OPTS} ${TOIL_R3_OPTS} --maskFilter ${MASK_LEN}
	 aws s3 cp  ${OUTPUT_NAME}.graphmap-split.log ${OUTPUT_BUCKET}/logs-${OUTPUT_NAME}/
fi

REFCONTIGS="${REFCONTIGS} chrOther"

# phase 3: align each chromosome with Cactus, producing output in both HAL and vg
if [[ $PHASE == "" || $PHASE == "map" || $PHASE == "split" || $PHASE == "align" ]]; then
	 aws s3 cp ${OUTPUT_BUCKET}/chroms-${OUTPUT_NAME}/chromfile.txt ./chromfile-${OUTPUT_NAME}.txt
	 aws s3 sync ${OUTPUT_BUCKET}/chroms-${OUTPUT_NAME}/seqfiles ./seqfiles-${OUTPUT_NAME} 
	 sed -i -e "s/seqfiles/seqfiles-${OUTPUT_NAME}/g" ./chromfile-${OUTPUT_NAME}.txt

	 cactus-align-batch $JOBSTORE ./chromfile-${OUTPUT_NAME}.txt ${OUTPUT_BUCKET}/align-batch-${OUTPUT_NAME} --alignCores 32 --alignOptions "--pafInput --pangenome --outVG --realTimeLogging --barMaskFilter ${MASK_LEN} --reference ${REFERENCE}" --logFile ${OUTPUT_NAME}.align.log ${TOIL_OPTS} ${TOIL_R3_OPTS}
	 aws s3 cp  ${OUTPUT_NAME}.align.log ${OUTPUT_BUCKET}/logs-${OUTPUT_NAME}/
fi

JOIN_OPTS="--clipLength ${MASK_LEN} --wlineSep . --indexCores 64"
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
do VGFILES="${VGFILES} ${OUTPUT_BUCKET}/align-batch-${OUTPUT_NAME}/${i}.vg"
done
HALFILES=""
for i in $REFCONTIGS
do HALFILES="${HALFILES} ${OUTPUT_BUCKET}/align-batch-${OUTPUT_NAME}/${i}.hal"
done
set -x

# phase 4: merge the chromosome output into whole genome HAL, GFA, VCF, XG, SNARLS and GBWT
cactus-graphmap-join $JOBSTORE --outDir $OUTPUT_BUCKET --outName $OUTPUT_NAME --reference $REFERENCE $JOIN_OPTS --vg $VGFILES --hal $HALFILES --logFile ${OUTPUT_NAME}.join.log ${TOIL_OPTS} ${TOIL_JOIN_OPTS}

aws s3 cp  ${OUTPUT_NAME}.join.log ${OUTPUT_BUCKET}/logs-${OUTPUT_NAME}/

date
