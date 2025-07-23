#!/bin/bash
# Job name:
#SBATCH --job-name=deconstruct
#
# Partition - This is the queue it goes in:
#SBATCH --partition=long
#
# Where to send email (optional)
#SBATCH --mail-user=glenn.hickey@gmail.com
#
# Number of nodes you need per job:
#SBATCH --nodes=1
#
# Memory needed for the jobs.  Try very hard to make this accurate.  DEFAULT = 4gb
#SBATCH --mem=250gb
#
# Number of tasks (one for each CPU desired for use case) (example):
#SBATCH --ntasks=1
#
# Processors per task:
# At least eight times the number of GPUs needed for nVidia RTX A5500
#SBATCH --cpus-per-task=32
#
# Number of GPUs, this can be in the format of "--gres=gpu:[1-8]", or "--gres=gpu:A5500:[1-8]" with the type included (optional)
#
# Standard output and error log
#SBATCH --output=deconstruct_%j.log
#
# Wall clock limit in hrs:min:sec:
#SBATCH --time=16:00:00
#
## Command(s) to run (example):

VG=$1
#make sure these come in inside quotes!
OPTS=$2
VCF=$3
	 
vg deconstruct $VG $OPTS -t 32 | bgzip > $VCF
tabix -fp vcf $VCF
