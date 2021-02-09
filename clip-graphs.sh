#!/bin/bash

NPROC=4

ls *.clip | parallel -j $NPROC "vg convert -f {} -Q hg38 -B -w . | pigz -p $NPROC > {}.gfa.gz"

# for a CHM13 reference:
ls *.clip | parallel -j $NPROC "vg convert -f {} -Q CHM13 -B -w . | pigz -p $NPROC > {}.gfa.gz"

