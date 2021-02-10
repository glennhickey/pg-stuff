#!/bin/bash

PREFIX="HG002."

ls *.vg | parallel -j $NPROC "vg paths -Q ${PREFIX} -Ev {} > {}.paths"

cat *.paths | awk '{print $1 "\t0\t" $2}'

rm -f *.paths
