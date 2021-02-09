#!/usr/bin/env python3

"""convenience script to get the BED output of the preprocessor all together into one file
that can be used for clipping
"""

import os
from argparse import ArgumentParser
import copy
import subprocess
import sys

def main():
    parser = ArgumentParser()
    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("workDir", help = "Working directory")
    parser.add_argument("ref", help = "Reference genome name to ignore")

    options = parser.parse_args()

    # read the seqfile
    bed_files = {}
    with open(options.seqFile, 'r') as seqfile:
        for line in seqfile:
            toks = line.strip().split()
            if not line.strip().startswith('(') and len(toks) == 2 and toks[0] != options.ref and toks[0] != '_MINIGRAPH_':
                bed_files[toks[0]] = toks[1]

    # download the bed files
    if not os.path.isdir(options.workDir):
        os.makedirs(options.workDir)
    for event, fa_path in bed_files.items():
        bed_path = fa_path + '.mask.bed'
        local_path = os.path.join(options.workDir, event + '.bed')
        cp_cmd = ['cp', bed_path, local_path]
        if  bed_path.startswith('s3://'):
            cp_cmd = ['aws', 's3'] + cp_cmd
        subprocess.check_call(cp_cmd)

    # merge the bed files, prepending the name
    for event, fa_path in bed_files.items():
        local_path = os.path.join(options.workDir, event + '.bed')
        with open(local_path, 'r') as local_file:
            for line in local_file:
                sys.stdout.write('{}.{}'.format(event, line))
        subprocess.check_call(['rm', local_path])

    subprocess.call(['rm', '-r', options.workDir])

if __name__ == "__main__":
    main()
