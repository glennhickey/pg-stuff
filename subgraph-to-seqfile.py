#!/usr/bin/env python3
"""
Turn a subgraph (ie extracted with vg chunk -S) into a seqfile
for re-aligning with Cactus
"""

import os, sys
import subprocess
import pysam
import argparse
import shutil
from collections import defaultdict

def main(command_line=None):                     
    parser = argparse.ArgumentParser('Make seqfile from a subgraph (ie from vg chunk -S graph.vg)')
    parser.add_argument('--subgraph', required=True,
                        help='subgraph to extract from')
    parser.add_argument('--fa-dir', required=True,
                        help='directory to write output FASTA files to')
    parser.add_argument('--threads', type=int, default=8,
                        help='number of paths to do in parallel')
    args = parser.parse_args(command_line)

    try:
        os.makedirs(args.fa_dir)
    except:
        pass
    assert os.path.isdir(args.fa_dir)

    # get the paths
    path_list = os.path.join(args.fa_dir, 'paths.list')
    if not os.path.isfile(path_list):
        with open(path_list, 'w') as path_file:
            subprocess.check_call(['vg', 'paths', '-Lx', args.subgraph], stdout=path_file)
    else:
        sys.stderr.write(f'Loading already-existing paths file {path_list}\n')

    haplotypes = set()
    # get the fasta
    with open(path_list, 'r') as path_file:
        for line in path_file:
            name = line.strip()            
            toks = name.split('#')
            sample = toks[0]
            hap = toks[1]
            haplotypes.add(f'{sample}#{hap}')

    commands = os.path.join(args.fa_dir, 'commands.txt')
    with open(commands, 'w') as commands_file:
        for hap in haplotypes:
            filename = os.path.join(args.fa_dir, hap.replace('#', '.') + '.fa.gz')
            commands_file.write(f'vg paths -x {args.subgraph} -F -Q "{hap}#" | sed -e "s/#0$//g" | bgzip > {filename}\n')
            print('{}\t{}'.format(hap.replace('#', '.'), filename))

    subprocess.check_call(f'cat {commands} | parallel -j {args.threads} {{}}', shell=True) 

            
if __name__ == '__main__':
    main()
