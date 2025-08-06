#!/usr/bin/env python3
"""
Download some bed files of hprc annotations
"""

import os, sys
import subprocess
import argparse
import shutil
import gzip
from collections import defaultdict

rm_idx_url='https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/repeat_masker/repeat_masker_bed_pre_release_v0.2.index.csv'
censat_idx_url='https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/censat/censat_centromeres_pre_release_v0.3.index.csv'
sd_idx_url='https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/segdups/segdups_v1.0.csv'
cat_idx_url='https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/annotation/cat/cat_genes_v1.index.csv'

# these were downloaded from the table browser (though would be really nice to have a url)
hg38_sd_path='/private/home/ghickey/dev/work/nesting-july/segdupes/grch38-segdupes-ucsc.bed'
hg38_rm_path='/private/home/ghickey/dev/work/nesting-july/segdupes/grch38-rm-ucsc.bed'
hg38_genes_path='/private/home/ghickey/dev/work/nesting-july/segdupes/grch38-gencode-ucsc.bed'

def download_from_table(table_url, column, threads, out_annot_path, gff=False):

    if os.path.isfile(out_annot_path):
        sys.stderr.write(f'{out_annot_path} exists, skipping re-gen\n')
        return
    
    subprocess.check_call(['wget', table_url, '-O', os.path.basename(table_url)])
    subprocess.check_call(f'grep s3 {os.path.basename(table_url)} | dos2unix | awk -F "," \'{{print ${column}}}\' | parallel -j {threads} "aws s3 cp {{}} ."', shell=True)

    subprocess.check_call(['rm', '-f', out_annot_path])
    with open(os.path.basename(table_url)) as table_file:
        for line in table_file:
            if 's3' in line:
                annot_path = os.path.basename(line.split(',')[column-1].strip())
                gff_cmd = ' | cut -f1,4,5' if gff else ''
                subprocess.check_call(f'cat {annot_path} {gff_cmd} >> {out_annot_path}', shell=True)
                subprocess.check_call(['rm', '-f', annot_path])

    subprocess.check_call(f'bedtools sort -i {out_annot_path} > {out_annot_path}.tmp', shell=True)
    subprocess.check_call(['mv', out_annot_path + '.tmp', out_annot_path])
    
def add_local(local_path, annot_path, out_path, prefix, merge=False, max_col=-1):
    with open(out_path, 'w') as out_file:
        with open(local_path, 'r') as local_file:
            for line in local_file:
                capped_line = '\t'.join(line.strip().split()[:max_col])
                out_file.write(f'{prefix}{capped_line}\n')
        with open(annot_path, 'r') as annot_file:
            for line in annot_file:
                capped_line = '\t'.join(line.strip().split()[:max_col])
                out_file.write(capped_line + '\n')
    merge_cmd = ' | bedtools merge' if merge else ''
    subprocess.check_call(f'bedtools sort -i {out_path} {merge_cmd} > {out_path}.tmp', shell=True)
    subprocess.check_call(['mv', out_path + '.tmp', out_path])

def main(command_line=None):                     
    parser = argparse.ArgumentParser('Download some bed files of hprc. Make sure you have aws cli tools')
    parser.add_argument('--threads', type=int, default=8,
                        help='Threads for parallel downloads')
    options = parser.parse_args(command_line)

    download_from_table(cat_idx_url, 4, options.threads, 'hprc-v2-genes.bed', gff=True)
    add_local(hg38_genes_path, 'hprc-v2-sd.bed', 'hprc-v2-sd-grch38.bed', 'GRCh38#0', True, 3)
    
    download_from_table(rm_idx_url, 4, options.threads, 'hprc-v2-rm.bed')
    add_local('/home/hick, 'hprc-v2-rm.bed', 'hprc-v2-rm-grch38.bed', 'GRCh38#0', False, 6)

    download_from_table(sd_idx_url, 4, options.threads, 'hprc-v2-sd.bed')
    add_local(hg38_rm_path, 'hprc-v2-sd.bed', 'hprc-v2-sd-grch38.bed', 'GRCh38#0', True, 3)
    

    
if __name__ == '__main__':
    main()
