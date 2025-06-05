#!/usr/bin/env python3
"""
VCF comparison with hap.py and truvari, alongiside some postproecssing functions

Requirements:  * python3 environment with pysam installed
               * bcftools, samtools and tabix
               * docker

Example:

vcfcomp.py download-q100

vcfcomp.py happy --truth CHM13v2.0_HG2-T2TQ100-V1.1_smvar.vcf.gz \
 --calls hprc-v2.0-mc-chm13.wave.vcf.gz --ref hs1.fa \
 --regions CHM13v2.0_HG2-T2TQ100-V1.1_smvar.benchmark.bed \
 --out-dir happy-hprc-v2.0-mc-chm13.HG002.wave --sample HG002 \
 --threads 8 --exclude-y --max-length 50

vcfcomp.py happy-breakdown --vcf happy-hprc-v2.0-mc-chm13.HG002.wave/happy.output.vcf.gz




"""

import os, sys
import subprocess
import pysam
import argparse
import shutil
from collections import defaultdict

default_happy_docker = 'jmcdani20/hap.py:v0.3.12'
default_happy_options = '--gender=male --pass-only --engine=vcfeval'
default_happy_max_length = 100
default_truvari_docker = 'solyris/truvari:v5.3'
default_truvari_options = '-O 0.0 -r 1000 -p 0.0 -P 0.3 -C 1000 -s 50 -S 15 --sizemax 100000 --no-ref c'

def vcf_preprocess(input_vcf,
                   output_vcf,
                   remove_y,
                   haploid_x,
                   sample,
                   max_length,
                   bi_allelic):
    """
    hap.py will fail if chrY is diploid (ie .|1).  So we need to flatten it to haploid.  Will also extract
    A sample if given, or just remove Y entirely if specified. 
    """

    output_dir = os.path.dirname(output_vcf)
    try:
        os.makedirs(output_dir)
    except:
        pass
    assert os.path.isdir(output_dir)

    tmp_vcf = output_vcf.replace('.vcf', '.c1.vcf'.format('.' + sample if sample else ''))
    view_cmd = ['bcftools', 'view', input_vcf, '-c', '1', '-Oz']
    if max_length:
        view_cmd += ['-i', 'STRLEN(REF) <= {} && STRLEN(ALT) <= {}'.format(max_length, max_length)]
    if sample:
        view_cmd += ['-as', sample]
    with open(tmp_vcf, 'w') as tmp_vcf_file:
        sys.stderr.write('Running vcf preprocessing: {}\n'.format(' '.join(view_cmd)))
        subprocess.check_call(view_cmd, stdout=tmp_vcf_file)    
    if bi_allelic:
        sys.stderr.write('Splitting into bi-allelic variants\n')
        norm_vcf = tmp_vcf.replace('.vcf', '.biallelic.vcf')
        with open(norm_vcf, 'w') as norm_vcf_file:
            subprocess.check_call(['bcftools', 'norm', '-m', '-any', tmp_vcf, '-Oz'],
                                  stdout=norm_vcf_file)
        subprocess.check_call(['tabix', '-fp', 'vcf', norm_vcf])
        os.remove(tmp_vcf)
        tmp_vcf = norm_vcf
        
    subprocess.check_call(['tabix', '-fp', 'vcf', tmp_vcf])
    
    tmp_vcf_file = pysam.VariantFile(tmp_vcf, 'rb')
    output_vcf_file = pysam.VariantFile(output_vcf, 'w', header=tmp_vcf_file.header)

    sys.stderr.write('Correcting/removing chrX,chrY,chrM GTs\n')
    for var in tmp_vcf_file.fetch():
        skip = ('chry' in var.contig.lower() and remove_y) or 'chrm' in var.contig.lower()
        flatten = not skip and ('chry' in var.contig.lower() or (haploid_x and 'chrx' in var.contig.lower()))
        if flatten:
            for sample in var.samples.values():
                gt = sample['GT']
                assert len(gt) in [1,2]
                if len(gt) == 2:
                    hap_gt = gt[0]
                    if hap_gt in ['.', None] and gt[1] not in ['.', None]:
                        hap_gt = gt[1]
                    sample['GT'] = tuple([hap_gt])
                                     
        if not skip:
            output_vcf_file.write(var)

    tmp_vcf_file.close()
    output_vcf_file.close()
    
    subprocess.check_call(['tabix', '-fp', 'vcf', output_vcf])
    os.remove(tmp_vcf)
    os.remove(tmp_vcf + '.tbi')

def happy(truth_vcf,
          calls_vcf,
          ref_fasta,
          bed_regions,
          sample,
          output_dir,
          threads,
          options,
          docker_image):
    """
    Run hap.py on two vcfs, using the provided reference and bed
    """

    try:
        os.makedirs(output_dir)
    except:
        pass
    assert os.path.isdir(output_dir)

    # put everything in the same place
    for f in [truth_vcf, truth_vcf + '.tbi', calls_vcf, calls_vcf + '.tbi', ref_fasta, ref_fasta + '.fai', bed_regions]:
        if os.path.isfile(f) and os.path.dirname(f) != output_dir:
            subprocess.check_call(['ln', '-f', os.path.abspath(f), os.path.abspath(output_dir)])

    # run hap.py with docker (and since it's python2 only, it's pretty much the only way to run it these days)
    happy_cmd = ['/opt/hap.py/bin/hap.py',
                 '/data/' + os.path.basename(truth_vcf),
                 '/data/' + os.path.basename(calls_vcf),                 
                 '-f', '/data/' + os.path.basename(bed_regions),
                 '-r', '/data/' + os.path.basename(ref_fasta),
                 '--threads', str(threads),
                 '-o', '/data/' + os.path.basename(output_dir)]
    if options:
        happy_cmd += options.split()

    docker_cmd = ['docker', 'run', '-it', '--rm',
                  '-u', '{}:{}'.format(os.getuid(), os.getgid()),
                  '-v', os.path.abspath(output_dir) + ':/data', docker_image] + happy_cmd
    print(docker_cmd)
    subprocess.check_call(docker_cmd)

    # use Parsa's convention of having a happy.output.vcf.gz
    subprocess.check_call(['ln', '-sf', os.path.abspath(os.path.join(output_dir, os.path.basename(output_dir) + '.vcf.gz')),
                           os.path.abspath(os.path.join(output_dir, 'happy-output.vcf.gz'))])
    subprocess.check_call(['ln', '-sf', os.path.abspath(os.path.join(output_dir, os.path.basename(output_dir) + '.vcf.gz.tbi')),
                           os.path.abspath(os.path.join(output_dir, 'happy-output.vcf.gz.tbi'))])    

def happy_chromsplit(happy_vcf):
    """ breakdown the happy vcf into a per-chromosome table"""

    snp_fp = defaultdict(int)
    snp_fn = defaultdict(int)
    indel_fp = defaultdict(int)
    indel_fn = defaultdict(int)

    vcf_file = pysam.VariantFile(happy_vcf, 'rb')
    contigs = set(['_total_'])
    for var in vcf_file.fetch():
        assert len(var.samples.values()) == 2
        contigs.add(var.contig)
        for sample in var.samples.values():
            bd = sample['BD']
            bvt = sample['BVT']
            if bd == 'FP' and bvt == 'SNP':
                snp_fp[var.contig] += 1
                snp_fp['_total_'] += 1
            elif bd == 'FP' and bvt == 'INDEL':
                indel_fp[var.contig] += 1
                indel_fp['_total_'] += 1
            elif bd == 'FN' and bvt == 'SNP':
                snp_fn[var.contig] += 1
                snp_fn['_total_'] += 1                
            elif bd == 'FN' and bvt == 'INDEL':
                indel_fn[var.contig] += 1
                indel_fn['_total_'] += 1                
    vcf_file.close()

    output_table = [['type'] + sorted(contigs)]
    output_table += [['ERRORS'], ['SNP-FP'], ['SNP-FN'], ['INDEL-FP'], ['INDEL-FN']]
    for contig in sorted(contigs):
        output_table[-5].append(str(snp_fp[contig] + snp_fn[contig] + indel_fp[contig] + indel_fn[contig]))
        output_table[-4].append(str(snp_fp[contig]))
        output_table[-3].append(str(snp_fn[contig]))
        output_table[-2].append(str(indel_fp[contig]))
        output_table[-1].append(str(indel_fn[contig]))
    
    return output_table

def truvari(truth_vcf,
            calls_vcf,
            ref_fasta,
            bed_regions,
            sample,
            output_dir,
            options,
            docker_image):
    """
    Run truvari on two vcfs, using the provided reference and bed
    """
    try:
        os.makedirs(output_dir)
    except:
        pass
    assert os.path.isdir(output_dir)

    # truvari needs a non-existent output directory
    tv_dir = os.path.join(output_dir, 'tv')
    if os.path.exists(tv_dir):
        shutil.rmtree(tv_dir)

    # put everything in the same place
    for f in [truth_vcf, truth_vcf + '.tbi', calls_vcf, calls_vcf + '.tbi', ref_fasta, ref_fasta + '.fai', bed_regions]:
        if os.path.isfile(f) and os.path.dirname(f) != output_dir:
            subprocess.check_call(['ln', '-f', os.path.abspath(f), os.path.abspath(output_dir)])

    # run truvari with docker (to be consistent with hap.py method, tho truvari easy to install with pip)
    truvari_cmd = ['truvari', 'bench',
                   '-c', '/data/' + os.path.basename(calls_vcf),
                   '-b', '/data/' + os.path.basename(truth_vcf),
                   '--includebed', '/data/' + os.path.basename(bed_regions),
                   '-f', '/data/' + os.path.basename(ref_fasta),
                   '-o', '/data/tv/']
    if options:
        truvari_cmd += options.split()

    subprocess.check_call(['docker', 'run', '-it', '--rm',
                           '-u', '{}:{}'.format(os.getuid(), os.getgid()),
                           '-v', os.path.abspath(output_dir) + ':/data', docker_image] + truvari_cmd)

    subprocess.check_call('mv {}/* {}'.format(tv_dir, output_dir), shell=True)
    shutil.rmtree(tv_dir)

def truvari_chromsplit(truvari_outdir):
    """ breakdown the truvari vcfs into a per-chromosome table"""

    fp = defaultdict(int)
    fn = defaultdict(int)

    fp_vcf = os.path.join(truvari_outdir, 'fp.vcf.gz')
    fp_vcf_file = pysam.VariantFile(fp_vcf, 'r')
    contigs = set(['_total_'])
    for var in fp_vcf_file.fetch():
        contigs.add(var.contig)
        fp[var.contig] += 1
        fp['_total_'] += 1
    fp_vcf_file.close()

    fn_vcf = os.path.join(truvari_outdir, 'fn.vcf.gz')
    fn_vcf_file = pysam.VariantFile(fn_vcf, 'r')
    contigs = set(['_total_'])
    for var in fn_vcf_file.fetch():
        contigs.add(var.contig)
        fn[var.contig] += 1
        fn['_total_'] += 1
    fn_vcf_file.close()

    output_table = [['type'] + sorted(contigs)]
    output_table += [['ERRORS'], ['SV-FP'], ['SV-FN']]
    for contig in sorted(contigs):
        output_table[-3].append(str(fp[contig] + fn[contig]))
        output_table[-2].append(str(fp[contig]))
        output_table[-1].append(str(fn[contig]))
    
    return output_table

def download_q100(dict_only=False):
    """
    Get all the benchmark data needed for the draft HG002 Q100 benchmarks
    """
    base_url = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/'

    name_dict = {}
    for ref in ['CHM13v2.0', 'GRCh38']:
        for vartype in ['smvar', 'stvar']:
            for ext in ['benchmark.bed', 'vcf.gz', 'vcf.gz.tbi']:
                name = '{}_HG2-T2TQ100-V1.1_{}.{}'.format(ref, vartype, ext)
                if not dict_only:
                    sys.stderr.write('Downloading {}\n'.format(name))
                    subprocess.check_call(['wget', '-q', os.path.join(base_url, name), '-O', name])
                    # hack to drop y from chm13 benchmark
                    if ref == 'CHM13v2.0' and ext == 'benchmark.bed':
                        no_y_name = name.replace('.benchmark.bed', '.no.y.benchmark.bed')
                        with open(no_y_name, 'w') as no_y_bed_file:
                            subprocess.check_call(['grep', '-v', '^chrY', name], stdout=no_y_bed_file)
                        name = no_y_name
                name_dict[(ref, vartype, ext)] = name

    if not dict_only:
        sys.stderr.write('Downlading hs1.fa.gz\n')
        subprocess.check_call(['wget', '-q', 'https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz', '-O', 'hs1.fa.gz'])

        sys.stderr.write('Downloading hg38.fa.gz\n')
        subprocess.check_call(['wget', '-q', 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz', '-O', 'hg38.fa.gz'])

    for ref in 'hs1', 'hg38':
        # note: I'm pretty sure hap.py doesn't work properly with compressed references
        if not dict_only:
            sys.stderr.write('Indexing {}.fa\n'.format(ref))
            subprocess.check_call(['gzip', '-fd', '{}.fa.gz'.format(ref)])
            subprocess.check_call(['samtools', 'faidx', '{}.fa'.format(ref)])
        name_dict[ref] = '{}.fa'.format(ref)

        #hack
        #name_dict['hg38'] = 'GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta'
    return name_dict
                                 
def eval_q100(out_dir, grch38_vcfs, chm13_vcfs, download, happy_max_length, happy_threads, eval_type):
    """ run q100 hap.py and truvari evaluation on set of input vcfs and
    tabulate the results """

    name_dict = download_q100(dict_only = not download)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    assert os.path.isdir(out_dir)

    results_table = {}

    for i, vcf in enumerate(grch38_vcfs + chm13_vcfs):
        is_grch38 = i < len(grch38_vcfs)
        happy_out_dir = os.path.join(out_dir, 'happy-' + os.path.basename(vcf.replace('.vcf.gz','')))
        # run hap.py preprocessing
        hap_calls = os.path.join(happy_out_dir, os.path.basename(vcf.replace('.vcf.gz', '.hap.vcf.gz')))
        if eval_type in ['happy', 'both']:
            vcf_preprocess(vcf,
                           hap_calls,
                           not is_grch38,
                           True,
                           'HG002',
                           happy_max_length,
                           False)        
            # run hap.py
            happy(name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'smvar', 'vcf.gz')],
                  hap_calls,
                  name_dict['hg38' if is_grch38 else 'hs1'],
                  name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'smvar', 'benchmark.bed')],
                  'HG002',
                  happy_out_dir,
                  happy_threads,
                  default_happy_options,
                  default_happy_docker)

            # run hap.py postprocessing
            results_table[os.path.basename(happy_out_dir)] = happy_chromsplit(os.path.join(happy_out_dir, 'happy-output.vcf.gz'))        

        if eval_type in ['truvari', 'both']:
            truvari_out_dir = os.path.join(out_dir, 'truvari-' + os.path.basename(vcf.replace('.vcf.gz','')))
            # run truvari preprocessing
            tru_calls = os.path.join(truvari_out_dir, os.path.basename(vcf.replace('.vcf.gz', '.tru.vcf.gz')))
            vcf_preprocess(vcf,
                           tru_calls,
                           not is_grch38,
                           True,
                           'HG002',
                           None,
                           True)
        
            # run truvari
            truvari(name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'stvar', 'vcf.gz')],
                    tru_calls,
                    name_dict['hg38' if is_grch38 else 'hs1'],
                    name_dict[('GRCh38' if is_grch38 else 'CHM13v2.0', 'stvar', 'benchmark.bed')],
                    'HG002',
                    truvari_out_dir,
                    default_truvari_options,
                    default_truvari_docker)

            # run truvari postprocessing
            results_table[os.path.basename(truvari_out_dir)] = truvari_chromsplit(truvari_out_dir)

    if eval_type in ['happy', 'both']:
        with open(os.path.join(out_dir, 'happy.tsv'), 'w') as happy_file:
            for k,v in results_table:
                if 'happy' in k:
                    for line in v:
                        happy_file.write('{}\t{}\n'.format(k, '\t'.join(line)))

    if eval_type in ['truvari', 'both']:
        with open(os.path.join(out_dir, 'truvari.tsv'), 'w') as truvari_file:
            for k,v in results_table:
                if 'truvari' in k:
                    for line in v:
                        truvari_file.write('{}\t{}\n'.format(k, '\t'.join(line)))
    
def main(command_line=None):                     
    parser = argparse.ArgumentParser('VCF comparison tools for evaluating pangenome graphs')
    subparsers = parser.add_subparsers(dest='command')
    
    happy_parser = subparsers.add_parser('happy', help='Run hap.py')
    happy_parser.add_argument('--truth', required=True,
                       help='true VCF')
    happy_parser.add_argument('--calls', required=True,
                       help='calls VCF')
    happy_parser.add_argument('--ref', required=True,
                       help='reference fasta')
    happy_parser.add_argument('--regions', required=True,
                       help='BED regions to evaluate')
    happy_parser.add_argument('--out-dir', required=True,
                       help='output directory')    
    happy_parser.add_argument('--sample',
                       help='subset to this sample')
    happy_parser.add_argument('--max-length', type=int, default=default_happy_max_length,
                       help='ignore sites with alleles longer than this')    
    happy_parser.add_argument('--docker', default=default_happy_docker,
                       help='use this docker image instead of the default')
    happy_parser.add_argument('--options', default=default_happy_options,
                       help='use these options instead of the default (surround in quotes)')
    happy_parser.add_argument('--threads', type=int, default=8,
                       help='number of threads (default=8)')
    happy_parser.add_argument('--haploid-x', action='store_true',
                              help='make sure chrX is haploid')    
    happy_parser.add_argument('--exclude-y', action='store_true',
                       help='completely ignore chrY')

    happy_breakdown_parser = subparsers.add_parser('happy-breakdown', help='Make chromosome-decomposed table of hap.py results')
    happy_breakdown_parser.add_argument('--vcf', required=True,
                                        help='happy output VCF')

    truvari_parser = subparsers.add_parser('truvari', help='Run truvari')
    truvari_parser.add_argument('--truth', required=True,
                       help='true VCF')
    truvari_parser.add_argument('--calls', required=True,
                       help='calls VCF')
    truvari_parser.add_argument('--ref', required=True,
                       help='reference fasta')
    truvari_parser.add_argument('--regions', required=True,
                       help='BED regions to evaluate')
    truvari_parser.add_argument('--out-dir', required=True,
                       help='output directory')    
    truvari_parser.add_argument('--sample',
                       help='subset to this sample')
    truvari_parser.add_argument('--docker', default=default_truvari_docker,
                                help='use this docker image instead of the default')
    truvari_parser.add_argument('--options', default=default_truvari_options,
                                help='use these options instead of the default (surround in quotes)')
    truvari_parser.add_argument('--haploid-x', action='store_true',
                              help='make sure chrX is haploid')    
    truvari_parser.add_argument('--exclude-y', action='store_true',
                              help='completely ignore chrY')

    truvari_breakdown_parser = subparsers.add_parser('truvari-breakdown', help='Make chromosome-decomposed table of truvari results')
    truvari_breakdown_parser.add_argument('--dir', required=True,
                                        help='truvari output directory')

    hg002_q100_download_parser = subparsers.add_parser('download-q100', help='Download all HG002 t2t-q100 truth set files')

    hg002_q100_eval_parser = subparsers.add_parser('eval-q100', help='Run t2t-q100 HG002 evaluation on some graphs')
    hg002_q100_eval_parser.add_argument('--download', action='store_true',
                                        help='automatically download benchmark data into current directory')
    hg002_q100_eval_parser.add_argument('--out-dir',  required=True,
                                        help='automatically download benchmark data into current directory')
    hg002_q100_eval_parser.add_argument('--happy-max-length', type=int, default=default_happy_max_length,
                                        help='ignore sites with alleles longer than this')
    hg002_q100_eval_parser.add_argument('--happy-threads', type=int, default=8,
                                        help='number of threads (default=8)')    
    hg002_q100_eval_parser.add_argument('--grch38-vcfs', nargs='*', default=[],
                                        help='vcfs for HG002-T2T-Q100-GRCh38 evaluation')
    hg002_q100_eval_parser.add_argument('--chm13-vcfs', nargs='*', default=[],
                                        help='vcfs for HG002-T2T-Q100-GRCh38 evaluation')
    hg002_q100_eval_parser.add_argument('--eval-type', default='both', choices = ['happy', 'truvari', 'both'],
                                        help='evaluation to run, can be happy, truvari or both')        
    

    args = parser.parse_args(command_line)

    if args.command in ['happy', 'truvari']:
        if not os.path.isdir(args.out_dir):
            os.makedirs(args.out_dir)

        assert args.truth.endswith('.vcf.gz')
        assert args.calls.endswith('.vcf.gz')
        assert args.truth != args.calls
        assert os.path.isfile(args.truth) and os.path.isfile(args.truth + '.tbi')
        assert os.path.isfile(args.calls) and os.path.isfile(args.calls + '.tbi')
        assert os.path.isfile(args.ref) and os.path.isfile(args.regions)

    if args.command == 'happy':
        assert not args.ref.endswith('.gz')
        # run some preprocessing (only on calls -- assume truth from giab is ready to go)
        hap_calls = os.path.join(args.out_dir, os.path.basename(args.calls.replace('.vcf.gz', '.hap.vcf.gz')))        
        vcf_preprocess(args.calls,
                       hap_calls,
                       args.exclude_y,
                       args.haploid_x,
                       args.sample,
                       args.max_length,
                       True)

        # run happy
        happy(args.truth, hap_calls, args.ref, args.regions, args.sample, args.out_dir,
              threads = args.threads,
              options=args.options,
              docker_image = args.docker)

    elif args.command == 'truvari':
        # run some preprocessing (only on calls -- assume truth from giab is ready to go)
        tv_calls = os.path.join(args.out_dir, os.path.basename(args.calls.replace('.vcf.gz', '.tv.vcf.gz')))        
        vcf_preprocess(args.calls,
                       tv_calls,
                       args.exclude_y,
                       args.haploid_x,
                       args.sample,
                       None,
                       True)

        # run truvari
        truvari(args.truth, tv_calls, args.ref, args.regions, args.sample, args.out_dir,
                options=args.options,
                docker_image = args.docker)

    elif args.command == 'happy-breakdown':
        table = happy_chromsplit(args.vcf)
        for row in table:
            print('\t'.join(row))

    elif args.command == 'truvari-breakdown':
        table = truvari_chromsplit(args.dir)
        for row in table:
            print('\t'.join(row))

    elif args.command == 'download-q100':
        download_q100()

    elif args.command == 'eval-q100':
        eval_q100(args.out_dir,
                  args.grch38_vcfs,
                  args.chm13_vcfs,
                  args.download,
                  args.happy_max_length,
                  args.happy_threads,
                  args.eval_type)

            
if __name__ == '__main__':
    main()
