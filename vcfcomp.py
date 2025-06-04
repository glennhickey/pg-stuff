#!/usr/bin/env python3
"""
VCF comparison with hap.py and truvari, alongiside some postproecssing functions
"""

import os, sys
import subprocess
import pysam
import argparse
from collections import defaultdict

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
                 '/data/' + os.path.basename(calls_vcf),
                 '/data/' + os.path.basename(truth_vcf),
                 '-f', '/data/' + os.path.basename(bed_regions),
                 '-r', '/data/' + os.path.basename(ref_fasta),
                 '--threads', str(threads),
                 '-o', '/data/' + os.path.basename(output_dir)]
    if options:
        happy_cmd += options.split()

    docker_cmd = ['docker', 'run', '-it', '--rm', '-v', os.path.abspath(output_dir) + ':/data', docker_image, ' '.join(happy_cmd)]
    print(docker_cmd)
    subprocess.check_call(['docker', 'run', '-it', '--rm', '-v', os.path.abspath(output_dir) + ':/data', docker_image] + happy_cmd)

    # use Parsa's convention of having a happy.output.vcf.gz
    subprocess.check_call(['ln', '-sf', os.path.abspath(os.path.join(output_dir, output_dir + '.vcf.gz')),
                           os.path.abspath(os.path.join(output_dir, 'happy-output.vcf.gz'))])
    subprocess.check_call(['ln', '-sf', os.path.abspath(os.path.join(output_dir, output_dir + '.vcf.gz.tbi')),
                           os.path.abspath(os.path.join(output_dir, 'happy-output.vcf.gz.tbi'))])    

def happy_preprocess(input_vcf,
                     output_vcf,
                     remove_y = False,
                     haploid_x = False,
                     sample = None):
    """
    hap.py will fail if chrY is diploid (ie .|1).  So we need to flatten it to haploid.  Will also extract
    A sample if given, or just remove Y entirely if specified. 
    """
    tmp_vcf = input_vcf.replace('.vcf', '.c1.vcf'.format('.' + sample if sample else ''))
    view_cmd = ['bcftools', 'view', input_vcf, '-c', '1', '-Oz']
    if sample:
        view_cmd += ['-as', sample]
    with open(tmp_vcf, 'w') as tmp_vcf_file:
        subprocess.check_call(view_cmd, stdout=tmp_vcf_file)
    subprocess.check_call(['tabix', '-fp', 'vcf', tmp_vcf])

    tmp_vcf_file = pysam.VariantFile(tmp_vcf, 'rb')
    output_vcf_file = pysam.VariantFile(output_vcf, 'w', header=tmp_vcf_file.header)

    for var in tmp_vcf_file.fetch():
        skip = 'chry' in var.contig.lower() and remove_y
        flatten = not skip and ('chry' in var.contig.lower() or (haploid_x and 'chrx' in var.contig.lower()))
        if flatten:
            for sample in var.samples.values():
                gt = sample['GT']
                assert len(gt) in [1,2]
                if len(gt) == 2:
                    hap_gt = gt[0]
                    if hap_gt == '.' and gt[1] != '.':
                        hap_gt = gt[1]
                    sample['GT'] = tuple([hap_gt])
        if not skip:
            output_vcf_file.write(var)

    tmp_vcf_file.close()
    output_vcf_file.close()
    
    subprocess.check_call(['tabix', '-fp', 'vcf', output_vcf])
    os.remove(tmp_vcf)

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
        
                                         
def main(command_line=None):                     
    parser = argparse.ArgumentParser('VCF Comparison Stuff')
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
    happy_parser.add_argument('--outDir', required=True,
                       help='output directory')    
    happy_parser.add_argument('--sample',
                       help='subset to this sample')
    happy_parser.add_argument('--docker', default='jmcdani20/hap.py:v0.3.12',
                       help='use this docker image instead of the default')
    happy_parser.add_argument('--options', default='--gender=male --pass-only --engine=vcfeval',
                       help='use these options instead of the default (surround in quotes)')
    happy_parser.add_argument('--threads', type=int, default=8,
                       help='number of threads (default=8)')
    happy_parser.add_argument('--excludeY', action='store_true',
                       help='completely ignore chrY')

    happy_breakdown_parser = subparsers.add_parser('happy-breakdown', help='Make chromosome-decomposed table of hap.py results')
    happy_breakdown_parser.add_argument('--vcf', required=True,
                                        help='happy output VCF')
    

    args = parser.parse_args(command_line)
    if args.command == 'happy':
        if not os.path.isdir(args.outDir):
            os.makedirs(args.outDir)

        assert args.truth.endswith('.vcf.gz')
        assert args.calls.endswith('.vcf.gz')
        assert args.truth != args.calls
        assert os.path.isfile(args.truth) and os.path.isfile(args.truth + '.tbi')
        assert os.path.isfile(args.calls) and os.path.isfile(args.calls + '.tbi')
        assert os.path.isfile(args.ref) and os.path.isfile(args.regions)

        # run some preprocessing
        hap_truth = os.path.join(args.outDir, os.path.basename(args.truth.replace('.vcf.gz', '.hap.vcf.gz')))
        hap_calls = os.path.join(args.outDir, os.path.basename(args.calls.replace('.vcf.gz', '.hap.vcf.gz')))        
        happy_preprocess(args.truth, hap_truth, remove_y = args.excludeY, sample = args.sample)
        happy_preprocess(args.calls, hap_calls, remove_y = args.excludeY, sample = args.sample)

        # run happy
        happy(hap_truth, hap_calls, args.ref, args.regions, args.sample, args.outDir,
              threads = args.threads,
              options=args.options,
              docker_image = args.docker)

    elif args.command == 'happy-breakdown':
        table = happy_chromsplit(args.vcf)
        for row in table:
            print('\t'.join(row))
        
        
if __name__ == '__main__':
    main()
