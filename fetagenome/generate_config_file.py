#!/usr/bin/env python

import os
import glob
import argparse


def main():
    parser = argparse.ArgumentParser(description='Helper script for generating FetaGenome config files. Given a '
                                                 'directory with FASTA files in it, will create a config file where '
                                                 'each FASTA file in that directory is part of the FetaGenome, all in '
                                                 'equal proportions.')
    parser.add_argument('-i', '--input_dir',
                        type=str,
                        required=True,
                        help='Directory containing the FASTA files you want to use to create your FetaGenome.')
    parser.add_argument('--fastq',
                        default=False,
                        action='store_true',
                        help='Activate this flag to create a FASTQ-formatted config file. This assumes that your '
                             'read files are paired and contain _R1/_R2 in the filenames.')
    parser.add_argument('-c', '--config_file',
                        type=str,
                        default='FetaGenomeConfig.csv',
                        help='Name of FetaGenome config file you want to create. Defaults to FetaGenomeConfig.csv')
    args = parser.parse_args()
    if args.fastq is False:
        with open(args.config_file, 'w') as outfile:
            outfile.write('Strain,Proportion\n')
            fasta_files = glob.glob(os.path.join(args.input_dir, '*.f*a'))
            for fasta in fasta_files:
                outfile.write('{},1\n'.format(os.path.abspath(fasta)))
    elif args.fastq is True:
        with open(args.config_file, 'w') as outfile:
            outfile.write('ForwardReads,ReverseReads,Proportion\n')
            fastq_files = glob.glob(os.path.join(args.input_dir, '*_R1_*.f*q*'))
            for fastq_file in fastq_files:
                reverse_reads = fastq_file.replace('_R1', '_R2')
                if os.path.isfile(reverse_reads):
                    outfile.write('{forward_reads},{reverse_reads},1'.format(forward_reads=os.path.abspath(fastq_file),
                                                                             reverse_reads=os.path.abspath(reverse_reads)))


if __name__ == '__main__':
    main()
