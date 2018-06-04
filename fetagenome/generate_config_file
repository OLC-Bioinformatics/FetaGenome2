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
    parser.add_argument('-c', '--config_file',
                        type=str,
                        default='FetaGenomeConfig.csv',
                        help='Name of FetaGenome config file you want to create. Defaults to FetaGenomeConfig.csv')
    args = parser.parse_args()
    with open(args.config_file, 'w') as outfile:
        outfile.write('Strain,Proportion\n')
        fasta_files = glob.glob(os.path.join(args.input_dir, '*.f*a'))
        for fasta in fasta_files:
            outfile.write('{},1\n'.format(os.path.abspath(fasta)))


if __name__ == '__main__':
    main()
