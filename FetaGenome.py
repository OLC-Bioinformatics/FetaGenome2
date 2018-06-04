#!/usr/bin/env python

import os
import csv
import logging
import argparse
import subprocess
from Bio import SeqIO


def parse_config_file(config_file):
    proportions = dict()
    with open(config_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            strain = row['Strain']
            proportion = row['Proportion']
            proportions[strain] = float(proportion)
    return proportions


def normalize_proportions(proportions):
    # Will allow user to put whatever numbers they want for proportions (don't have to add to 1 or 100).
    # Do normalization so that everything works out to adding to one here.
    # First, get total number specified.
    total_count = 0.0
    for strain in proportions:
        total_count += proportions[strain]
    # Now modify all counts by dividing by total number so that sum of proportions adds to 1.
    for strain in proportions:
        proportions[strain] = proportions[strain]/total_count
    return proportions


def find_genome_length(fasta_file):
    genome_length = 0
    for contig in SeqIO.parse(fasta_file, 'fasta'):
        genome_length += len(contig)
    return genome_length


def create_fastq_from_fasta(fasta_file, output_fastq, depth=20, read_length=250, insert_size=350, insert_std=10,
                            platform='MSv1'):
    # Use ART to simulate us some reads
    cmd = 'art_illumina -ss {platform} -i {input_fasta} -l {read_length} -na -p -f {depth} -m {insert_size}' \
          ' -s {insert_std} -o {output_fastq}'.format(input_fasta=fasta_file,
                                                      output_fastq=output_fastq,
                                                      depth=depth,
                                                      read_length=read_length,
                                                      insert_size=insert_size,
                                                      insert_std=insert_std,
                                                      platform=platform)
    subprocess.call(cmd, shell=True)
    # Rename the files and gzip them for space savings.
    cmd = 'mv {fastq_name}1.fq {fastq_name}_R1.fastq'.format(fastq_name=output_fastq)
    os.system(cmd)
    cmd = 'mv {fastq_name}2.fq {fastq_name}_R2.fastq'.format(fastq_name=output_fastq)
    os.system(cmd)


def simulate_reads(normalized_proportions, desired_number_reads, output_dir):
    # With normalized proportions figured out, need to do some read simulation.
    # Use knowledge of genome size to figure out what depth needs to be simulated to so that we don't have to subsample
    # afterwards.
    for strain in normalized_proportions:
        genome_length = find_genome_length(strain)
        genome_num_reads = normalized_proportions[strain] * desired_number_reads
        # Genome coverage level: Get the number of bases we want by multiplying number of reads needed for this genome
        # by read length (250), then dividing by genome length to get coverage
        genome_coverage = (genome_num_reads * 250)/genome_length
        output_fastq_name = os.path.join(output_dir, os.path.split(strain)[-1].split('.')[0])
        create_fastq_from_fasta(fasta_file=strain,
                                output_fastq=output_fastq_name,
                                depth=genome_coverage)


def concatenate_into_fetagenome(output_dir, fetagenome_name):
    # Take advantage of the fact that cat is super speedy and make a system call to it. To create merged file.
    forward_fetagenome = os.path.join(output_dir, fetagenome_name + '_R1.fastq')
    reverse_fetagenome = os.path.join(output_dir, fetagenome_name + '_R2.fastq')
    cmd = 'cat {forward_fastqs} > {merged_forward}'.format(forward_fastqs=os.path.join(output_dir, '*_R1.fastq'),
                                                           merged_forward=forward_fetagenome)
    os.system(cmd)
    cmd = 'cat {reverse_fastqs} > {merged_reverse}'.format(reverse_fastqs=os.path.join(output_dir, '*_R2.fastq'),
                                                           merged_reverse=reverse_fetagenome)
    os.system(cmd)
    # Now gzip the forward and reverse files, and then get rid of all the other files within the folder.
    cmd = 'gzip {forward_feta} {reverse_feta}'.format(forward_feta=forward_fetagenome,
                                                      reverse_feta=reverse_fetagenome)
    os.system(cmd)
    cmd = 'rm {fastq_files}'.format(fastq_files=os.path.join(output_dir, '*.fastq'))
    os.system(cmd)


def main():
    parser = argparse.ArgumentParser(description='Given a configuration file, will create a FetaGenome from FASTA files'
                                                 ' by simulating reads with ART and pasting reads together into a '
                                                 'FetaGenome.')
    parser.add_argument('-c', '--config_file',
                        type=str,
                        required=True,
                        help='Path to your configuration file for FetaGenome creation.')
    parser.add_argument('-o', '--output_dir',
                        type=str,
                        required=True,
                        help='Output directory for your FetaGenome files. Will be created if it does not'
                             ' already exist.')
    parser.add_argument('-n', '--number_reads',
                        type=int,
                        default=1000000,
                        help='Number of reads to include in FetaGenome. Defaults to 1000000.')
    parser.add_argument('-f', '--fetagenome_name',
                        type=str,
                        default='FetaGenome',
                        help='Name of your FetaGenome file. Defaults to FetaGenome (so reads will be called '
                             'FetaGenome_R1.fastq.gz and FetaGenome_R2.fastq.gz)')
    args = parser.parse_args()
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    proportions = parse_config_file(args.config_file)
    normalized_proportions = normalize_proportions(proportions)
    simulate_reads(normalized_proportions=normalized_proportions,
                   desired_number_reads=args.number_reads,
                   output_dir=args.output_dir)
    concatenate_into_fetagenome(output_dir=args.output_dir,
                                fetagenome_name=args.fetagenome_name)


if __name__ == '__main__':
    main()
