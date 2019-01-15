#!/usr/bin/env python

import os
import csv
import shutil
import logging
import argparse
import subprocess
from Bio import SeqIO
import multiprocessing


def parse_config_file(config_file):
    """
    Looks through a specified config file (which is a csv with 2 columns) and creates a dictionary with strain and
    proportion information.
    :param config_file: Path to CSV-formatted configuration file.
    :return: Dictionary with path to strain file (uncompresed fasta/multifasta) as key and proportion (as float)
    as the value
    """
    # TODO: Implement some checks that config file is properly formatted/will work and give nice error messages if not
    proportions = dict()
    with open(config_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            strain = row['Strain']
            proportion = row['Proportion']
            proportions[strain] = float(proportion)
    return proportions


def normalize_proportions(proportions):
    """
    Proportions specified in config file are all just relative - this normalizes them to a proportion out of 1.0
    :param proportions: The proportions dictionary generated by parse_config_file
    :return: A dictionary in the same format as that generated by parse_config_file, but with normalized values.
    """
    # First, get total number specified.
    total_count = 0.0
    for strain in proportions:
        total_count += proportions[strain]
    # Now modify all counts by dividing by total number so that sum of proportions adds to 1.
    for strain in proportions:
        proportions[strain] = proportions[strain]/total_count
    return proportions


def find_genome_length(fasta_file):
    """
    Given a FASTA file, finds the total length of all contigs within that file.
    :param fasta_file: Path to a fasta-formatted, uncompressed file.
    :return: Total length contigs in that fasta file, as an integer.
    """
    genome_length = 0
    for contig in SeqIO.parse(fasta_file, 'fasta'):
        genome_length += len(contig)
    return genome_length


def create_fastq_from_fasta(fasta_file, output_fastq, depth=20, read_length=250, insert_size=350, insert_std=10,
                            platform='MSv1', quality_shift=0):
    """
    Creates a command to send to art_illumina to create simulated FASTQ files. Also does some renaming,
    since I don't like ART's naming scheme for paired-end reads.
    :param fasta_file: Path to uncompressed fasta file to simulate reads from.
    :param output_fastq: Base name for output fastq file.
    :param depth: Desired coverage depth.
    :param read_length: Read length - allowable values for this change based on the platform chosen.
    :param insert_size: Desired insert size (by insert, ART means fragment length)
    :param insert_std: Standard deviation of insert size.
    :param platform: Sequencing platform. Defaults to MSv1, since MSv3 quality profiles look very wrong.
    :param quality_shift: Amount to quality shift bases by. Positive value for better quality, negative for lower.
    :return: cmd - the command to run which will run ART, and rename files to make
    """
    # Use ART to simulate us some reads
    cmd = 'art_illumina -ss {platform} -i {input_fasta} -l {read_length} -na -p -f {depth} -m {insert_size}' \
          ' -s {insert_std} -o {output_fastq} -qs {quality_shift} -qs2 {quality_shift}'.format(input_fasta=fasta_file,
                                                                                               output_fastq=output_fastq,
                                                                                               depth=depth,
                                                                                               read_length=read_length,
                                                                                               insert_size=insert_size,
                                                                                               insert_std=insert_std,
                                                                                               platform=platform,
                                                                                               quality_shift=quality_shift)
    # Rename the files to use the _R1/_R2 convention.
    cmd += ' && mv {fastq_name}1.fq {fastq_name}_R1.fastq'.format(fastq_name=output_fastq)
    cmd += ' && mv {fastq_name}2.fq {fastq_name}_R2.fastq'.format(fastq_name=output_fastq)
    return cmd


def run_cmd(cmd, log_message):
    logging.info(log_message)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    messages = cmd + '\n' + stdout.decode() + '\n' + stderr.decode() + '\n'
    return messages


def simulate_reads(normalized_proportions, desired_number_reads, output_dir, quality_shift, read_length, platform, threads=1,
                   logfile='fetagenome.log', insert_size=200):
    """
    Simulates reads in proportions specified.
    :param normalized_proportions: Normalized proportion dictionary created by normalize_proportions
    :param desired_number_reads: Number of total reads in MetaGenome
    :param output_dir: Directory to store created fastq files.
    :param quality_shift: Amount to quality shift simulated reads by.
    """
    # Use knowledge of genome size to figure out what depth needs to be simulated to so that we don't have to subsample
    # afterwards.
    commands = list()
    log_messages = list()
    for strain in normalized_proportions:
        log_messages.append('Simulating reads for {}'.format(strain))
        genome_length = find_genome_length(strain)
        genome_num_reads = normalized_proportions[strain] * desired_number_reads
        # Genome coverage level: Get the number of bases we want by multiplying number of reads needed for this genome
        # by read length (250), then dividing by genome length to get coverage
        genome_coverage = (genome_num_reads * read_length)/genome_length
        output_fastq_name = os.path.join(output_dir, os.path.split(strain)[-1].split('.')[0])
        commands.append(create_fastq_from_fasta(fasta_file=strain,
                                                output_fastq=output_fastq_name,
                                                depth=genome_coverage,
                                                quality_shift=quality_shift,
                                                read_length=read_length,
                                                platform=platform,
                                                insert_size=insert_size))
    pool = multiprocessing.Pool(processes=threads)
    messages = pool.starmap(run_cmd, zip(commands, log_messages))
    pool.close()
    pool.join()
    with open(logfile, 'a+') as f:
        for message in messages:
            f.write(message)


def concatenate_into_fetagenome(output_dir, fetagenome_name, logfile='fetagenome.log'):
    """
    Takes individual fastq files created by simulate_reads and merges them into a FetaGenome, then cleans up
    individual files and gzips the FetaGenome for space savings.
    :param output_dir: Directory where simulate_reads put individual fastq files, and the merged fetagenome will go.
    :param fetagenome_name: Base name for FetaGenome. _R1/_R2 suffixes/file extensions will be added to this.
    :param logfile: Path to logfile to store commands used.
    """
    # Take advantage of the fact that cat is super speedy and make a system call to it. To create merged file.
    forward_fetagenome = os.path.join(output_dir, fetagenome_name + '_R1.fastq')
    reverse_fetagenome = os.path.join(output_dir, fetagenome_name + '_R2.fastq')
    with open(logfile, 'a+') as f:
        cmd = 'cat {forward_fastqs} > {merged_forward}'.format(forward_fastqs=os.path.join(output_dir, '*_R1.fastq'),
                                                               merged_forward=forward_fetagenome)
        os.system(cmd)
        f.write(cmd + '\n')
        cmd = 'cat {reverse_fastqs} > {merged_reverse}'.format(reverse_fastqs=os.path.join(output_dir, '*_R2.fastq'),
                                                               merged_reverse=reverse_fetagenome)
        os.system(cmd)
        f.write(cmd + '\n')
        # Now gzip the forward and reverse files, and then get rid of all the other files within the folder.
        cmd = 'gzip {forward_feta} {reverse_feta}'.format(forward_feta=forward_fetagenome,
                                                          reverse_feta=reverse_fetagenome)
        os.system(cmd)
        f.write(cmd + '\n')
        cmd = 'rm {fastq_files}'.format(fastq_files=os.path.join(output_dir, '*.fastq'))
        os.system(cmd)
        f.write(cmd + '\n')


def main():
    num_cpu_cores = multiprocessing.cpu_count()
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
    parser.add_argument('-q', '--quality_shift',
                        type=int,
                        default=0,
                        help='By default, ART will simulate Illumina reads with fairly high quality. If you want '
                             'this changed, you can make them even higher quality with a positive integer (to shift up'
                             ' by 2 on average, enter 2) or make them lower quality with a negative number.')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=num_cpu_cores,
                        help='Number of threads to run, allows for much faster simulation of reads. Defaults to '
                             'number of cores on your machine.')
    parser.add_argument('-l', '--read_length',
                        type=int,
                        default=250,
                        help='Read length. Defaults to 250.')
    parser.add_argument('-i', '--insert_size',
                        type=int,
                        default=300,
                        help='Insert size. Defaults to 300.')
    parser.add_argument('-p', '--platform',
                        type=str,
                        default='MSv1',
                        choices=['MSv1', 'HS25'],
                        help='Sequencing platform to simulate from. Choices are MSv1 (MiSeq) or '
                             'HS25 (HiSeq 2500). Defaults to MiSeq')
    args = parser.parse_args()

    # Setup the logger.
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')

    # Check that ART is available on the $PATH, and boot user if it isn't
    if shutil.which('art_illumina') is None:
        logging.error('ERROR: Could not find \'art_illumina\' executable on $PATH. Please install ART and'
                      ' then try running FetaGenome')

    # Make our fetagenomes.
    if not os.path.isdir(args.output_dir):
        logging.info('Created output directory {}'.format(args.output_dir))
        os.makedirs(args.output_dir)
    logging.info('Parsing configuration file.')
    proportions = parse_config_file(args.config_file)
    normalized_proportions = normalize_proportions(proportions)
    logging.info('Proportions put in were:')
    for proportion in normalized_proportions:
        logging.info('{}:{}'.format(proportion, normalized_proportions[proportion]))
    simulate_reads(normalized_proportions=normalized_proportions,
                   desired_number_reads=args.number_reads,
                   output_dir=args.output_dir,
                   quality_shift=args.quality_shift,
                   threads=args.threads,
                   logfile=os.path.join(args.output_dir, 'fetagenome.log'),
                   read_length=args.read_length,
                   platform=args.platform,
                   insert_size=args.insert_size)
    logging.info('Merging created reads into FetaGenome and cleaning up.')
    concatenate_into_fetagenome(output_dir=args.output_dir,
                                fetagenome_name=args.fetagenome_name,
                                logfile=os.path.join(args.output_dir, 'fetagenome.log'))
    logging.info('FetaGenome run complete.')


if __name__ == '__main__':
    main()
