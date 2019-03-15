#!/usr/bin/env python

import os
import csv
import logging
import argparse
import subprocess

import pysam
import tempfile
import numpy as np
from Bio import SeqIO


def main():
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    parser = argparse.ArgumentParser(description='Given a configuration file, will create a FetaGenome from FASTA files'
                                                 ' by simulating reads with ART and pasting reads together into a '
                                                 'FetaGenome.')
    parser.add_argument('-c', '--config_file',
                        type=str,
                        required=True,
                        help='Path to your configuration file for FetaGenome creation.')
    parser.add_argument('-o', '--output_file',
                        type=str,
                        required=True,
                        help='Output file for your FetaGenome.')
    parser.add_argument('-n', '--number_reads',
                        type=int,
                        default=1000000,
                        help='Number of reads to include in FetaGenome. Defaults to 1000000.')
    parser.add_argument('-l', '--read_length',
                        type=int,
                        default=150,
                        help='Read length. Defaults to 150.')
    parser.add_argument('-i', '--insert_size',
                        type=int,
                        default=250,
                        help='Insert size. Defaults to 250.')
    parser.add_argument('-p', '--platform',
                        type=str,
                        default='HS25',
                        choices=['MSv1', 'HS25'],
                        help='Sequencing platform to simulate from. Choices are MSv1 (MiSeq) or '
                             'HS25 (HiSeq 2500). Defaults to HiSeq.')
    args = parser.parse_args()
    fetastrains = read_config_file(config_file=args.config_file,
                                   total_number_reads=args.number_reads,
                                   output_file=args.output_file,
                                   platform=args.platform,
                                   read_length=args.read_length,
                                   insert_size=args.insert_size)
    for fetastrain in fetastrains:
        logging.info('Simulating for {}'.format(fetastrain.assembly))
        fetastrain.calculate_contig_depths()
        fetastrain.find_reads_per_contig()
        fetastrain.simulate_reads()
    logging.info('Done!')


def read_config_file(config_file, total_number_reads, output_file, platform, read_length, insert_size):
    proportions = dict()
    fetastrains = list()
    # First read through - figure out what total is so we can proportion properly
    with open(config_file) as csvfile:
        total = 0
        reader = csv.DictReader(csvfile)
        for row in reader:
            assembly = row['Strain']
            proportion = float(row['Proportion'])
            proportions[assembly] = proportion
            total += proportion
    # Second read through - create a FetaStrain object for each of our strains now that we know how much of each genome
    # should be present.
    with open(config_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            assembly = row['Strain']
            forward_reads = row['ForwardReads']
            reverse_reads = row['ReverseReads']
            fetastrain = FetaStrain(forward_reads=forward_reads,
                                    reverse_reads=reverse_reads,
                                    assembly=assembly,
                                    number_reads=total_number_reads * (proportions[assembly]/total),
                                    metagenome_file=output_file,
                                    platform=platform,
                                    read_length=read_length,
                                    insert_size=insert_size)
            fetastrains.append(fetastrain)
    return fetastrains


def run_cmd(cmd):
    """
    Runs a command using subprocess, and returns both the stdout and stderr from that command
    If exit code from command is non-zero, raises subproess.CalledProcessError
    :param cmd: command to run as a string, as it would be called on the command line
    :return: out, err: Strings that are the stdout and stderr from the command called.
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, cmd=cmd)
    return out, err


class FetaStrain:
    def __init__(self, forward_reads, reverse_reads, assembly, number_reads, metagenome_file, platform='HS25', read_length=150, insert_size=250):
        self.forward_reads = forward_reads
        self.reverse_reads = reverse_reads
        self.assembly = assembly
        self.number_reads = number_reads
        self.contig_depths = dict()
        self.reads_per_contig = dict()
        self.platform = platform
        self.read_length = read_length
        self.insert_size = insert_size
        self.metagenome_file = metagenome_file

    def calculate_contig_depths(self):
        # First, create a bamfile using bbmap.
        with tempfile.TemporaryDirectory() as tmpdir:
            sorted_bam = os.path.join(tmpdir, 'sorted_bamfile.bam')
            cmd = 'bbmap.sh ref={assembly} in={forward_reads} in2={reverse_reads} out=stdout.bam nodisk | ' \
                  'samtools sort > {sorted_bam}'.format(forward_reads=self.forward_reads,
                                                        reverse_reads=self.reverse_reads,
                                                        sorted_bam=sorted_bam,
                                                        assembly=self.assembly)
            logging.info('Aligning reads to assembly...')
            out, err = run_cmd(cmd)
            pysam.index(sorted_bam)
            logging.info('Calculating coverage...')
            bamfile = pysam.AlignmentFile(sorted_bam, 'rb')
            for contig in SeqIO.parse(self.assembly, 'fasta'):
                a, c, g, t = bamfile.count_coverage(contig=contig.id)
                a, c, g, t = np.array(a), np.array(c), np.array(g), np.array(t)
                self.contig_depths[contig.id] = np.mean(a + c + g + t)
            bamfile.close()

    def find_reads_per_contig(self):
        sum_of_all = 0  # Bad variable name!
        proportions = dict()
        for contig in SeqIO.parse(self.assembly, 'fasta'):
            depth_length_product = len(contig.seq) * self.contig_depths[contig.id]
            sum_of_all += depth_length_product
            proportions[contig.id] = depth_length_product

        for contig in proportions:
            # Multiply by 0.5 since we'll simulate forward and reverse reads separately.
            self.reads_per_contig[contig] = round(0.5 * self.number_reads * (proportions[contig]/sum_of_all))

    def simulate_reads(self):
        if os.path.isfile(self.metagenome_file + '_R1.fastq'):
            logging.warning('WARNING: Output FASTQ file already exists!')
        with tempfile.TemporaryDirectory() as tmpdir:
            for contig in SeqIO.parse(self.assembly, 'fasta'):
                temp_fasta = os.path.join(tmpdir, 'tmp_fasta.fasta')
                SeqIO.write([contig], temp_fasta, 'fasta')
                output_fastq = os.path.join(tmpdir, contig.id)
                cmd = 'art_illumina -ss {platform} -i {temp_fasta} -l {read_length} -na -p -c {num_reads} ' \
                      '-m {insert_size} -s 10 -o {output_fastq}'.format(platform=self.platform,
                                                                        temp_fasta=temp_fasta,
                                                                        read_length=self.read_length,
                                                                        num_reads=self.reads_per_contig[contig.id],
                                                                        insert_size=self.insert_size,
                                                                        output_fastq=output_fastq)
                out, err = run_cmd(cmd)
                os.system('cat {} >> {}'.format(output_fastq + '1.fq', self.metagenome_file + '_R1.fastq'))
                os.system('cat {} >> {}'.format(output_fastq + '2.fq', self.metagenome_file + '_R2.fastq'))


if __name__ == '__main__':
    main()
