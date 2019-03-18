#!/usr/bin/env python

import os
import csv
import shutil
import logging
import argparse
import datetime
import subprocess
import pkg_resources
import multiprocessing

import pysam
import tempfile
import numpy as np
from Bio import SeqIO


def main():
    version = get_version()
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
    parser.add_argument('-f', '--fetagenome_name',
                        type=str,
                        default='FetaGenome',
                        help='Name of your FetaGenome file. (so reads will be called FetaGenome_R1.fastq.gz and ' \
                             'FetaGenome_R2.fastq.gz)')
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
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=multiprocessing.cpu_count(),
                        help='Number of threads to run. Defaults to all cores on your machine.')
    parser.add_argument('--log',
                        type=str,
                        default='FetaGenome_log_{}.txt'.format(str(datetime.datetime.now()).replace(' ', '_').split('.')[0]),
                        help='Name for log file. Defaults to FetaGenome_log_YYYY-MM-DD_HH:MM:SS.txt')
    parser.add_argument('--coverage_report',
                        default=None,
                        type=str,
                        help='If you want a report on coverages for all of your contigs, specify a file name here '
                             'for them to be written to. This file will be in CSV format.')
    parser.add_argument('-v', '--version',
                        action='version',
                        version=version)
    parser.add_argument('-p', '--platform',
                        type=str,
                        default='HS25',
                        choices=['MSv1', 'HS25'],
                        help='Sequencing platform to simulate from. Choices are MSv1 (MiSeq) or '
                             'HS25 (HiSeq 2500). Defaults to HiSeq.')
    args = parser.parse_args()
    logging.info('Welcome to {}. Beginning metagenome simulation...'.format(version))
    dependencies_present = dependency_check()
    if dependencies_present is False:
        logging.error('ERROR: One or more dependencies not found. Please make sure all dependencies are installed and '
                      'available on your $PATH. Exiting...')
        quit(code=1)
    fetastrains = read_config_file(config_file=args.config_file,
                                   total_number_reads=args.number_reads,
                                   output_file=args.fetagenome_name,
                                   platform=args.platform,
                                   read_length=args.read_length,
                                   insert_size=args.insert_size,
                                   logfile=args.log,
                                   coverage_report=args.coverage_report,
                                   threads=args.threads)
    if os.path.isfile(args.fetagenome_name + '_R1.fastq.gz') or os.path.isfile(args.fetagenome_name + '_R2.fastq.gz'):
        logging.warning('WARNING: Your output file already exists. You\'ll probably want to quit and make a new file.')
    for fetastrain in fetastrains:
        logging.info('Simulating for {}'.format(fetastrain.assembly))
        fetastrain.calculate_contig_depths()
        fetastrain.find_reads_per_contig()
        fetastrain.simulate_reads()
    logging.info('Compressing output files...')
    cmd = 'gzip {} {}'.format(args.fetagenome_name + '_R1.fastq', args.fetagenome_name + '_R2.fastq')
    out, err = run_cmd(cmd)
    if args.log is not None:
        write_to_logfile(args.log, out, err, cmd)
    logging.info('Done!')


def get_version():
    try:
        version = 'FetaGenome {}'.format(pkg_resources.get_distribution('fetagenome').version)
    except pkg_resources.DistributionNotFound:
        version = 'FetaGenome (Unknown version)'
    return version


def dependency_check():
    dependencies_present = True
    dependencies = ['art_illumina', 'bbmap.sh', 'gzip', 'cat']
    for dep in dependencies:
        if shutil.which(dep) is None:
            logging.info('Checking for {}...\033[0;31;31m NOT FOUND'.format(dep))
            dependencies_present = False
        else:
            logging.info('Checking for {}...\033[0;32;32m Good'.format(dep))
    return dependencies_present


def write_to_logfile(logfile, out, err, cmd):
    """
    Writes stdout, stderr, and a command to a logfile
    :param logfile: Path to file to write output to.
    :param out: Stdout of program called, as a string
    :param err: Stderr of program called, as a string
    :param cmd: command that was used
    """
    with open(logfile, 'a+') as outfile:
        outfile.write('Command used: {}\n\n'.format(cmd))
        outfile.write('STDOUT: {}\n\n'.format(out))
        outfile.write('STDERR: {}\n\n'.format(err))


def read_config_file(config_file, total_number_reads, output_file, platform, read_length, insert_size,
                     logfile=None, coverage_report=None, threads=1):
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
                                    insert_size=insert_size,
                                    logfile=logfile,
                                    coverage_report=coverage_report,
                                    threads=threads)
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
    def __init__(self, forward_reads, reverse_reads, assembly, number_reads, metagenome_file, platform='HS25',
                 read_length=150, insert_size=250, logfile=None, coverage_report=None, threads=1):
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
        self.logfile = logfile
        self.coverage_report = coverage_report
        self.threads = threads

    def calculate_contig_depths(self):
        # First, create a bamfile using bbmap.
        with tempfile.TemporaryDirectory() as tmpdir:
            sorted_bam = os.path.join(tmpdir, 'sorted_bamfile.bam')
            cmd = 'bbmap.sh threads={threads} ref={assembly} in={forward_reads} in2={reverse_reads} out=stdout.bam nodisk | ' \
                  'samtools sort > {sorted_bam}'.format(forward_reads=self.forward_reads,
                                                        reverse_reads=self.reverse_reads,
                                                        sorted_bam=sorted_bam,
                                                        assembly=self.assembly,
                                                        threads=self.threads)
            logging.info('Aligning reads to assembly...')
            out, err = run_cmd(cmd)
            if self.logfile is not None:
                write_to_logfile(self.logfile, out, err, cmd)
            # Index bamfile or pysam can't figure out coverage.
            pysam.index(sorted_bam)
            logging.info('Calculating coverage...')
            bamfile = pysam.AlignmentFile(sorted_bam, 'rb')
            # Get coverage, and store in contig_depths dictionary
            for contig in SeqIO.parse(self.assembly, 'fasta'):
                a, c, g, t = bamfile.count_coverage(contig=contig.id)
                a, c, g, t = np.array(a), np.array(c), np.array(g), np.array(t)
                self.contig_depths[contig.id] = np.mean(a + c + g + t)
            bamfile.close()

    def find_reads_per_contig(self):
        logging.info('Finding number of reads to generate per contig...')
        sum_of_all = 0  # Bad variable name!
        proportions = dict()
        for contig in SeqIO.parse(self.assembly, 'fasta'):
            depth_length_product = len(contig.seq) * self.contig_depths[contig.id]
            sum_of_all += depth_length_product
            proportions[contig.id] = depth_length_product

        for contig in proportions:
            # Multiply by 0.5 since we'll simulate forward and reverse reads separately.
            self.reads_per_contig[contig] = round(0.5 * self.number_reads * (proportions[contig]/sum_of_all))

        # Also write out a report on coverage if desired.
        if self.coverage_report is not None:
            if not os.path.isfile(self.coverage_report):
                with open(self.coverage_report, 'w') as f:
                    f.write('Assembly,Contig,ContigDepth,AssemblyAverageDepth\n')

            # Calculate average depth across all contigs.
            average_assembly_depth = np.mean(list(self.contig_depths.values()))
            with open(self.coverage_report, 'a+') as f:
                for contig in self.contig_depths:
                    f.write('{assembly},{contig},{contigdepth},{averagedepth}\n'.format(assembly=self.assembly,
                                                                                        contig=contig,
                                                                                        contigdepth='%.2f' % self.contig_depths[contig],
                                                                                        averagedepth='%.2f' % average_assembly_depth))

    def simulate_reads(self):
        logging.info('Simulating reads!')
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
                if self.logfile is not None:
                    write_to_logfile(self.logfile, out, err, cmd)
                cmd = 'cat {} >> {}'.format(output_fastq + '1.fq', self.metagenome_file + '_R1.fastq')
                out, err = run_cmd(cmd)
                if self.logfile is not None:
                    write_to_logfile(self.logfile, out, err, cmd)
                cmd = 'cat {} >> {}'.format(output_fastq + '2.fq', self.metagenome_file + '_R2.fastq')
                out, err = run_cmd(cmd)
                if self.logfile is not None:
                    write_to_logfile(self.logfile, out, err, cmd)


if __name__ == '__main__':
    main()
