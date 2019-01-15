#!/usr/bin/env python

import os
import shutil
import logging
import tempfile
import argparse
from io import StringIO
from fetagenome import fetagenome
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline


class GeneLocation:
    def __init__(self, fasta_file, contig_name, start_position, end_position):
        self.fasta_file = fasta_file
        self.contig_name = contig_name
        self.start_position = start_position
        self.end_position = end_position


def find_target_gene_locations(proportions_dictionary, target_fasta):
    # BLAST the target genes against the genome, finding location for each, store in GeneLocation object.
    gene_locations = list()
    # TODO: Get rid of os.systems and use something nicer.
    with tempfile.TemporaryDirectory() as tmpdir:
        for genome_fasta in proportions_dictionary:
            linked_fasta = os.path.join(tmpdir, os.path.split(genome_fasta)[1])
            os.symlink(genome_fasta, linked_fasta)
            cmd = 'makeblastdb -dbtype nucl -in {}'.format(linked_fasta)
            print(cmd)
            os.system(cmd)
            blastn = NcbiblastnCommandline(query=target_fasta, db=linked_fasta, outfmt=5)
            out, err = blastn()
            for record in NCBIXML.parse(StringIO(out)):
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        subject_length = len(hsp.sbjct)
                        # Allow hits that are at least 90 percent identical over 90 percent of length.
                        # TODO: Let the user parameters trickle through to here.
                        if hsp.align_length >= subject_length * 0.9 and hsp.positives >= subject_length * 0.9:
                            gene_loc = GeneLocation(fasta_file=genome_fasta,
                                                    )
    return gene_locations


def find_proportion_target_bases_each_genome(proportions_dictionary, gene_locations):
    # For each genome, figure out how much in terms of proportion is covered by target genes.
    # This information will get used to figure out how much of each gene to simulate.
    pass


def extract_gene_sequences(gene_locations, fragment_size, fragment_stdev, output_fasta):
    # For each gene, need to pull out the gene itself. Choose X random starting points within the gene,
    # and make the fragments we'd actually simulate from based on them. This way, should end up with most sequences
    # coming from inside the gene, and some from the tails around the gene, up to roughly the fragment size.
    pass


def dependency_check():
    dependencies = ['blastn', 'makeblastdb', 'art_illumina']
    all_deps_present = True
    for dependency in dependencies:
        if shutil.which(dependency) is None:
            logging.error('ERROR: Could not find dependency {}. Please make sure it is installed and on '
                          'your $PATH.'.format(dependency))
            all_deps_present = False
    if all_deps_present is False:
        quit()


def main():
    parser = argparse.ArgumentParser(description='Simulated a baited metagenome.')
    parser.add_argument('-t', '--targets',
                        required=True,
                        type=str,
                        help='Path to a FASTA-formatted file of targets that will be baited out in your metagenome.')
    parser.add_argument('-c', '--config_file',
                        required=True,
                        type=str,
                        help='Path to configuration file. Can be generated with generate_config_file from this package.')
    parser.add_argument('--off_target_fraction',
                        type=float,
                        default=0.2,
                        help='Fraction of reads that are off-target and not from your baited targets. Must be between '
                             '0 and 1, defaults to 0.2.')
    parser.add_argument('--fragment_size',
                        type=int,
                        default=400,
                        help='Average fragment size of the library being sequenced. Defaults to 400. This number '
                             'also affects how far past the ends of target sequences you\'ll see reads from.')
    parser.add_argument('--fragment_stdev',
                        type=int,
                        default=20,
                        help='Standard deviation in fragment size, defaults to 20.')
    parser.add_argument('--bait_percent_identity',
                        type=int,
                        default=90,
                        help='Percent identity to targets required for a gene to be baited. Defaults to 90 percent.')
    parser.add_argument('-n', '--number_reads',
                        default=10000000,
                        type=int,
                        help='Number of reads you want to generate for the metagenome. Defaults to 10 million.')
    args = parser.parse_args()

    # Setup the logger.
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')

    dependency_check()

    proportions_dictionary = fetagenome.parse_config_file(args.config_file)
    normalized_proportions = fetagenome.normalize_proportions(proportions_dictionary)

    # For baiting out metagenome sequences, we need to do the following:
    # 1) Figure out where in each genome each target is, if it's there at all. Use BLAST for this.
    gene_locations = find_target_gene_locations(proportions_dictionary=normalized_proportions,
                                                target_fasta=args.targets)
    # 2) Once we know where each gene is within each genome, can then pull out fasta sequences of the gene + upstream
    # and downstream regions roughly equal to fragment size. Will likely need to pull out more than once, with the gene
    # sequences themselves being represented lots of times, and the upstream/downstream regions not as much.
    # 3) Once all of that stuff is pulled out, just need to simulate the FASTA file using ART (or whatever else).
    # 4) Also simulate some stuff against all the background genomes for the off_target_fraction.


if __name__ == '__main__':
    main()
