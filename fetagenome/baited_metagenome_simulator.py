#!/usr/bin/env python

import os
import glob
import argparse
from fetagenome import fetagenome


class GeneLocation:
    def __init__(self, fasta_file, contig_name, start_position, end_position):
        self.fasta_file = fasta_file
        self.contig_name = contig_name
        self.start_position = start_position
        self.end_position = end_position


def find_target_gene_locations(target_gene_fasta, genome_dir):
    genome_fastas = sorted(glob.glob(os.path.join(genome_dir, '*.fasta')))
    for genome_fasta in genome_fastas:
        # BLAST the target genes against the genome, finding location for each, store in GeneLocation object.
        pass
    # return list of GeneLocation objects?


def extract_gene_sequences(gene_locations, fragment_size, fragment_stdev, output_fasta):
    # For each gene, need to pull out the gene itself. Choose X random starting points within the gene,
    # and make the fragments we'd actually simulate from based on them. This way, should end up with most sequences
    # coming from inside the gene, and some from the tails around the gene, up to roughly the fragment size.
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulated a baited metagenome.')
    parser.add_argument('-t', '--targets',
                        required=True,
                        type=str,
                        help='Path to a FASTA-formatted file of targets that will be baited out in your metagenome.')
    parser.add_argument('-c', '--config_file',
                        required=True,
                        help='Path to a folder containing the genomes that your metagenome will be simulated from.')
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

    dependencies = ['blastn', 'makeblastdb', 'art_illumina']  # TODO: Check for these
    # TODO: Proportions? I'm somewhat unclear on exactly how the amount of each genome and amount of each gene within
    # each genome interacts.

    # For baiting out metagenome sequences, we need to do the following:
    # 1) Figure out where in each genome each target is, if it's there at all. Use BLAST for this.
    # 2) Once we know where each gene is within each genome, can then pull out fasta sequences of the gene + upstream
    # and downstream regions roughly equal to fragment size. Will likely need to pull out more than once, with the gene
    # sequences themselves being represented lots of times, and the upstream/downstream regions not as much.
    # 3) Once all of that stuff is pulled out, just need to simulate the FASTA file using ART (or whatever else).
    # 4) Also simulate some stuff against all the background genomes for the off_target_fraction.
