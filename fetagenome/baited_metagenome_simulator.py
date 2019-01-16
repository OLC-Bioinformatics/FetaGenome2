#!/usr/bin/env python

# Core library
import os
import random
import shutil
import logging
import tempfile
import argparse
from io import StringIO

# Established packages
import numpy as np
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

# My own, project specific stuff
from fetagenome import fetagenome


# TODO: ALL THE DOCSTRINGS


class GeneLocation:
    def __init__(self, fasta_file, contig_name, start_position, end_position):
        self.fasta_file = fasta_file
        self.contig_name = contig_name
        self.start_position = start_position
        self.end_position = end_position
        # Give each gene location a name that should be unique
        self.name = '{}_{}_{}_{}'.format(self.fasta_file, self.contig_name, self.start_position, self.end_position)
        self.gene_length = abs(self.end_position - self.start_position)


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
                            contig_name = alignment.title.split()[1]
                            gene_loc = GeneLocation(fasta_file=genome_fasta,
                                                    contig_name=contig_name,
                                                    start_position=hsp.sbjct_start,
                                                    end_position=hsp.sbjct_end)
                            gene_locations.append(gene_loc)
    return gene_locations


def find_proportion_target_bases_each_gene(proportions_dict, gene_locations):
    # First, find how many total bases of target gene we have across all genomes.
    total_target_length = 0
    for gene_location in gene_locations:
        total_target_length += gene_location.gene_length
    # Now that we know how many bases of target gene total we have, go through and find the proportion for each target
    # gene by multiplying the proportion of that gene compared to all genes against the genome proportion. Normalize
    # after that.
    target_gene_proportion_dictionary = dict()
    for gene_location in gene_locations:
        target_gene_proportion_dictionary[gene_location.name] = proportions_dict[gene_location.fasta_file] * gene_location.gene_length/total_target_length

    normalized_target_gene_proportions = fetagenome.normalize_proportions(target_gene_proportion_dictionary)
    return normalized_target_gene_proportions


def create_fragment(gene_location, fragment_size, fragment_stdev):
    # Index the fasta file
    genome_index = SeqIO.index(gene_location.fasta_file, 'fasta')

    # Draw a fragment size from normal distribution.
    fragment_size = int(np.random.normal(fragment_size, fragment_stdev))

    # Figure out where in the gene we'll start from.
    if gene_location.end_position > gene_location.start_position:
        starting_base = random.randint(gene_location.start_position, gene_location.end_position)
    else:
        starting_base = random.randint(gene_location.end_position, gene_location.start_position)

    # Figure out whether to go in forward or reverse direction.
    if random.random() >= 0.5:  # Forward direction
        sequence = genome_index[gene_location.contig_name].seq[starting_base: starting_base + fragment_size]
    else:
        if starting_base - fragment_size < 0:
            start = 0
        else:
            start = starting_base - fragment_size
        sequence = genome_index[gene_location.contig_name].seq[start: starting_base]
    return sequence


def extract_gene_sequences(gene_locations, gene_proportions, fragment_size, fragment_stdev, output_fasta):
    total_genes_to_write = 10000
    genes_written = 1
    for gene_location in gene_locations:
        number_fastas_to_create = int(total_genes_to_write * gene_proportions[gene_location.name])
        for i in range(number_fastas_to_create):
            dna_fragment = create_fragment(gene_location=gene_location,
                                           fragment_size=fragment_size,
                                           fragment_stdev=fragment_stdev)
            with open(output_fasta, 'a+') as f:
                f.write('>seq{}\n'.format(genes_written))
                f.write('{}\n'.format(dna_fragment))
            genes_written += 1


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


def simulate_target_reads(extracted_targets_fasta, num_reads, fragment_size, fragment_stdev, output_dir):
    read_length = 150
    logging.info('Aiming for {} reads'.format(num_reads))
    # Find length of extracted_targets_fasta so we know what fold coverage we're shooting for
    fasta_length = fetagenome.find_genome_length(extracted_targets_fasta)
    logging.info('Total fasta length is {}'.format(fasta_length))

    # Get fold coverage desired by multiplying read length by num_reads to find how many bases we want, then
    # divide by how long the fasta actually is.
    # TODO: This is somehow getting only 3/4 of the number of reads I want for reasons I don't think I understand...
    fold_coverage = (num_reads * read_length)/fasta_length
    art_cmd = 'art_illumina -ss HS25 -i {} -l 150 -p -f {} -m {} -s {} -o {}'.format(extracted_targets_fasta,
                                                                                     fold_coverage,
                                                                                     fragment_size,
                                                                                     fragment_stdev,
                                                                                     os.path.join(output_dir, 'target'))
    logging.info(art_cmd)
    os.system(art_cmd)


def simulate_off_target_reads(normalized_proportions, number_reads, read_length, fragment_size, fragment_stdev, output_dir):
    fetagenome.simulate_reads(normalized_proportions=normalized_proportions,
                              desired_number_reads=number_reads,
                              read_length=read_length,
                              quality_shift=0,
                              platform='HS25',
                              insert_size=fragment_size,
                              output_dir=output_dir)


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
    parser.add_argument('-o', '--output_dir',
                        required=True,
                        type=str,
                        help='Output directory for reads. Will be created.')
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
                        help='Percent identity to targets required for a gene to be baited. Defaults to 90 percent. '
                             'WARNING: Does not actually do anything at this point.')
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

    # With the gene locations known, we can find how much of each gene should be represented in the baited
    # metagenome by weighting the gene length by genome proportion.
    gene_proportions = find_proportion_target_bases_each_gene(proportions_dict=normalized_proportions,
                                                              gene_locations=gene_locations)

    # Once that's done, extract genes from FASTA files, in proportions expected for each genome, with some end fragments
    # hanging around. Then, simulate 1 - off_target_fraction * number_reads reads from that file.
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    extracted_targets = os.path.join(args.output_dir, 'extracted_targets.fasta')
    extract_gene_sequences(gene_locations=gene_locations,
                           gene_proportions=gene_proportions,
                           fragment_size=args.fragment_size,
                           fragment_stdev=args.fragment_stdev,
                           output_fasta=extracted_targets)

    # Use the extracted targets fasta to simulate ART reads.
    simulate_target_reads(extracted_targets_fasta=extracted_targets,
                          num_reads=int((1 - args.off_target_fraction) * args.number_reads),
                          fragment_size=args.fragment_size,
                          fragment_stdev=args.fragment_stdev,
                          output_dir=args.output_dir)

    # Simulate off target reads by putting all genomes involved into a file and simulating off_target_fraction*number_reads
    simulate_off_target_reads(normalized_proportions=normalized_proportions,
                              number_reads=int(args.off_target_fraction * args.number_reads),
                              read_length=150,
                              fragment_size=args.fragment_size,
                              fragment_stdev=args.fragment_stdev,
                              output_dir=args.output_dir)


    # ???

    # Profit


if __name__ == '__main__':
    main()
