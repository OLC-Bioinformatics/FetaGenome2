# FetaGenome Generator

Scripts in this repository allow for the creation of 'Fake Metagenomes' (aka FetaGenomes) with known community
compositions.

In order to have these scripts work, you'll need to following installed and available on your `$PATH`:

- [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm)
- [Python3](https://www.python.org/downloads/)
- Some sort of Unix-based system with basic utilties (gzip, cat, rm) installed.

### Installing FetaGenome Generator

Ideally, create a virtualenv (some instructions can be found [here](https://packaging.python.org/guides/installing-using-pip-and-virtualenv/)),
and then install FetaGenome via pip: `pip install fetagenome`

### Making FetaGenomes

To make FetaGenomes, run `FetaGenome` - this will create a mock fetagenome with known proportions as defined
in a config file that is given as one of the arguments to the script.

The config file should have two columns, comma-separated - the first should be `Strain`, and the second should be `Proportion`.
Put the absolute path to FASTA files you want to use in the first column and the relative proportion in the second.
If you want to have this process simplified, put all the FASTA files you're interested in having as part of your
FetaGenome into a folder and run `generate_contig_file -i /path/to/folder`. This will create a file called
`FetaGenomeConfig.csv` in your current working directory with all FASTA-formatted files in the folder specified as
strains in equal proportions.

```
usage: FetaGenome [-h] -c CONFIG_FILE -o OUTPUT_DIR [-n NUMBER_READS]
                  [-f FETAGENOME_NAME] [-q QUALITY_SHIFT]

Given a configuration file, will create a FetaGenome from FASTA files by
simulating reads with ART and pasting reads together into a FetaGenome.

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG_FILE, --config_file CONFIG_FILE
                        Path to your configuration file for FetaGenome
                        creation.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory for your FetaGenome files. Will be
                        created if it does not already exist.
  -n NUMBER_READS, --number_reads NUMBER_READS
                        Number of reads to include in FetaGenome. Defaults to
                        1000000.
  -f FETAGENOME_NAME, --fetagenome_name FETAGENOME_NAME
                        Name of your FetaGenome file. Defaults to FetaGenome
                        (so reads will be called FetaGenome_R1.fastq.gz and
                        FetaGenome_R2.fastq.gz)
  -q QUALITY_SHIFT, --quality_shift QUALITY_SHIFT
                        By default, ART will simulate Illumina reads with
                        fairly high quality. If you want this changed, you can
                        make them even higher quality with a positive integer
                        (to shift up by 2 on average, enter 2) or make them
                        lower quality with a negative number.
```
