#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="fetagenome",
    version="0.1.6",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'FetaGenome = fetagenome.fetagenome:main',
            'generate_config_file = fetagenome.generate_config_file:main',
            'baited_metagenome_simulator = fetagenome.baited_metagenome_simulator:main'
        ],
    },
    author="Andrew Low",
    author_email='andrew.low@canada.ca',
    url='https://github.com/lowandrew/FetaGenome2',
    install_requires=['biopython',
                      'numpy']
)
