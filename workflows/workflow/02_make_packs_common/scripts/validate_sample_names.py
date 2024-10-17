#!/usr/bin/env python3

####################################################################################################
# This script validates the sample names in the input file
# It compares them to the provided phenotypes csv file,
# Validates the format of the sample names and checks for duplicates.
# It also checks the format of the phenotype values.
#
# The script takes two input files:
# 1. samples_file: a csv file containing the following columns:
#   - sample: the sample name
#   - fq1: the path to the forward fastq file
#   - fq2: the path to the reverse fastq file, is optional
# 2. phenotypes_file: a csv file containing the following columns:
#   - Empty, sample or id: the sample/individual name
#   - Multiple columns with the phenotype name. It is expected to be numeric, or NA
#
# Usage: python validate_sample_names.py -s <samples_file> -p <phenotypes_file>
####################################################################################################

import os
import sys
import pandas as pd
import argparse

def validate_sample_names(samples_file, phenotypes_file):
    # Read the samples file
    samples = pd.read_csv(samples_file)
    # Read the phenotypes file
    phenotypes = pd.read_csv(phenotypes_file)

    # Check if the sample names are in the correct format
    # They should be alphanumeric and contain no special characters
    # except for underscores and hyphens
    for sample in samples['sample']:
        if not sample.replace('_', '').replace('-', '').isalnum():
            print(f"Error: The sample name {sample} is not in the correct format")
            sys.exit(1)
    print ('✓ Sample names are in the correct format')
    

    # Check if the sample names in the samples file are the same as in the phenotypes file
    samples_not_in_phenotypes = set(samples['sample']) - set(phenotypes['sample'])
    if samples_not_in_phenotypes:
        print(f"Error: The following samples are not in the phenotypes file: {samples_not_in_phenotypes}")
        sys.exit(1)
    print ('✓ Sample names in the phenotypes file are the same as in the samples file')


    # Check if the fastq files exist
    for fq in samples['fq1']:
        if not os.path.exists(fq):
            print(f"Error: The fastq file {fq} does not exist")
            sys.exit(1)
    # If fq2 is provided, check if it exists
    if 'fq2' in samples.columns:
        for fq in samples['fq2']:
            if not os.path.exists(fq):
                print(f"Error: The fastq file {fq} does not exist")
                sys.exit(1)
    print ('✓ Fastq files exist')

    
    # Check for duplicates in the samples file
    duplicates = samples[samples.duplicated('sample')]
    if not duplicates.empty:
        print(f"Error: The following samples are duplicated in the samples file: {duplicates['sample'].tolist()}")
        sys.exit(1)
    print ('✓ No duplicates in the samples file')


    # Check the format of the phenotype values
    for col in phenotypes.columns[1:]:
        if not phenotypes[col].apply(lambda x: str(x).replace('.', '').replace('E', '').replace('-', '').replace('NA', '').isnumeric()).all():
            print(f"Error: The phenotype values in column {col} are not numeric")
            sys.exit(1)
    print('✓ Phenotype values are numeric')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Validate the sample and phenotypes input files')
    parser.add_argument('-s', '--samples_file', help='The samples file')
    parser.add_argument('-p', '--phenotypes_file', help='The phenotypes file')
    args = parser.parse_args()

    validate_sample_names(args.samples_file, args.phenotypes_file)
