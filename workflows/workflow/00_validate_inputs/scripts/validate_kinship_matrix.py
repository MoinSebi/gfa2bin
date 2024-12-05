#!/usr/bin/env python3

####################################################################################################
# This scripts validates the kinship matrix file, by checking the following:
# 1. The file exists
# 2. The number of rows and columns is the same, and matches the number of samples in the samples file
# 3. The entries are numeric, either float or integer
#
# Usage: python validate_kinship_matrix.py -k <kinship_matrix> -s <samples_file>
####################################################################################################

import os
import sys
import pandas as pd
import argparse

def validate_kinship_matrix(kinship_matrix, samples_file):
    # Check if the kinship matrix file exists
    if not os.path.isfile(kinship_matrix):
        print(f"Error: The kinship matrix file {kinship_matrix} does not exist")
        sys.exit(1)

    # Check if the samples file exists
    if not os.path.isfile(samples_file):
        print(f"Error: The samples file {samples_file} does not exist")
        sys.exit(1)

    # Read the samples file
    samples_df = pd.read_csv(samples_file, sep="\t|;|,", dtype=str, engine='python')
    num_samples = samples_df.shape[0]

    # Read the kinship matrix file
    kinship_df = pd.read_csv(kinship_matrix, delim_whitespace=True, header=None)

    # Check if the number of rows and columns in the kinship matrix matches the number of samples
    if kinship_df.shape[0] != num_samples or kinship_df.shape[1] != num_samples:
        print(f"Error: The number of rows and columns in the kinship matrix file {kinship_matrix} does not match the number of samples in the samples file {samples_file}")
        print (f"Number of samples in the samples file: {num_samples}")
        print (f"Number of rows and columns in the kinship matrix: {kinship_df.shape[0]}, {kinship_df.shape[1]}")
        sys.exit(1)
    print ('✓ The number of rows and columns in the kinship matrix matches the number of samples in the samples file')

    # Check if the entries in the kinship matrix are numeric
    if not kinship_df.map(lambda x: isinstance(x, (int, float))).all().all():
        print(f"Error: The entries in the kinship matrix file {kinship_matrix} are not numeric")
        sys.exit(1)
    print ('✓ The entries in the kinship matrix are numeric')

    # Print a reminder that the order of the samples in the kinship matrix should match the order in the samples file
    print(f'! Important: Make sure the order of the samples in the kinship matrix should match the order in the samples file {samples_file}!')
    print('! The pipeline cannot check this, so please verify it manually.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Validate the kinship matrix file')
    parser.add_argument('-k', '--kinship_matrix', help='The kinship matrix file', required=True)
    parser.add_argument('-s', '--samples_file', help='The samples file', required=True)
    args = parser.parse_args()
    validate_kinship_matrix(args.kinship_matrix, args.samples_file)