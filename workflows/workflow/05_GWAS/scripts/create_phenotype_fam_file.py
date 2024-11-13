#!/usr/bin/env python3
import csv
import sys

def main(fam_file_name, phenotype_file_name, output_fam_file_name):
    # Read the fam file
    with open(fam_file_name, 'r') as fam_file:
        reader = csv.reader(fam_file, delimiter='\t')
        fam_data = [list(row) for row in reader]

    # Read the phenotype file into a dictionary
    phenotype_data = {}
    with open(phenotype_file_name, 'r') as phenotype_file:
        reader = csv.DictReader(phenotype_file, delimiter='\t')
        for row in reader: 
            phenotype_data[row['sample']] = row['phenotype_value']
    # Replace the 6th column in the fam file based on the phenotype file
    for row in fam_data:
        #if row[0] in phenotype_data:
        #    # Check if value is numeric -- if so, write a float, otherwise write a string
        #    try:
        #        row[5] = float(phenotype_data[row[0]])
        #    except ValueError:
        row[5] = phenotype_data[row[0]]

    # Write the updated data back to a new fam file
    with open(output_fam_file_name, 'w', newline='') as updated_fam_file:
        writer = csv.writer(updated_fam_file, delimiter=' ')
        writer.writerows(fam_data)

    # # Print the updated fam file data for verification
    # for row in fam_data:
    #     print(row)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: create_phenotype_fam_file.py <fam_file> <phenotype_file> <output_fam_file>")
        sys.exit(1)
    
    fam_file_name = sys.argv[1]
    phenotype_file_name = sys.argv[2]
    output_fam_file_name = sys.argv[3]
    
    main(fam_file_name, phenotype_file_name, output_fam_file_name)
