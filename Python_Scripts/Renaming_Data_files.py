'''
#########################################################################################################################
This script is meant to rename files with a csv containing the old- and new names of the files.
To see the example see Renmae_00_data_files.csv.
Additional information in the csv file is species names, to be used later (maybe).
This script is meant to be run before the workflow.
When calling this script, three arguments are needed; csv complete path, path to not-renamed files, path to renamed files
#########################################################################################################################
'''

import argparse, subprocess, csv

parser = argparse.ArgumentParser()
parser.add_argument("csv_file")
parser.add_argument("path_in")
parser.add_argument("path_out")
args = parser.parse_args()
csv_file = str(args.csv_file)
path_in = str(args.path_in)
path_out = str(args.path_out)

file = open(csv_file)
csvreader = csv.reader(file)
header = next(csvreader)
rows = []
for row in csvreader:
    rows.append(row)
file.close()

endings = ["_clean-Read12-single.fastq", "_clean-Read1.fastq", "_clean-Read1-single.fastq", "_clean-Read2.fastq", "_clean-Read2-single.fastq"]

for i in range(0, len(rows)):
    for j in range(0, len(endings)):

        old_file_name = rows[i][0] + endings[j]
        new_file_name = rows[i][1] + endings[j]

        cmd = 'mv '+path_in+old_file_name+' '+path_out+new_file_name
        subprocess.call(cmd,shell=True)