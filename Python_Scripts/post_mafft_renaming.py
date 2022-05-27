import argparse, subprocess, csv
import csv
import os
import os.path

parser = argparse.ArgumentParser()
parser.add_argument("in_file")
args = parser.parse_args()
in_file = str(args.in_file)

def read_csv(file_name, file_delimiter):
    id_list = []
    name_list = []
    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=file_delimiter)
        for row in csv_reader:
            id_list.append(row[1])
            name_list.append(row[2])
    return (id_list, name_list)

def rename_mafft(id_list, name_list, concat1_path, in_file_path):
    file = open(in_file, 'r')
    file_line_list = file.readlines()
    for line in file_line_list:
        if line[0] == '>':
            sp_id = line[1:6]
            for i in range(0, len(id_list)):
                if sp_id == id_list[i]:
                    new_file.write('>'+name_list[i])
                    new_file.write('\n')
                else:
                    new_file.write(line)
    file.close()
    return ':-S'

rename = "/home/sarahe/GitHub/BSc/Renaming_csv_files/Rename_00_data_files.csv"
concat1_path = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/06_concatenate/"
in_file_path = in_file

id_list, name_list = read_csv(rename, ';')
#print(id_list, name_list)

rename_mafft(id_list, name_list, concat1_path, in_file_path)