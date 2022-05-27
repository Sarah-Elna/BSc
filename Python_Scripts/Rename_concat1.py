import csv
import os
import os.path

def read_csv(file_name, file_delimiter):
    id_list = []
    name_list = []
    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=file_delimiter)
        for row in csv_reader:
            id_list.append(row[1])
            name_list.append(row[2])
    return (id_list, name_list)


def rename_concat(id_list, name_list, concat1_path):
    new_file = open(concat1_path+'concat1_rename.fasta', 'a')
    original_file = open(concat1_path+'concat1.fasta', 'r')
    file_line_list = original_file.readlines()
    for line in file_line_list:
        if line[0] == '>':
            for i in range(0, len(id_list)):
                if id_list[i] == line[1:6]:
                    new_file.write('>'+name_list[i])
                    new_file.write('\n')
        if line[0] != '>':
            new_file.write(line)
    return ':-D'

## Run

rename = "/home/sarahe/GitHub/BSc/Renaming_csv_files/Rename_00_data_files.csv"
concat1_path = '/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/06_concatenate/'

id_list, name_list = read_csv(rename, ';')
print(id_list, name_list)

print(rename_concat(id_list, name_list, concat1_path))