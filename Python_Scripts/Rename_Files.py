## This code has been written to rename files.
# These files have four distinct cifers at the begining of their name before they are renamed.
# These four cifers correspond to a scientific name in a csv file called 'rename.csv'

## Reading in the cvs file
import csv
import os
import os.path

def read_csv(file_name, file_delimiter):
    number_list = []
    name_list = []
    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=file_delimiter)
        for row in csv_reader:
            number_list.append(row[0])
            name_list.append(row[1])
    return (number_list, name_list)

def rename_files(number_list, name_list, ending_list, path_in):
    for i in range(len(number_list)):
        for j in range(len(ending_list)):
            old_ending = ending_list[j].replace('Read', 'READ')
            path = ("{}{}{}.fastq".format(path_in, number_list[i], old_ending))
            isFile = os.path.isfile(path)
            if (isFile == True):
                print(True)
                old_name = r"{}/{}{}.fastq".format(path_in, number_list[i], old_ending)
                new_name = r"{}/{}{}.fastq".format(path_in, name_list[i], ending_list[j])
                os.rename(old_name, new_name)
            else:
                continue
    return ':-D'

## Run

rename = "/home/sarahe/GitHub/BSc/Renaming_csv_files/Rename_Files.csv"
ending_list = ['_clean-Read1', '_clean-Read1-single', '_clean-Read2', '_clean-Read2-single', '_clean-Read12-single']
path_in = '/home/sarahe/BSc/00_data/'

number_list, name_list = read_csv(rename, ';')
print(number_list, name_list)

print(rename_files(number_list, name_list, ending_list, path_in))
