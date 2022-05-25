import argparse, subprocess, csv

parser = argparse.ArgumentParser()
parser.add_argument("in_file")
args = parser.parse_args()
in_file = str(args.in_file)

def rename_targets(in_file_path):
    file = open(in_file, 'r')
    file_line_list = file.readlines()
    for line in file_line_list:
        if line != '' and line[0] == '>':
            new_line = line[0:6]
            new_line.strip()
            print(new_line)
        else:
            line.strip()
            print(line)
    file.close()
    return ':-S'

in_file_path = in_file

rename_targets(in_file_path)