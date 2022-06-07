import argparse, subprocess, csv

parser = argparse.ArgumentParser()
parser.add_argument("in_file")
args = parser.parse_args()
in_file = str(args.in_file)

duplicats_startswith_0 = ['DHET0', 'DHET1', 'MRDA0', 'MRDA1', 'MRDA2', 'MRDA3', 'MSKO0', 'MSKO1', 'DLEU0', 'DLEU1', 'MRIN0', 'MRIN1', 'MRIN2']
duplicats_other = ['DLAN1', 'DLAN2', 'DLAN3']

def rename_targets(in_file_path):
    file = open(in_file, 'r')
    file_line_list = file.readlines()
    for line in file_line_list:
        if line != '' and line[0] == '>' and line[1] == '_'
            new_line = line[0] + line[4:9]
            for i in range(0, len(duplicats_startswith_0)):
                if duplicats_startswith_0[i] in new_line:
                    new_line = line[0:5]+'0'
            for j in range(0, len(duplicats_other)):
                if duplicats_other[j] in new_line:
                    new_line = line[0:5]+'1'
            new_line.strip()
            print(new_line, end='\n')
        if line != '' and line[0] == '>':
            new_line = line[0:6]
            for i in range(0, len(duplicats_startswith_0)):
                if duplicats_startswith_0[i] in new_line:
                    new_line = line[0:5]+'0'
            for j in range(0, len(duplicats_other)):
                if duplicats_other[j] in new_line:
                    new_line = line[0:5]+'1'
            new_line.strip()
            print(new_line, end='\n')
        else:
            line.strip()
            print(line, end='')
    file.close()
    return ':-S'

in_file_path = in_file

rename_targets(in_file_path)