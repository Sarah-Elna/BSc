## This script has been written with the purpouse of taking a Fasta file with multiple sequences and renaming them so the file can be used as a target file for HybPiper.

import os
import csv

gene_name_list = ['accD', 'atpA', 'atpB', 'atpF', 'atpH', 'atpI', 'ccsA', 'clpP', 'matK', 'ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'petA', 'petB', 'petD', 'petL', 'petN', 'psaA', 'psaB', 'psaC', 'psaI', 'psbA', 'psbB', 'psbC', 'psbD', 'psbK', 'psbL', 'psbZ', 'psI', 'rbcL', 'rpl16', 'rpl2', 'rpl22', 'rpl23', 'rpl32', 'rpoA', 'rpoB', 'rpoC1', 'rpoC2', 'rps11', 'rps12', 'rps15', 'rps16', 'rps18', 'rps2', 'rps3', 'rps4', 'rps7', 'rps8', 'rrn16', 'rrn23', 'rrn4,5', 'rrn5', 'rrn5 Dio', 'trnA', 'trnC', 'trnD', 'trnfM', 'trnG', 'trnH', 'trnI', 'trnK', 'trnL', 'trnN', 'trnP', 'trnQ', 'trnS', 'trnT', 'trnV', 'ycf1', 'ycf2', 'ycf3', 'ycf4', 'ycl2']

def get_headers(original_file_path):
    header_list = []
    original_file = open(original_file_path, 'r')
    file_line_list = original_file.readlines()
    for line in file_line_list:
        if line != '' and line[0] == '>':
            header_list.append(line)
        else:
            continue
    return header_list

def list_to_string(the_list): 
    string = ""
    for i in range(0, len(the_list)):
        string += the_list[i]
    return string

def get_gene_header_list(header_list):
    gene_header_list = []
    gene = 'other'
    for line in header_list:
        header_genes = []
        for g in gene_name_list:
            if g in line:
                header_genes.append(g)
        if len(header_genes) == 1:
            gene = header_genes[0]
            if gene not in gene_header_list:
                gene_header_list.append(gene)
        else:
            gene = list_to_string(header_genes)
            if gene not in gene_header_list:
                gene_header_list.append(gene)
    return gene_header_list

def get_gene_number_list(gene_header_list):
    gene_number_list = []
    for i in range(1, len(gene_header_list)+1):
        gene_number = str(i)
        gene_number = gene_number.zfill(3)
        gene_number_list.append(gene_number)
    return gene_number_list

def get_csv_file(gene_header_list, gene_number_list, csv_file_path):
    csv_file = open(csv_file_path, "w")
    writer = csv.writer(csv_file)
    for i in range(0, len(gene_number_list)):
        number = gene_number_list[i]
        header = gene_header_list[i]
        row = number,header
        writer.writerow(row)
    csv_file.close()
    return 'csv file created'

def format_header(line, gene_header_list, gene_number_list):
    AccessNumber = line[1:9]
    header_genes = []
    gene_number = 'other'
    for g in gene_name_list:
        if g in line:
            header_genes.append(g)
        if len(header_genes) == 1:
            gene = header_genes[0]
            if gene in gene_header_list:
                i = gene_header_list.index(gene)
                gene_number = gene_number_list[i]
        else:
            gene = list_to_string(header_genes)
            if gene in gene_header_list:
                i = gene_header_list.index(gene)
                gene_number = gene_number_list[i]
    return ('>' + AccessNumber + '-gene{}'.format(gene_number))


def rename_targets(origninal_file_path, csv_file_path):
    header_list = get_headers(original_file_path)
    gene_header_list = get_gene_header_list(header_list)
    gene_number_list = get_gene_number_list(gene_header_list)
    get_csv_file(gene_header_list, gene_number_list, csv_file_path)
    new_file = open('Renames_Target_file2.fasta', 'a')
    original_file = open(original_file_path, 'r')
    file_line_list = original_file.readlines()
    for line in file_line_list:
        if line != '' and line[0] == '>':
            new_header = format_header(line, gene_header_list, gene_number_list)
            new_file.write(new_header)
            new_file.write('\n')
        else:
            new_file.write(line)
    original_file.close()
    return 'Done :-D'

original_file_path = "C:\\Users\\Sarah\\Documents\\AU\\6\\Bachelor\\GitHub\\BSc\\Target_filer\\Original_Target_file.fasta"
csv_file_path = "C:\\Users\\Sarah\\Documents\\AU\\6\\Bachelor\\GitHub\\BSc\\Renaming_csv_files\\Rename_Targets.csv"

print(rename_targets(original_file_path, csv_file_path))