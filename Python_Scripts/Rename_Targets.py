## This script has been written with the purpouse of taking a Fasta file with multiple sequences and renaming them so the file can be used as a target file for HybPiper.

gene_list = ['accD', 'atpA', 'atpB', 'atpF', 'atpH', 'atpI', 'ccsA', 'clpP', 'matK', 'ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'petA', 'petB', 'petD', 'petL', 'petN', 'psaA', 'psaB', 'psaC', 'psaI', 'psbA', 'psbB', 'psbC', 'psbD', 'psbK', 'psbL', 'psbZ', 'psI', 'rbcL', 'rpl16', 'rpl2', 'rpl22', 'rpl23', 'rpl32', 'rpoA', 'rpoB', 'rpoC1', 'rpoC2', 'rps11', 'rps12', 'rps15', 'rps16', 'rps18', 'rps2', 'rps3', 'rps4', 'rps7', 'rps8', 'rrn16', 'rrn23', 'rrn4,5', 'rrn5', 'rrn5 Dio', 'trnA', 'trnC', 'trnD', 'trnfM', 'trnG', 'trnH', 'trnI', 'trnK', 'trnL', 'trnN', 'trnP', 'trnQ', 'trnS', 'trnT', 'trnV', 'ycf1', 'ycf2', 'ycf3', 'ycf4', 'ycl2']

def format_header(line):
    AccessNumber = line[1:9]
    header_genes = []
    for gene in gene_list:
        if gene in line:
            header_genes.append(gene)
    if len(header_genes) == 1:
        return ('>' + AccessNumber + '-' + header_genes[0])
    if len(header_genes) >= 1:
        genes = ''
        for e in range(len(header_genes)):
            if e == 0:
                genes = header_genes[e]
            else:
                genes = genes + '&' + header_genes[e]
        return ('>' + AccessNumber + '-' + genes)
    return ('>' + AccessNumber + '-other')

def rename_targets(origninal_file_path):
    new_file = open('Renames_Target_file.fasta', 'a')
    original_file = open(original_file_path, 'r')
    file_line_list = original_file.readlines()
    for line in file_line_list:
        if line != '' and line[0] == '>':
            new_header = format_header(line)
            new_file.write(new_header)
            new_file.write('\n')
        else:
            new_file.write(line)
    original_file.close()
    return ':-S'

original_file_path = "C:\\Users\\Sarah\\Documents\\AU\\6\\Bachelor\\GitHub\\BSc\\Target_filer\\Original_Target_file.fasta"

print(rename_targets(original_file_path))