gene_recovered_species_one = ['DACA0','DACU0','DAFF0','DALB0','DAMB0','DAMB1','DAMB2','DAND0','DAND1','DAND2','DANG0','DANG1','DANK0','DANT0','DAQU0','DARE0','DBAR0','DBAS0','DBEE0','DBEJ0','DBER0','DBET0','DBET1','DBET2','DBOI0','DBON0','DBOS0','DBRE0','DBRI0','DCAB0','DCAN0','DCAR0','DCAT0','DCAU0','DCER0','DCOM0','DCON0','DCON1','DCOO0','DCOR0','DCOR1','DCOU0','DCRI0','DCUL0','DCUR0','DDEC0','DDEC1','DDEL0','DDIG0','DDRA0','DDRA1','DELE0','DERI0','DFAN0','DFAN1','DFAS0','DFIB0','DFOR0','DGEU0','DGLA0','DGRO0','DHEN0','DHET0','DHET1','DHET2','DHET3','DHIA0','DHIL0','DHOV0','DHUM0','DHUM1','DHUM2','DIFA0','DINT0','DINT1','DINT2','DJER0','DJUM0','DJUR0','DLAE0','DLAN0','DLAN1','DLAN2','DLAN3','DLAN4','DLAS0','DLEP0','DLEU0','DLEU1','DLIL0','DLIN0','DLOK0','DLOU0','DLUT0','DLUT1','DMAD0','DMAD1','DMAH0','DMAK0','DMAL0','DMAN0','DMAN1','DMAR0','DMCD0','DMET0','DMIN0','DMIR0','DMIR1','DMOC0','DMON0','DMOO0','DNAU0','DNOD0','DNOS0','DOCC0','DONI0','DONI1','DORE0','DORO0','DOVO0','DOVO1','DPAC0','DPAL0','DPER0','DPER1','DPIL0','DPIN0','DPLU0','DPOI0','DPRE0','DPRO0','DPRO1','DPSA0','DPUL0','DPUM0','DPUS0','DPUS1','DRAB0','DRAK0','DRAM0','DREF0','DREM0','DRIV0','DROB0','DROS0','DSAH0','DSAI0','DSAI1','DSAN0','DSAN1','DSCA0','DSCH0','DSCO0','DSER0','DSIM0','DSIN0','DSPB0','DSPI0','DSPN0','DSPN1','DSUB0','DTAN0','DTEN0','DTHE0','DTHI0','DTOK0','DTRA0','DTSA0','DTUR0','DUTI0','DVIR0','DVON0','LEHA0','LESP0','LORU0','LORU1','MRDA0','MRDA1','MRDA2','MRDA3','MRIN0','MRIN1','MRIN2','MRSP0','MSKO0','MSKO1','MSMA0']

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

rename = "/mnt/c/Users/Sarah/Documents/AU/6/Bachelor/GitHub/BSc/Renaming_csv_files/Rename_00_data_files.csv"
id_list, name_list = read_csv(rename, ';')
gene_recovered_species_one_translated = []

for i in range(0, len(gene_recovered_species_one)):
    for j in range(0, len(id_list)):
        if gene_recovered_species_one[i] == id_list[j]:
            if name_list[j] not in gene_recovered_species_one_translated:
                gene_recovered_species_one_translated.append(name_list[j])

print(gene_recovered_species_one_translated, len(gene_recovered_species_one_translated))
