'''
------------------------------------------------------------------------------------------------------------------------
This workflow is used to run HybPiper on GenomeDK to find chloroplast DNA in the data. 
------------------------------------------------------------------------------------------------------------------------
This code is a variant of oscar_eksempel_species_workflow.py by Oscar Wrisberg
Eddited by Sarah E.K. Kessel 
------------------------------------------------------------------------------------------------------------------------
'''

from os import O_SYNC, name
from gwf import Workflow
import os.path
import csv

gwf = Workflow()

########################################################################################################################
################################################---- Hybpiper ----######################################################
########################################################################################################################
def hybpiper(species, p1, p2, un, path_out, path_in, done):
    """Hybpiper."""
    inputs = [path_in + species + p1, path_in + species + p2, path_in + species + un] # The files which the job will look for before it runs
    outputs = [path_out + species, done] # The files which will have to be created in order for the job to be "completed"
    options = {'cores': 1, 'memory': "20g", 'walltime': "100:00:00"} #'account':"Coryphoideae"} #Slurm commands

    spec = """
    source activate base

    cd {out}
        
    /home/sarahe/HybPiper/reads_first.py --cpu 16 --readfiles {p1} {p2} --unpaired {un} -b /home/sarahe/GitHub/BSc/Target_filer/Wolf_Target.fasta --prefix {species} --bwa
    touch {done}
    """.format(species=species, p1 = path_in + species + p1, p2 = path_in + species + p2, un = path_in + species + un , out = path_out, done = done)


    return (inputs, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################

sp = []

def read_csv(file_name, file_delimiter):
    name_list = []
    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=file_delimiter)
        for row in csv_reader:
            for i in range(0, len(ending)):
                file_path = ("/home/sarahe/BSc/00_data/"+str(row[1])+"_clean-Read1.fastq")
                isFile = os.path.isfile(file_path)
                if (isFile == True) and row[1] not in name_list:
                    name_list.append(row[1])
                else:
                    print('could not find' + file_path)
    return name_list

rename = "/home/sarahe/GitHub/BSc/Renaming_csv_files/Rename_Files.csv"

sp = read_csv(rename, ';')
print(len(sp))
print(sp)

#for i in range(len(sp)):
#    gwf.target_from_template('Hybpiper_'+str(i), hybpiper(species = sp[i],
#                                                        p1 = "_clean-Read1.fastq",
#                                                        p2 = "_clean-Read2.fastq",
#                                                        un = "_clean-Read12-single.fastq",
#                                                        path_out = "/home/sarahe/BSc/01_HybPiper_wolf_test/",
#                                                        path_in = "/home/sarahe/BSc/00_data/",
#                                                        done = "/home/sarahe/BSc/01_HybPiper_wolf_test/done/Hybpiper/"+sp[i]))