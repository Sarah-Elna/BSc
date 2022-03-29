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

gwf = Workflow()

########################################################################################################################
################################################---- Hybpiper ----######################################################
########################################################################################################################
def hybpiper(species, p1, p2, un, path_out, path_in, done):
    """Hybpiper."""
    inputs = [path_in + species +p1, path_in + species + p2, path_in + species + un] # The files which the job will look for before it runs
    outputs = [path_out + species, done] # The files which will have to be created in order for the job to be "completed"
    options = {'cores': 1, 'memory': "20g", 'walltime': "100:00:00"} #'account':"Coryphoideae"} #Slurm commands

    spec = """
    source activate base

    cd {out}
        
    /home/sarahe/HybPiper/reads_first.py --cpu 1 --readfiles {p1} {p2} --unpaired {un} -b /home/sarahe/GitHub/BSc/matK_rbcL_psbA_target.fasta --prefix {species} --bwa

    touch {done}
    """.format(species=species, p1 = path_in + species + p1,p2 = path_in + species + p2, un = path_in + species + un , out = path_out, done = done)


    return (inputs, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################

sp = []

#def read_csv(file_name, file_delimiter):
#    name_list = []
#    with open(file_name) as csv_file:
#        csv_reader = csv.reader(csv_file, delimiter=file_delimiter)
#        for row in csv_reader:
#            name_list.append(row[1])
#    return (name_list)

#rename = "C://Users//Sarah//Documents//AU//6//Bachelor//GenomeDK//rename.csv"

## sp = read_csv(rename, ';')

for i in range(len(sp)):
    gwf.target_from_template('Hybpiper_'+sp[i], hybpiper(species = sp[i],
                                                        paired_1 = "_clean-READ1.fastq",
                                                        paired_2 = "_clean-READ2.fastq",
                                                        unpaired = "_clean-READ12-single.fastq",
                                                        path_out = "/home/sarahe/BSc/01_HybPiper/",
                                                        path_in = "/home/sarahe/BSc/00_data/",
                                                        done_file = "/home/sarahe/BSc/01_HybPiper/done/Hybpiper/"+sp[i]))