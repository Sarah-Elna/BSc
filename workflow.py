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
def hybpiper(species, paired_1, paired_2, unpaired, path_out, path_in, done_file):
    """Hybpiper."""
    inputs = [path_in + species + paired_1, path_in + species + paired_2, path_in + species + unpaired] # The files which the job will look for before it runs
    outputs = [path_out + species, done_file] # The files which will have to be created in order for the job to be "completed"
    options = {'cores': 1, 'memory': "20g", 'walltime': "100:00:00"} #Slurm commands

    spec = """
    source activate base

    cd {out}
        
/home/sarahe/HybPiper/reads_first.py --cpu 16 -r {p1} {p2} --unpaired {un} -b /home/sarahe/GitHub/BSc/matK_rbcL_psbA_target.fasta --prefix {sp} --bwa

    touch {done}
    """.format(sp=species, p1 = path_in + species + paired_1, p2 = path_in + species + paired_2, un = path_in + species + unpaired , out = path_out, done = done_file)


    return (inputs, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################

#sp = ['0001', '0002', '0003', '0004', '0005', '0006', '0007', '0008', '0009', '0010', '0011', '0012', '0013', '0014', '0015', '0016', '0017', '0018', '0019', '0020']

sp =[]

names = open('names.txt', 'r')
for line in names:
    line_stripped = line.strip()
    if line_stripped not in sp and (line_stripped != 'file') and (line_stripped != 'READ'):
        sp.append(line_stripped)

print(sp)

#for i in range(len(sp)):
#    gwf.target_from_template('Hybpiper_'+sp[i], hybpiper(species = sp[i],
#                                                        paired_1 = "_clean-READ1.fastq",
#                                                        paired_2 = "_clean-READ2.fastq",
#                                                        unpaired = "_clean-READ12-single.fastq",
#                                                        path_out = "/home/sarahe/BSc/01_HybPiper/",
#                                                        path_in = "/home/sarahe/BSc/00_data/",
#                                                        done_file = "/home/sarahe/BSc/01_HybPiper/done/Hybpiper/"+sp[i]))