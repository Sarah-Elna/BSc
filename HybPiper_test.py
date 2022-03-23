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
        
"""/home/sarahe/HybPiper/reads_first.py --cpu 1 --readfiles {p1} {p2} --unpaired {un} -b ##/home/owrisberg/Coryphoideae/target_sequence/PhyloPalms_loci_renamed_794-176_HEYcorrected.fasta --prefix {species} --bwa
"""
    touch {done}
    """.format(species=species, p1 = path_in + species + paired_1, p2 = path_in + species + paired_2, un = path_in + species + unpaired , out = path_out, done = done_file)


    return (inputs, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################

sp = []

for line in names.txt:
    sp.append(line)
    return sp

print(sp)

#for i in range(len(sp)):
    #### Running Hybpiper
#    gwf.target_from_template('Hybpiper_'+sp[i], hybpiper(species = sp[i],
#                                                        p1 = #"_1P.fastq",
#                                                        p2 = #"_2P.fastq",
#                                                        un = #"_UN.fastq",
#                                                        path_out= #"/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/",
#                                                        path_in = #"/home/owrisberg/Coryphoideae/work_flow/02_trimmed/",
#                                                        done = #"/home/owrisberg/Coryphoideae/work_flow/03_hybpiper/done/Hybpiper/"+sp[i]))