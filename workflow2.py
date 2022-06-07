'''
------------------------------------------------------------------------------------------------------------------------
This workflow is used to run HybPiper on GenomeDK to find chloroplast DNA in the data. 
------------------------------------------------------------------------------------------------------------------------
This code is a variant of oscar_eksempel_species_workflow.py by Oscar Wrisberg
Eddited by Sarah E. Valentin
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
    options = {'cores': 1, 'memory': "20g", 'walltime': "100:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"} #Slurm commands

    spec = """
    source /home/sarahe/miniconda3/etc/profile.d/conda.sh
    source activate base

    cd {out}

    /home/sarahe/HybPiper/reads_first.py --cpu 16 --readfiles {p1} {p2} --unpaired {un} -b /home/sarahe/GitHub/BSc/Target_filer/Renamed_Target_file4.fasta --prefix {species} --bwa
    mkdir -p done
    cd done
    touch {species}
    """.format(species=species, p1 = path_in + species + p1, p2 = path_in + species + p2, un = path_in + species + un , out = path_out, done = done)


    return (inputs, outputs, options, spec)

########################################################################################################################
#############################################---- Name list ----########################################################
########################################################################################################################
def get_namelist(done_path, name_path):
    """Get namelist file."""
    inputs = []
    outputs = [name_path + 'name_list.txt']
    options = {'cores': 1, 'memory': "20g", 'walltime': "1:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    ls {done_path} > {name_path}name_list.txt
    """.format(done_path = done_path, name_path = name_path)

    return (inputs, outputs, options, spec)    

########################################################################################################################
#############################################---- Seq lengths ----######################################################
########################################################################################################################
def seq_lenghts(path):
    """Get the sequence lengths."""
    inputs = [path + 'name_list.txt']
    outputs = [path + 'seq_lengths.txt']
    options = {'cores': 1, 'memory': "20g", 'walltime': "1:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    cd /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper
    python /home/sarahe/HybPiper/get_seq_lengths.py /home/sarahe/GitHub/BSc/Target_filer/Renamed_Target_file4.fasta {path}name_list.txt dna > {path}seq_lengths.txt
    """.format(path = path)

    return (inputs, outputs, options, spec)

########################################################################################################################
#############################################---- Stats summary ----####################################################
########################################################################################################################
def stats_summary(path):
    """Get statistical summary."""
    inputs = [path + 'seq_lengths.txt', path + 'name_list.txt']
    outputs = [path + 'stats.txt']
    options = {'cores': 1, 'memory': "20g", 'walltime': "1:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    python /home/sarahe/HybPiper/hybpiper_stats.py {path}seq_lengths.txt {path}name_list.txt > {path}stats.txt
    """.format(path = path)

    return (inputs, outputs, options, spec)

########################################################################################################################
#############################################---- Intronerate ----######################################################
########################################################################################################################

def intronerate(species, path_in, done):
    """Intronerate the sequencec from hybpiper."""
    inputs = [path_in + species]
    outputs = [done]
    options = {'cores': 4, 'memory': "20g", 'walltime': "16:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    source /home/sarahe/miniconda3/etc/profile.d/conda.sh
    source activate base

    cd {path_in}

    python3 /home/sarahe/HybPiper/intronerate.py --prefix {sp} &>> intronerate_out.txt
        
    
    touch {done}
    """.format(sp = species, done = done, path_in = path_in)

    return (inputs, outputs, options, spec)

########################################################################################################################
#############################################---- Coverage ----#########################################################
########################################################################################################################
def coverage(species, dir_in, dir_out, dir_wrk, path_in, path_out, done, all_bam, all_sorted_bam, all_sorted_bam_bai, bam, cov, fasta, fasta_amb, fasta_ann, fasta_bwt, fasta_pac, fasta_sa, trimmed_fasta, up_bam):
    """Calculating coverage of sequences."""
    inputs = [path_in + 'seq_lengths.txt', path_in + 'name_list.txt']
    outputs = [path_out+species+all_bam, path_out+species+all_sorted_bam, path_out+species+all_sorted_bam_bai, path_out+species+bam,
    path_out+species+cov, path_out+species+fasta, path_out+species+fasta_amb, path_out+species+fasta_ann, path_out+species+fasta_bwt,
    path_out+species+fasta_pac, path_out+species+fasta_sa, path_out+species+trimmed_fasta, path_out+species+up_bam, done]
    options = {'cores': 4, 'memory': "20g", 'walltime': "08:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    source /home/sarahe/miniconda3/etc/profile.d/conda.sh
    source activate base
    
    cd {path_in}

    python3 /home/sarahe/GitHub/BSc/Python_Scripts/coverage_eddit.py {sp} {dir_in} {dir_out} {dir_wrk}
    
    touch {done}

    """.format(sp = species, done = done, path_in = path_in, dir_in = dir_in, dir_out = dir_out, dir_wrk = dir_wrk)

    return (inputs, outputs, options, spec)

##########################################################################################################################
###########################################---- Retrieve Sequences ----###################################################
##########################################################################################################################

def retrieve(path_in, path_python, path_out, done):
    """Retrieve gene sequences from all the species and create an unaligned multifasta for each gene."""
    inputs = []
    outputs = [done]
    options = {'cores': 10, 'memory': "20g", 'walltime': "12:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    source /home/sarahe/miniconda3/etc/profile.d/conda.sh
    source activate base

    cd {path_in}

    ls *trimmed.fasta > filelist.txt

    python {psth_python}samples2genes.py > {path_out}outstats.csv

    touch {path_out}Retrieve_all_done.txt

    """.format(path_in = path_in, path_python = path_python, path_out = path_out, done = done)

    return (inputs, outputs, options, spec)

##########################################################################################################################
###############################################---- MAFFT ----#############################################################
##########################################################################################################################

def mafft(gene, path_in, path_out, done):
    """Aligning all the sequences for each gene."""
    inputs = ["/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/done/Retrieve_Genes/Retrieve_all_done.txt", path_in+gene+".FNA"]
    outputs = [done,path_out+gene+"_aligned.fasta"] 
    options = {'cores': 1, 'memory': "500g", 'walltime': "10:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """

    source activate mafft_env

    cd {path_in}

    mafft --globalpair --large --adjustdirectionaccurately --thread 1 {gene}.FNA > {path_out}{gene}_aligned.fasta
    
    touch {done}

    """.format(gene = gene, done = done, path_in = path_in, path_out=path_out)

    return (inputs, outputs, options, spec)

##########################################################################################################################
#######################################---- Post MAFFT Renaming ----######################################################
##########################################################################################################################

def post_mafft(gene, path_in, path_out, path_python, done):
    """Aligning all the sequences for each gene."""
    inputs = [path_in+gene+"_aligned.fasta"]
    outputs = [done+gene, path_out+gene+"_renamed_aligned.fasta"] 
    options = {'cores': 1, 'memory': "20g", 'walltime': "0:15:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    source /home/sarahe/miniconda3/etc/profile.d/conda.sh
    source activate base

    python {path_python}post_mafft_renaming.py {path_in}{gene}_aligned.fasta > {path_out}{gene}_renamed_aligned.fasta

    touch {done}{gene}

    """.format(gene = gene, done = done, path_in = path_in, path_out=path_out, path_python=path_python)

    return (inputs, outputs, options, spec)

##########################################################################################################################
###############################################---- IQtree ----###########################################################
##########################################################################################################################

def iqtree(done, path_out, path_in):
    """Runs IQtree from the concoctanated file from Trim"""
    inputs = []
    outputs = [done+'iqtree.txt'] 
    options = {'cores': 5, 'memory': "10g", 'walltime': "24:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    source /home/sarahe/miniconda3/etc/profile.d/conda.sh
    source activate base

    cd {path_out}

    iqtree -s {path_in} -m MFP -nt 12 -T AUTO -B 1000 --redo

    touch {done}iqtree.txt

    """.format(done=done, path_out=path_out, path_in=path_in)

    return (inputs, outputs, options, spec)

########################################################################################################################
######################################################---- RUN ----#####################################################
########################################################################################################################

# Create empty species list
sp = []

# function that creates a full namelist from csv file
def read_csv(file_name, file_delimiter):
    name_list = []
    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=file_delimiter)
        for row in csv_reader:
            file_path = ("/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/00_data/"+str(row[1])+"_clean-Read1.fastq")
            isFile = os.path.isfile(file_path)
            if (isFile == True) and row[1] not in name_list and row[1]:
                name_list.append(row[1])
    return name_list

# define csv file placement
rename = "/home/sarahe/GitHub/BSc/Renaming_csv_files/Rename_00_data_files.csv"

# create full species list
sp = read_csv(rename, ';')
#print(sp)

# # run hybpiper
# for i in range(0, len(sp)):
#     gwf.target_from_template('hybpiper_'+sp[i], hybpiper(species = sp[i],
#                                                         p1 = "_clean-Read1.fastq",
#                                                         p2 = "_clean-Read2.fastq",
#                                                         un = "_clean-Read12-single.fastq",
#                                                         path_out = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/",
#                                                         path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/00_data/",
#                                                         done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/done/"+sp[i]))

# # get name list from hybpiper run
# gwf.target_from_template('name_list', get_namelist(done_path = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/done/",
#                                                     name_path = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/"))

# # get sequence length file necessary for statistical summary
# gwf.target_from_template('sequence_length', seq_lenghts(path = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/"))

# # get statistics, which can be used with gene_coverage_gg_plot.R to visualise results
# gwf.target_from_template('statistics', stats_summary(path = '/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/'))

# # run intronerate
# for i in range(0, len(sp)):
#     gwf.target_from_template('intronerate_'+sp[i], intronerate(species= sp[i],
#                                                         path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/",
#                                                         done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/done/Intronerate/"+sp[i]))

# run Coverage to estimate the significance of the contigs found by hybpiper
# for i in range(0, len(sp)):
#     gwf.target_from_template('coverage_'+sp[i], coverage(species = sp[i],
#                                                         dir_in = '/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/00_data/',
#                                                         dir_out = '/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/',
#                                                         dir_wrk = '/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/',
#                                                         path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/",
#                                                         all_bam = "_all.bam",
#                                                         all_sorted_bam ="_all_sorted.bam",
#                                                         all_sorted_bam_bai="_all_sorted.bam.bai",
#                                                         bam =".bam",
#                                                         cov=".cov",
#                                                         fasta = ".fasta",
#                                                         fasta_amb = ".fasta.amb",
#                                                         fasta_ann = ".fasta.ann",
#                                                         fasta_bwt = ".fasta.bwt",
#                                                         fasta_pac = ".fasta.pac",
#                                                         fasta_sa = ".fasta.sa",
#                                                         trimmed_fasta = "_trimmed.fasta",
#                                                         up_bam = "_up.bam",
#                                                         path_out = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/",
#                                                         done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/done/Coverage/"+sp[i]))
# Define genes
genes = ['accD', 'atpA', 'atpB', 'atpF', 'atpH', 'atpI', 'ccsA', 'clpP', 'matK', 'ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'petA', 'petB', 'petD', 'petL', 'petN', 'psaA', 'psaB', 'psaC', 'psaI', 'psbA', 'psbB', 'psbC', 'psbD', 'psbK', 'psbL', 'psbZ', 'psI', 'rbcL', 'rpl16', 'rpl2', 'rpl22', 'rpl23', 'rpl32', 'rpoA', 'rpoB', 'rpoC1', 'rpoC2', 'rps11', 'rps12', 'rps15', 'rps16', 'rps18', 'rps2', 'rps3', 'rps4', 'rps7', 'rps8', 'rrn16', 'rrn23', 'rrn4,5', 'rrn5', 'rrn5 Dio', 'trnA', 'trnC', 'trnD', 'trnfM', 'trnG', 'trnH', 'trnI', 'trnK', 'trnL', 'trnN', 'trnP', 'trnQ', 'trnS', 'trnT', 'trnV', 'ycf1', 'ycf2', 'ycf3', 'ycf4', 'ycl2']

# Define gt values
gt_values =["0.1","0.15","0.2","0.25","0.3","0.33","0.4","0.45","0.5","0.55","0.6","0.67","0.7","0.75","0.8","0.85","0.9","0.95"]

# # Retrieve sequences and sort into files with gene names
#gwf.target_from_template('Retrieve_genes', retrieve(path_in="/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/"))

# Running MAFFT
# for i in range(len(genes)):
#     pth = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/03_blacklisting/"+genes[i]+'.FNA'
#     if os.path.isfile(pth):
#         gwf.target_from_template('Mafft_'+str(i), mafft(gene = genes[i],
#         path_out= "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/04_mafft/",
#         path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/03_blacklisting/",
#         done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/04_mafft/done/"+genes[i]))

# Renaming mafft output files to be ready for IQtree
# for i in range(0, len(genes)):
#     pth = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/04_mafft/"+genes[i]+'_aligned.fasta'
#     if os.path.isfile(pth):
#         gwf.target_from_template('post_mafft_'+str(i), post_mafft(gene = genes[i],
#                                                                     path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/04_mafft/",
#                                                                     path_out = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/05_post_mafft/",
#                                                                     path_python = "/home/sarahe/GitHub/BSc/Python_Scripts/",
#                                                                     done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/05_post_mafft/done/"))

# Running IQtree
gwf.target_from_template('IQtree', iqtree(path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/05_post_mafft/fasta/",
                                        path_out = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/06_iqtree/",
                                        done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/06_iqtree/done/"))

#################################################################################################################################
########################################### Running workflow for good gene recovery #############################################
#################################################################################################################################

sp2 = ['LORU1', 'DSCO0', 'DSCA0', 'DRAB0', 'DPUS1', 'DPRE0', 'DOVO1', 'DMET0', 'DHIL0', 'DHIA0', 'DDRA1', 'DCON1', 'DBON0', 'DBAR0']

# # run hybpiper
# for i in range(0, len(sp2)):
#     gwf.target_from_template('hybpiper_'+sp2[i], hybpiper(species = sp2[i],
#                                                         p1 = "_clean-Read1.fastq",
#                                                         p2 = "_clean-Read2.fastq",
#                                                         un = "_clean-Read12-single.fastq",
#                                                         path_out = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/07_HybPiper2/",
#                                                         path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/00_data/",
#                                                         done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/07_HybPiper2/done/"+sp2[i]))

# # get name list from hybpiper run
# gwf.target_from_template('name_list', get_namelist(done_path = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/07_HybPiper2/done/",
#                                                     name_path = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/07_HybPiper2/"))

# # get sequence length file necessary for statistical summary
# gwf.target_from_template('sequence_length', seq_lenghts(path = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/07_HybPiper2/"))

# # get statistics, which can be used with gene_coverage_gg_plot.R to visualise results
# gwf.target_from_template('statistics', stats_summary(path = '/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/07_HybPiper2/'))

# # run intronerate
# for i in range(0, len(sp2)):
#     gwf.target_from_template('intronerate_'+sp2[i], intronerate(species= sp2[i],
#                                                         path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/07_HybPiper2/",
#                                                         done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/07_HybPiper2/done/Intronerate/"+sp2[i]))

# # run Coverage to estimate the significance of the contigs found by hybpiper
# for i in range(0, len(sp2)):
#     gwf.target_from_template('coverage_'+sp2[i], coverage(species = sp2[i],
#                                                         dir_in = '/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/00_data/',
#                                                         dir_out = '/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/08_Coverage2/',
#                                                         dir_wrk = '/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/',
#                                                         path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/07_HybPiper2/",
#                                                         all_bam = "_all.bam",
#                                                         all_sorted_bam ="_all_sorted.bam",
#                                                         all_sorted_bam_bai="_all_sorted.bam.bai",
#                                                         bam =".bam",
#                                                         cov=".cov",
#                                                         fasta = ".fasta",
#                                                         fasta_amb = ".fasta.amb",
#                                                         fasta_ann = ".fasta.ann",
#                                                         fasta_bwt = ".fasta.bwt",
#                                                         fasta_pac = ".fasta.pac",
#                                                         fasta_sa = ".fasta.sa",
#                                                         trimmed_fasta = "_trimmed.fasta",
#                                                         up_bam = "_up.bam",
#                                                         path_out = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/08_Coverage2/",
#                                                         done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/08_Coverage2/done/Coverage/"+sp2[i]))

# Define genes
genes = ['accD', 'atpA', 'atpB', 'atpF', 'atpH', 'atpI', 'ccsA', 'clpP', 'matK', 'ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhF', 'ndhG', 'ndhH', 'ndhI', 'ndhJ', 'petA', 'petB', 'petD', 'petL', 'petN', 'psaA', 'psaB', 'psaC', 'psaI', 'psbA', 'psbB', 'psbC', 'psbD', 'psbK', 'psbL', 'psbZ', 'psI', 'rbcL', 'rpl16', 'rpl2', 'rpl22', 'rpl23', 'rpl32', 'rpoA', 'rpoB', 'rpoC1', 'rpoC2', 'rps11', 'rps12', 'rps15', 'rps16', 'rps18', 'rps2', 'rps3', 'rps4', 'rps7', 'rps8', 'rrn16', 'rrn23', 'rrn4,5', 'rrn5', 'rrn5 Dio', 'trnA', 'trnC', 'trnD', 'trnfM', 'trnG', 'trnH', 'trnI', 'trnK', 'trnL', 'trnN', 'trnP', 'trnQ', 'trnS', 'trnT', 'trnV', 'ycf1', 'ycf2', 'ycf3', 'ycf4', 'ycl2']

# Define gt values
gt_values =["0.1","0.15","0.2","0.25","0.3","0.33","0.4","0.45","0.5","0.55","0.6","0.67","0.7","0.75","0.8","0.85","0.9","0.95"]

# Retrieve sequences and sort into files with gene names
gwf.target_from_template('Retrieve_genes', retrieve(path_in="/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/08_Coverage2/",
                                                        path_python = '/home/sarahe/GitHub/BSc/Python_Scripts/', 
                                                        path_out = '/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/done/Retrieve_Genes/',
                                                        done = path_out+'Retrieve_all_done.txt'))

# Running MAFFT
# for i in range(len(genes)):
#     pth = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/09_blacklist2/"+genes[i]+'.FNA'
#     if os.path.isfile(pth):
#         gwf.target_from_template('Mafft_'+str(i), mafft(gene = genes[i],
#         path_out= "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/10_mafft2/",
#         path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/09_blacklist2/",
#         done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/10_mafft2/done/"+genes[i]))

# # Renaming mafft output files to be ready for IQtree
# for i in range(len(genes)):
#     pth = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/10_mafft2/"+genes[i]+'_aligned.fasta'
#     if os.path.isfile(pth):
#         gwf.target_from_template('post_mafft_'+str(i), post_mafft(gene = genes[i],
#                                                                     path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/10_mafft2/",
#                                                                     path_out = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/11_post_mafft2/",
#                                                                     path_python = "/home/sarahe/GitHub/BSc/Python_Scripts/",
#                                                                     done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/11_post_mafft2/done/"))

# # Running IQtree
# gwf.target_from_template('IQtree', iqtree(path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/11_post_mafft2/fasta/",
#                                         path_out = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/12_iqtree2/",
#                                         done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/12_iqtree2/done/"))
