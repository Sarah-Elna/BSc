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
    inputs = [done_path+'DACA0',done_path+'DANG1',done_path+'DBET1',done_path+'DCAU0',done_path+'DCUR0',done_path+'DFAS0',done_path+'DHIA0',done_path+'DJUM0',done_path+'DLEU1',done_path+'DMAL0',done_path+'DMOO0',done_path+'DPAC0',done_path+'DPSA0',done_path+'DROB0',done_path+'DSIM0',done_path+'DTOK0',done_path+'MRDA0',done_path+'DACU0',done_path+'DANK0',done_path+'DBET2',done_path+'DCER0',done_path+'DDEC0',done_path+'DFIB0',done_path+'DHIL0',done_path+'DJUR0',done_path+'DLIL0',done_path+'DMAN0',done_path+'DNAU0',done_path+'DPAL0',done_path+'DPUL0',done_path+'DROS0',done_path+'DSIN0',done_path+'DTRA0',done_path+'MRDA1',done_path+'DAFF0',done_path+'DANT0',done_path+'DBOI0',done_path+'DCOM0',done_path+'DDEC1',done_path+'DFOR0',done_path+'DHOV0',done_path+'DLAE0',done_path+'DLIN0',done_path+'DMAN1',done_path+'DNOD0',done_path+'DPER0',done_path+'DPUM0',done_path+'DSAH0',done_path+'DSPB0',done_path+'DTSA0',done_path+'MRDA2',done_path+'DALB0',done_path+'DAQU0',done_path+'DBON0',done_path+'DCON0',done_path+'DDEL0',done_path+'DGEU0',done_path+'DHUM0',done_path+'DLAN0',done_path+'DLOK0',done_path+'DMAR0',done_path+'DNOS0',done_path+'DPER1',done_path+'DPUS0',done_path+'DSAI0',done_path+'DSPI0',done_path+'DTUR0',done_path+'MRDA3',done_path+'DAMB0',done_path+'DARE0',done_path+'DBOS0',done_path+'DCON1',done_path+'DDIG0',done_path+'DGLA0',done_path+'DHUM1',done_path+'DLAN1',done_path+'DLOU0',done_path+'DMCD0',done_path+'DOCC0',done_path+'DPIL0',done_path+'DPUS1',done_path+'DSAI1',done_path+'DSPN0',done_path+'DUTI0',done_path+'MRIN0',done_path+'DAMB1',done_path+'DBAR0',done_path+'DBRE0',done_path+'DCOO0',done_path+'DDRA0',done_path+'DGRO0',done_path+'DHUM2',done_path+'DLAN2',done_path+'DLUT0',done_path+'DMET0',done_path+'DONI0',done_path+'DPIN0',done_path+'DRAB0',done_path+'DSAN0',done_path+'DSPN1',done_path+'DVIR0',done_path+'MRIN1',done_path+'DAMB2',done_path+'DBAS0',done_path+'DBRI0',done_path+'DCOR0',done_path+'DDRA1',done_path+'DHEN0',done_path+'DIFA0',done_path+'DLAN3',done_path+'DLUT1',done_path+'DMIN0',done_path+'DONI1',done_path+'DPLU0',done_path+'DRAK0',done_path+'DSAN1',done_path+'DSUB0',done_path+'DVON0',done_path+'MRIN2',done_path+'DAND0',done_path+'DBEE0',done_path+'DCAB0',done_path+'DCOR1',done_path+'DELE0',done_path+'DHET0',done_path+'DINT0',done_path+'DLAN4',done_path+'DMAD0',done_path+'DMIR0',done_path+'DORE0',done_path+'DPOI0',done_path+'DRAM0',done_path+'DSCA0',done_path+'DTAN0',done_path+'LEHA0',done_path+'MRSP0',done_path+'DAND1',done_path+'DBEJ0',done_path+'DCAN0',done_path+'DCOU0',done_path+'DERI0',done_path+'DHET1',done_path+'DINT1',done_path+'DLAS0',done_path+'DMAD1',done_path+'DMIR1',done_path+'DORO0',done_path+'DPRE0',done_path+'DREF0',done_path+'DSCH0',done_path+'DTEN0',done_path+'LESP0',done_path+'MSKO0',done_path+'DAND2',done_path+'DBER0',done_path+'DCAR0',done_path+'DCRI0',done_path+'DFAN0',done_path+'DHET2',done_path+'DINT2',done_path+'DLEP0',done_path+'DMAH0',done_path+'DMOC0',done_path+'DOVO0',done_path+'DPRO0',done_path+'DREM0',done_path+'DSCO0',done_path+'DTHE0',done_path+'LORU0',done_path+'MSKO1',done_path+'DANG0',done_path+'DBET0',done_path+'DCAT0',done_path+'DCUL0',done_path+'DFAN1',done_path+'DHET3',done_path+'DJER0',done_path+'DLEU0',done_path+'DMAK0',done_path+'DMON0',done_path+'DOVO1',done_path+'DPRO1',done_path+'DRIV0',done_path+'DSER0',done_path+'DTHI0',done_path+'LORU1',done_path+'MSMA0']
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

def retrieve(path_in):
    """Retrieve gene sequences from all the species and create an unaligned multifasta for each gene."""
    inputs = [path_in + 'done/Coverage/DACA0', path_in + 'done/Coverage/DANG1', path_in + 'done/Coverage/DBET1', path_in + 'done/Coverage/DCAU0',path_in + 'done/Coverage/DCUR0', path_in + 'done/Coverage/DFAS0', path_in + 'done/Coverage/DHIA0', path_in + 'done/Coverage/DJUM0', path_in + 'done/Coverage/DLEU1', path_in + 'done/Coverage/DMAL0', path_in + 'done/Coverage/DMOO0', path_in + 'done/Coverage/DPAC0', path_in + 'done/Coverage/DPSA0', path_in + 'done/Coverage/DROB0', path_in + 'done/Coverage/DSIM0', path_in + 'done/Coverage/DTOK0', path_in + 'done/Coverage/MRDA0', path_in + 'done/Coverage/DACU0', path_in + 'done/Coverage/DANK0', path_in + 'done/Coverage/DBET2', path_in + 'done/Coverage/DCER0', path_in + 'done/Coverage/DDEC0', path_in + 'done/Coverage/DFIB0', path_in + 'done/Coverage/DHIL0', path_in + 'done/Coverage/DJUR0', path_in + 'done/Coverage/DLIL0', path_in + 'done/Coverage/DMAN0', path_in + 'done/Coverage/DNAU0', path_in + 'done/Coverage/DPAL0', path_in + 'done/Coverage/DPUL0', path_in + 'done/Coverage/DROS0', path_in + 'done/Coverage/DSIN0', path_in + 'done/Coverage/DTRA0', path_in + 'done/Coverage/MRDA1', path_in + 'done/Coverage/DAFF0', path_in + 'done/Coverage/DANT0', path_in + 'done/Coverage/DBOI0', path_in + 'done/Coverage/DCOM0', path_in + 'done/Coverage/DDEC1', path_in + 'done/Coverage/DFOR0', path_in + 'done/Coverage/DHOV0', path_in + 'done/Coverage/DLAE0', path_in + 'done/Coverage/DLIN0', path_in + 'done/Coverage/DMAN1', path_in + 'done/Coverage/DNOD0', path_in + 'done/Coverage/DPER0', path_in + 'done/Coverage/DPUM0', path_in + 'done/Coverage/DSAH0', path_in + 'done/Coverage/DSPB0', path_in + 'done/Coverage/DTSA0', path_in + 'done/Coverage/MRDA2', path_in + 'done/Coverage/DALB0', path_in + 'done/Coverage/DAQU0', path_in + 'done/Coverage/DBON0', path_in + 'done/Coverage/DCON0', path_in + 'done/Coverage/DDEL0', path_in + 'done/Coverage/DGEU0', path_in + 'done/Coverage/DHUM0', path_in + 'done/Coverage/DLAN0', path_in + 'done/Coverage/DLOK0', path_in + 'done/Coverage/DMAR0', path_in + 'done/Coverage/DNOS0', path_in + 'done/Coverage/DPER1', path_in + 'done/Coverage/DPUS0', path_in + 'done/Coverage/DSAI0', path_in + 'done/Coverage/DSPI0', path_in + 'done/Coverage/DTUR0', path_in + 'done/Coverage/MRDA3', path_in + 'done/Coverage/DAMB0', path_in + 'done/Coverage/DARE0', path_in + 'done/Coverage/DBOS0', path_in + 'done/Coverage/DCON1', path_in + 'done/Coverage/DDIG0', path_in + 'done/Coverage/DGLA0', path_in + 'done/Coverage/DHUM1', path_in + 'done/Coverage/DLAN1', path_in + 'done/Coverage/DLOU0', path_in + 'done/Coverage/DMCD0', path_in + 'done/Coverage/DOCC0', path_in + 'done/Coverage/DPIL0', path_in + 'done/Coverage/DPUS1', path_in + 'done/Coverage/DSAI1', path_in + 'done/Coverage/DSPN0', path_in + 'done/Coverage/DUTI0', path_in + 'done/Coverage/MRIN0', path_in + 'done/Coverage/DAMB1', path_in + 'done/Coverage/DBAR0', path_in + 'done/Coverage/DBRE0', path_in + 'done/Coverage/DCOO0', path_in + 'done/Coverage/DDRA0', path_in + 'done/Coverage/DGRO0', path_in + 'done/Coverage/DHUM2', path_in + 'done/Coverage/DLAN2', path_in + 'done/Coverage/DLUT0', path_in + 'done/Coverage/DMET0', path_in + 'done/Coverage/DONI0', path_in + 'done/Coverage/DPIN0', path_in + 'done/Coverage/DRAB0', path_in + 'done/Coverage/DSAN0', path_in + 'done/Coverage/DSPN1', path_in + 'done/Coverage/DVIR0', path_in + 'done/Coverage/MRIN1', path_in + 'done/Coverage/DAMB2', path_in + 'done/Coverage/DBAS0', path_in + 'done/Coverage/DBRI0', path_in + 'done/Coverage/DCOR0', path_in + 'done/Coverage/DDRA1', path_in + 'done/Coverage/DHEN0',  path_in + 'done/Coverage/DIFA0', path_in + 'done/Coverage/DLAN3', path_in + 'done/Coverage/DLUT1', path_in + 'done/Coverage/DMIN0', path_in + 'done/Coverage/DONI1', path_in + 'done/Coverage/DPLU0', path_in + 'done/Coverage/DRAK0', path_in + 'done/Coverage/DSAN1', path_in + 'done/Coverage/DSUB0', path_in + 'done/Coverage/DVON0', path_in + 'done/Coverage/MRIN2', path_in + 'done/Coverage/DAND0', path_in + 'done/Coverage/DBEE0', path_in + 'done/Coverage/DCAB0', path_in + 'done/Coverage/DCOR1', path_in + 'done/Coverage/DELE0', path_in + 'done/Coverage/DHET0', path_in + 'done/Coverage/DINT0', path_in + 'done/Coverage/DLAN4', path_in + 'done/Coverage/DMAD0', path_in + 'done/Coverage/DMIR0', path_in + 'done/Coverage/DORE0', path_in + 'done/Coverage/DPOI0', path_in + 'done/Coverage/DRAM0', path_in + 'done/Coverage/DSCA0', path_in + 'done/Coverage/DTAN0', path_in + 'done/Coverage/LEHA0', path_in + 'done/Coverage/MRSP0', path_in + 'done/Coverage/DAND1', path_in + 'done/Coverage/DBEJ0', path_in + 'done/Coverage/DCAN0', path_in + 'done/Coverage/DCOU0', path_in + 'done/Coverage/DERI0', path_in + 'done/Coverage/DHET1', path_in + 'done/Coverage/DINT1', path_in + 'done/Coverage/DLAS0', path_in + 'done/Coverage/DMAD1', path_in + 'done/Coverage/DMIR1', path_in + 'done/Coverage/DORO0', path_in + 'done/Coverage/DPRE0', path_in + 'done/Coverage/DREF0', path_in + 'done/Coverage/DSCH0', path_in + 'done/Coverage/DTEN0', path_in + 'done/Coverage/LESP0', path_in + 'done/Coverage/MSKO0', path_in + 'done/Coverage/DAND2', path_in + 'done/Coverage/DBER0', path_in + 'done/Coverage/DCAR0', path_in + 'done/Coverage/DCRI0', path_in + 'done/Coverage/DFAN0', path_in + 'done/Coverage/DHET2', path_in + 'done/Coverage/DINT2', path_in + 'done/Coverage/DLEP0', path_in + 'done/Coverage/DMAH0', path_in + 'done/Coverage/DMOC0', path_in + 'done/Coverage/DOVO0', path_in + 'done/Coverage/DPRO0', path_in + 'done/Coverage/DREM0', path_in + 'done/Coverage/DSCO0', path_in + 'done/Coverage/DTHE0', path_in + 'done/Coverage/LORU0', path_in + 'done/Coverage/MSKO1', path_in + 'done/Coverage/DANG0', path_in + 'done/Coverage/DBET0', path_in + 'done/Coverage/DCAT0', path_in + 'done/Coverage/DCUL0', path_in + 'done/Coverage/DFAN1', path_in + 'done/Coverage/DHET3', path_in + 'done/Coverage/DJER0', path_in + 'done/Coverage/DLEU0', path_in + 'done/Coverage/DMAK0', path_in + 'done/Coverage/DMON0', path_in + 'done/Coverage/DOVO1', path_in + 'done/Coverage/DPRO1', path_in + 'done/Coverage/DRIV0', path_in + 'done/Coverage/DSER0', path_in + 'done/Coverage/DTHI0', path_in + 'done/Coverage/LORU1', path_in + 'done/Coverage/MSMA0']
    outputs = [path_in + 'done/Retrieve_Genes/Retrieve_all_done.txt']
    options = {'cores': 10, 'memory': "20g", 'walltime': "12:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    source /home/sarahe/miniconda3/etc/profile.d/conda.sh
    source activate base

    cd {path_in}

    ls *trimmed.fasta > filelist.txt

    python /home/sarahe/GitHub/BSc/Python_Scripts/samples2genes.py > /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/done/Retrieve_Genes/outstats.csv

    touch /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/done/Retrieve_Genes/Retrieve_all_done.txt

    """.format(path_in = path_in)

    return (inputs, outputs, options, spec)

##########################################################################################################################
###############################################---- MAFT ----#############################################################
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
    options = {'cores': 1, 'memory': "20g", 'walltime': "1:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    source /home/sarahe/miniconda3/etc/profile.d/conda.sh
    source activate base

    python {path_python}post_mafft_renaming.py {path_in}{gene}_aligned.fasta > {path_out}{gene}_renamed_aligned.fasta

    touch {done}{gene}

    """.format(gene = gene, done = done, path_in = path_in, path_out=path_out, path_python=path_python)

    return (inputs, outputs, options, spec)


##########################################################################################################################
###############################################---- Trim ----#############################################################
##########################################################################################################################

def trim(path_python, path_in, path_out, done):
    """Trim the mafft files"""
    inputs = [path_in+'done/accD', path_in+'done/atpF', path_in+'done/clpP', path_in+'done/ndhB', path_in+'done/petA', path_in+'done/psaB', path_in+'done/psbA', path_in+'done/psbD', path_in+'done/rpl2', path_in+'done/rps12', path_in+'done/rps2', path_in+'done/rrn5', path_in+'done/trnI', path_in+'done/trnP', path_in+'done/ycf4', path_in+'done/atpA', path_in+'done/atpI', path_in+'done/matK', path_in+'done/ndhD', path_in+'done/petB', path_in+'done/psaC', path_in+'done/psbB', path_in+'done/rbcL', path_in+'done/rpoB', path_in+'done/rps16', path_in+'done/rps3', path_in+'done/trnA', path_in+'done/trnL', path_in+'done/trnS', path_in+'done/atpB', path_in+'done/ccsA', path_in+'done/ndhA', path_in+'done/ndhH', path_in+'done/petD', path_in+'done/psaI', path_in+'done/psbC', path_in+'done/rpl16', path_in+'done/rpoC1', path_in+'done/rps18', path_in+'done/rrn16', path_in+'done/trnG', path_in+'done/trnN', path_in+'done/trnV']
    outputs = [done, path_out+'concat1.fasta'] 
    options = {'cores': 1, 'memory': "20g", 'walltime': "5:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    source /home/sarahe/miniconda3/etc/profile.d/conda.sh
    source activate base

    cd {path_in}

    python {path_python}ConcatFasta.py --files *.fasta --dir . --outfile {path_out}concat1.fasta --part

    touch {done}

    """.format(path_python = path_python, path_in = path_in, path_out = path_out, done = done)

    return (inputs, outputs, options, spec)

##########################################################################################################################
###############################################---- IQtree ----###########################################################
##########################################################################################################################

def iqtree(inputs, done, path_out, path_in, path_part):
    """Runs IQtree from the concoctanated file from Trim"""
    inputs = [inputs]
    outputs = [done+'iqtree.txt'] 
    options = {'cores': 5, 'memory': "10g", 'walltime': "12:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

    spec = """
    source /home/sarahe/miniconda3/etc/profile.d/conda.sh
    source activate base

    cd {path_out}

    iqtree -s {path_in}concat1.fasta -T AUTO -B 1000 --redo -o LORU1

    touch {done}iqtree.txt

    """.format(inputs=inputs, done=done, path_out=path_out, path_in=path_in, path_part=path_part)

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

# Renaming mafft output files to be ready for Trim
for i in range(len(genes)):
    pth = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/04_mafft/"+genes[i]+'_aligned.fasta'
    if os.path.isfile(pth):
        gwf.target_from_template('post_mafft_'+str(i), post_mafft(gene = genes[i],
                                                                    path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/04_mafft/",
                                                                    path_out = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/05_post_mafft/",
                                                                    path_python = "/home/sarahe/GitHub/BSc/Python_Scripts/",
                                                                    done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/05_post_mafft/done/"))


# Running Trim
# gwf.target_from_template('trim', trim(path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/05_post_mafft/",
#                                         path_out = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/06_concatenate/",
#                                         path_python = "/home/sarahe/GitHub/BSc/Python_Scripts/",
#                                         done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/06_concatenate/done/concatenate_done"))

# After running Trim, you must look at the gwf logs Trim and create part.txt manually, and transfer the information about the alignment regions start and end to part.txt

# # Running IQtree
# gwf.target_from_template('IQtree', iqtree(inputs = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/06_concatenate/done/concatenate_done",
#                                         path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/06_concatenate/",
#                                         path_out = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/07_iqtree/",
#                                         path_part = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/05_post_mafft/",
#                                         done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/07_iqtree/done/"))