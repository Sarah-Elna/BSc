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
    inputs = [done_path + 'Dypsis_acaulis_161017-1_S21_L001', done_path + 'Dypsis_lutescens_45152_S58_L001', done_path + 'Dypsis_acuminum_161213-3_S28_L001', done_path + 'Dypsis-madagascar-foxtail_S43_L001', done_path + 'Dypsis-aff-bejofo_S3_L001', done_path + 'Dypsis_madagascariensis_45090_S49_L001', done_path + 'Dypsis_albofarinosa_161213-10_S35_L001', done_path + 'Dypsis_mahia_161011-16_S11_L001', done_path + 'Dypsis_ambanjae_45083_S47_L001', done_path + 'Dypsis-makirae_S31_L001', done_path + 'Dypsis_ambilaensis_45078_S46_L001', done_path + 'Dypsis-malcomberi_S58_L001', done_path + 'Dypsis-ambositrae-SBL175-S45', done_path + 'Dypsis-mananjarensis_S1_L001', done_path + 'Dypsis-andapae-SBL237', done_path + 'Dypsis-mangorensis-SBL563', done_path + 'Dypsis-andilamenensis-SBL276', done_path + 'Dypsis-marojejyi_S67_L001', done_path + 'Dypsis_andrianatonga_161213-5_S30_L001', done_path + 'Dypsis-mcdonaldiana-SBL550', done_path + 'Dypsis-angusta-SBL238', done_path + 'Dypsis_metallica_161115-6_S20_L001', done_path + 'Dypsis_angustifolia_161213-12_S37_L001', done_path + 'Dypsis-minuta_S77_L001', done_path + 'Dypsis_ankirindro_45121_S55_L001', done_path + 'Dypsis-mirabilis-florencei_S64_L001', done_path + 'Dypsis_antanambensis_161213-6_S31_L001', done_path + 'Dypsis-mirabilis_S65_L001', done_path + 'Dypsis-aquatilis_S70_L001', done_path + 'Dypsis_mocquerysiana_161017-5_S4_L001', done_path + 'Dypsis_arenarum_45089_S48_L001', done_path + 'Dypsis-montana-SBL540', done_path + 'Dypsis_baronii_161011-13_S8_L001', done_path + 'Dypsis-moorei-SBL547', done_path + 'Dypsis_basilonga_161026-2_S6_L001', done_path + 'Dypsis-nauseosa_S57_L001', done_path + 'Dypsis_beentjei_161213-7_S32_L001', done_path + 'Dypsis_nodifera_45151_S77_L001', done_path + 'Dypsis_bejofo_45093_S50_L001', done_path + 'Dypsis-nossibensis_S13_L001', done_path + 'Dypsis_bernierana_161213-8_S33_L001', done_path + 'Dypsis_occidentalis_45080_S66_L001', done_path + 'Dypsis-betamponensis-BKL107', done_path + 'Dypsis_onilahensis_45015_S64_L001', done_path + 'Dypsis-betamponensis-SBL310', done_path + 'Dypsis-onilahensis-Isalo-form_S56_L001', done_path + 'Dypsis_betsimisarakae_45011_S39_L001', done_path + 'Dypsis_oreophila_45091_S69_L001', done_path + 'Dypsis_boiviniana_161213-11_S36_L001', done_path + 'Dypsis-oropedionis-SBL178', done_path + 'Dypsis_bonsai_45108_S52_L001', done_path + 'Dypsis-ovobontsira_S35_L001', done_path + 'Dypsis_bosseri_45019_S42_L001', done_path + 'Dypsis_ovojavavy_161115-5_S19_L001', done_path + 'Dypsis-brevicaulis_S36_L001', done_path + 'Dypsis-pachyramea-SBL316', done_path + 'Dypsis-brittiana-SBL277', done_path + 'Dypsis_paludosa_45084_S67_L001', done_path + 'Dypsis-cabadae-SBL194', done_path + 'Dypsis_perrieri_161115-1_S15_L001', done_path + 'Dypsis-canaliculata-SBL275', done_path + 'Dypsis-pervillei-SBL281', done_path + 'Dypsis-carlsmithii_S41_L001', done_path + 'Dypsis-pilulifera-SBL252', done_path + 'Dypsis-catatiana_S47_L001', done_path + 'Dypsis-pinnatifrons-SBL253', done_path + 'Dypsis-caudata-SBL539', done_path + 'Dypsis-plumosa_S32_L001', done_path + 'Dypsis_ceracea_45140_S56_L001', done_path + 'Dypsis-poivreana_S22_L001', done_path + 'Dypsis-commersoniana-SBL240', done_path + 'Dypsis_prestoniana_45129_S74_L001', done_path + 'Dypsis_concinna_45155_S60_L001', done_path + 'Dypsis_procera_161017-8_S23_L001', done_path + 'Dypsis_confusa_161115-3_S17_L001', done_path + 'Dypsis-procumbens-SBL183', done_path + 'Dypsis_cookei_161213-4_S29_L001', done_path + 'Dypsis_psammophila_45116_S73_L001', done_path + 'Dypsis-coriacea-SBL311', done_path + 'Dypsis-pulchella-SBL564', done_path + 'Dypsis-corniculata-SBL525', done_path + 'Dypsis-pumila_S75_L001', done_path + 'Dypsis_coursii_45109_S53_L001', done_path + 'Dypsis-pusilla_S17_L001', done_path + 'Dypsis_crinita_161213-9_S34_L001', done_path + 'Dypsis_pustulata_45114_S72_L001', done_path + 'Dypsis_culminis_161213-2_S27_L001', done_path + 'Dypsis_rabepierrei_161017-4_S14_L001', done_path + 'Dypsis-curtisii-SBL279', done_path + 'Dypsis-rakotonasoloi-SBL287', done_path + 'Dypsis-decaryi-SBL193-S47', done_path + 'Dypsis-ramentacea-SBL196', done_path + 'Dypsis_decipiens_45008_S38_L001', done_path + 'Dypsis_reflexa_140823-9_S2_L001', done_path + 'Dypsis_delicatula_45068_S44_L001', done_path + 'Dypsis-remotiflora-SBL286', done_path + 'Dypsis_digitata_161011-15_S10_L001', done_path + 'Dypsis-rivularis_S45_L001', done_path + 'Dypsis-dracaenoides_S19_L001', done_path + 'Dypsis_robusta_161011-14_S9_L001', done_path + 'Dypsis_dransfieldii_161017-6_S5_L001', done_path + 'Dypsis-rosea_S54_L001', done_path + 'Dypsis_elegans_161011-20_S13_L001', done_path + 'Dypsis-sahanofensis_S59_L001', done_path + 'Dypsis_eriostachys_45117_S54_L001', done_path + 'Dypsis_saintelucei_45139_S75_L001', done_path + 'Dypsis-faneva_S66_L001', done_path + 'Dypsis-saintelucei-SBL284', done_path + 'Dypsis-fanjana_S76_L001', done_path + 'Dypsis-sanctaemariae_S26_L001', done_path + 'Dypsis-fasciculata_S8_L001', done_path + 'Dypsis-sancta-SBL288', done_path + 'Dypsis-fibrosa_S71_L001', done_path + 'Dypsis_scandens_161026-1_S26_L001', done_path + 'Dypsis_forficifolia_161017-2_S22_L001', done_path + 'Dypsis-schatzii_S21_L001', done_path + 'Dypsis-gautieri_S61_L001', done_path + 'Dypsis_scottiana_45070_S65_L001', done_path + 'Dypsis-glabrescens_S18_L001', done_path + 'Dypsis_serpentina_45143_S76_L001', done_path + 'Dypsis-gronophyllum-SBL523', done_path + 'Dypsis-simianensis_S38_L001', done_path + 'Dypsis-henrici-SBL278', done_path + 'Dypsis_singularis_161011-18_S12_L001', done_path + 'Dypsis-heteromorpha-JD7822_S62_L001', done_path + 'Dypsis-sp-bef_S74_L001', done_path + 'Dypsis-heteromorpha-WLE111_S28_L001', done_path + 'Dypsis_spicata_45085_S68_L001', done_path + 'Dypsis-heterophylla-SBL179', done_path + 'Dypsis-sp-nov-aff-lutea_S46_L001', done_path + 'Dypsis-heterophylla-SBL179-repooled', done_path + 'Dypsis-sp-nov-sira_S7_L001', done_path + 'Dypsis_hiarakae_161017-3_S3_L001', done_path + 'Dypsis-subacaulis_S52_L001', done_path + 'Dypsis_hildebrandtii_45148_S57_L001', done_path + 'Dypsis-tanalensis-SBL180', done_path + 'Dypsis-hovomantsina_S9_L001', done_path + 'Dypsis-tenuissima_S29_L001', done_path + 'Dypsis-humbertii-SBL313', done_path + 'Dypsis_thermarum_161026-3_S7_L001', done_path + 'Dypsis-humblotiana-SBL562', done_path + 'Dypsis-thiryana_S60_L001', done_path + 'Dypsis_humilis_45094_S51_L001', done_path + 'Dypsis_tokoravina_45010_S63_L001', done_path + 'Dypsis-ifanadianae_S37_L001', done_path + 'Dypsis-trapezoidea_S68_L001', done_path + 'Dypsis_integra_45014_S41_L001', done_path + 'Dypsis_tsaravoasira_45098_S70_L001', done_path + 'Dypsis-intermedia-SBL314', done_path + 'Dypsis-turkii-SBL176', done_path + 'Dypsis-interrupta-SBL177', done_path + 'Dypsis-utilis_S27_L001', done_path + 'Dypsis-jeremiei-SBL245', done_path + 'Dypsis_viridis_45112_S71_L001', done_path + 'Dypsis_jumelleana_45157_S61_L001', done_path + 'Dypsis_vonitrandambo_161115-2_S16_L001', done_path + 'Dypsis-jurassic-park_S23_L001', done_path + 'Lemurophoenix-halleuxii-SBL195', done_path + 'Dypsis-laevis-SBL246', done_path + 'Lemurophoenix-sp-Nov-SBL471', done_path + 'Dypsis-lanceolata-SBL529', done_path + 'Loxococcus-rupicola-SBL234-S35', done_path + 'Dypsis-lantzeana-BKL103', done_path + 'Loxococcus-rupicola-SBL8-S7', done_path + 'Dypsis-lantzeana-BKL111', done_path + 'Marojejya_darianii_161115-4_S18_L001', done_path + 'Dypsis-lantzeana_S78_L001', done_path + 'Marojejya-darianii-BKL104', done_path + 'Dypsis-lanuginosa_S12_L001', done_path + 'Marojejya-darianii-BKL106', done_path + 'Dypsis_lastelliana_161017-17_S25_L001', done_path + 'Marojejya-darianii-SBL173-S39', done_path + 'Dypsis_leptocheilos_45022_S43_L001', done_path + 'Marojejya-insignis-BKL105', done_path + 'Dypsis-leucomalla-BKL108', done_path + 'Marojejya-insignis-BKL109', done_path + 'Dypsis-leucomalla-SBL247', done_path + 'Marojejya-insignis-SBL317', done_path + 'Dypsis_lilacina_45153_S59_L001', done_path + 'Marojejya-sp-BKL110', done_path + 'Dypsis-linearis_S2_L001', done_path + 'Masoala_kona_140823-8_S1_L001', done_path + 'Dypsis_lokohoensis_45012_S40_L001', done_path + 'Masoala-kona-SBL172', done_path + 'Dypsis-louvelii_S48_L001', done_path + 'Masoala_madagascariensis_161017-14_S24_L001', done_path + 'Dypsis-lutea-SBL282']
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

#def retrieve(path_in):
#    """Retrieve gene sequences from all the species and create an unaligned multifasta for each gene."""
#    inputs = [#! Coverage output]
#    outputs = ["#!/home/owrisberg/Coryphoideae/work_flow/04_coverage/done/Retrieve_Genes/Retrieve_all_done.txt"]
#    options = {'cores': 10, 'memory': "20g", 'walltime': "12:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

#    spec = """
#    source /home/sarahe/miniconda3/etc/profile.d/conda.sh
#    source activate base

#    cd {path_in}

#    ls *trimmed.fasta > filelist.txt

#    python #!/home/owrisberg/Coryphoideae/github_code/coryphoideae_species_tree/samples2genes.py > outstats.csv

#    touch #!/home/owrisberg/Coryphoideae/work_flow/04_coverage/done/Retrieve_Genes/Retrieve_all_done.txt

#    """.format(path_in = path_in)

#    return (inputs, outputs, options, spec)

##########################################################################################################################
###############################################---- MAFT ----#############################################################
##########################################################################################################################

#def mafft(gene, path_in, path_out, done):
#    """Aligning all the sequences for each gene."""
#    inputs = [#!"/home/owrisberg/Coryphoideae/work_flow/04_coverage/done/Retrieve_Genes/Retrieve_all_done.txt",path_in+gene+".FNA"]
#    outputs = [done,path_out+gene+"_aligned.fasta"] 
#    options = {'cores': 1, 'memory': "500g", 'walltime': "48:00:00", 'account':"Dypsis_Chloroplast_Phylogeny"}

#    spec = """

#    source activate #!mafft_env

#    cd {path_in}

#    mafft #!--globalpair #!--large #!--adjustdirectionaccurately #!--thread 1 {gene}.FNA > {path_out}{gene}_aligned.fasta
    
#    touch {done}

#    """.format(gene = gene, done = done, path_in = path_in, path_out=path_out)

#    return (inputs, outputs, options, spec)

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

# #get name list from hybpiper run
# gwf.target_from_template('name_list', get_namelist(done_path = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/done/",
#                                                     name_path = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/"))

# #get sequence length file necessary for statistical summary
# gwf.target_from_template('sequence_length', seq_lenghts(path = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/"))

# #get statistics, which can be used with gene_coverage_gg_plot.R to visualise results
# gwf.target_from_template('statistics', stats_summary(path = '/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/'))

#run intronerate
for i in [42, 66, 103]:
    gwf.target_from_template('Intronerate_'+sp[i], intronerate(species= sp[i],
                                                        path_in = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/",
                                                        done = "/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/01_HybPiper/done/Intronerate/"+sp[i]))

# #run Coverage to estimate the significance of the contigs found by hybpiper
# for i in range(0, len(sp)):
#     gwf.target_from_template('Coverage_'+sp[i], coverage(species = sp[i],
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
################################################################################################################################
#genes = ["gene00"]
#gt_values =["0.1","0.15","0.2","0.25","0.3","0.33","0.4","0.45","0.5","0.55","0.6","0.67","0.7","0.75","0.8","0.85","0.9","0.95"]

#gwf.target_from_template('Retrieve_genes', retrieve(path_in="/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/"))


#for i in range(len(genes)):
#    gwf.target_from_template('Mafft_'+genes[i], mafft(gene = genes[i],
#                                                        path_out= #!"/home/owrisberg/Coryphoideae/work_flow/06_alignment/",
#                                                        path_in = #!"/home/owrisberg/Coryphoideae/work_flow/05_blacklisting/",
#                                                        done = #!"/home/owrisberg/Coryphoideae/work_flow/06_alignment/done/"+genes[i]))