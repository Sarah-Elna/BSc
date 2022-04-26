'''
------------------------------------------------------------------------------------------------------------------------
This workflow is used in a Workflow to estimate coverage. 
------------------------------------------------------------------------------------------------------------------------
This code is a variant of coverage.py by Oscar Wrisberg
Eddited by Sarah E.K. Kessel 
------------------------------------------------------------------------------------------------------------------------
'''

import os, argparse, subprocess
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("sample")
args = parser.parse_args()
sample = str(args.sample)

# depth required to KEEP (i.e. anything <trshld will be discarded)
trshld = 2

# Get all subdirectories in the current working directory. these are the loci recovered by hybpiper
loci = next(os.walk(sample))[1]
sequences = {}

for locus in loci: 
	#supercontig
	pth = sample+'/'+locus+'/'+sample+'/sequences/intron/'+locus+'_supercontig.fasta'	
	
	if os.path.isfile(pth):
		for record in SeqIO.parse(pth, "fasta"):
			record.id = record.id+'_'+locus
			#sequences.append(record)
			sequences[record.id] = record
			
with open('/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'.fasta', "w") as outfile:
 	SeqIO.write(list(sequences.values()), outfile, "fasta")


### /home/owrisberg/Coryphoideae/work_flow/02_trimmed/

print(sample+'.fasta generated')
	
# BWA index targets
cmd = 'bwa index /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'.fasta'
subprocess.call(cmd,shell=True)
print(sample+'.fasta indexed')

# BWA mem paired reads
cmd = 'bwa mem /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'.fasta /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/00_data/'+sample+'_clean-Read1.fastq /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/00_data/'+sample+'_clean-Read2.fastq | samtools view -b -o /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'.bam'
subprocess.call(cmd,shell=True)
print('paired reads mapped to '+sample+'.fasta')

# BWA mem unpaired reads
cmd = 'bwa mem /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'.fasta /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/00_data/'+sample+'_clean-Read12-single.fastq | samtools view -b -o /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'_up.bam'
subprocess.call(cmd,shell=True)
print('unpaired reads mapped to '+sample+'.fasta')

# merge BAM files
cmd = 'samtools merge -f /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'_all.bam /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'.bam /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'_up.bam'
subprocess.call(cmd,shell=True)
print('BAMs merged')

# sort and index BAM files
cmd = 'samtools sort /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'_all.bam -o /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'_all_sorted.bam'
subprocess.call(cmd,shell=True)
cmd = 'samtools index /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'_all_sorted.bam'
subprocess.call(cmd,shell=True)
print('BAM indexed and sorted')

# remove duplicates
cmd = 'picard MarkDuplicates I=/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'_all_sorted.bam O=/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'_all_sorted_deduplicated.bam M=/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'marked_dup_metrics.txt REMOVE_DUPLICATES=true'
subprocess.call(cmd,shell=True)
print('reads deduplicated for sample '+sample)

# calculate coverage
cmd = 'samtools depth /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'_all_sorted_deduplicated.bam > /home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/
print('coverage calculated for sample '+sample)

# define function to replace nth position of sequence with N
def n2N(sqnc, pstn):
	sqnc = list(sqnc)
	sqnc[int(pstn)-1] = "N"
	return "".join(sqnc)

# process coverage
with open('/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'.cov', "r") as covfile:
	for line in covfile:
		line = line.strip()
		LINE = line.split("\t")
		if int(LINE[2]) < trshld:
			sequences[LINE[0]].seq = n2N(sequences[LINE[0]].seq,LINE[1])

# remove unnecessary leading and trailing Ns
for nm in sequences.keys():
	sequences[nm].seq = sequences[nm].seq.strip("N")
	if isinstance(sequences[nm].seq, str):
		sequences[nm].seq = Seq(sequences[nm].seq)

print('coverage trimming completed, keeping only positions with coverage of '+str(trshld)+' or above')

# write outfile
with open('/home/sarahe/Dypsis_Chloroplast_Phylogeny/BSc/02_Coverage/'+sample+'_trimmed.fasta', "w") as outfile:
	SeqIO.write(list(sequences.values()), outfile, "fasta")
print('trimmed seqs written to '+sample+'_trimmed.fasta')

# remove unnecessary files
#cmd = "find ../coverage -type f ! -name '*.fasta' -delete"
#subprocess.call(cmd, shell=True)
#print('tidied up.')