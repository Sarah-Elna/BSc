'''
------------------------------------------------------------------------------------------------------------------------
This workflow is used in a Workflow to estimate coverage. 
------------------------------------------------------------------------------------------------------------------------
This code is a variant of coverage.py by Wolf Eiserhardt
Eddited by Sarah E.K. Kessel 
------------------------------------------------------------------------------------------------------------------------
'''

import os, argparse, subprocess
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument("sample")
parser.add_argument("directory_in")
parser.add_argument("directory_out")
parser.add_argument("directory_wrk")
args = parser.parse_args()
sample = str(args.sample)
directory_in = str(args.directory_in)
directory_out = str(args.directory_out)
directory_wrk = str(args.directory_wrk)

# depth required to KEEP (i.e. anything <trshld will be discarded)
trshld = 2

# Go to working directory
cmd = 'cd '+directory_wrk
subprocess.call(cmd,shell=True)

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
			
with open(directory_out+sample+'.fasta', "w") as outfile:
 	SeqIO.write(list(sequences.values()), outfile, "fasta")

print(sample+'.fasta generated')
	
# BWA index targets
cmd = 'bwa index '+directory_out+sample+'.fasta'
subprocess.call(cmd,shell=True)
print(sample+'.fasta indexed')

# BWA mem paired reads and @HD tag
cmd = 'bwa mem '+directory_out+sample+'.fasta '+directory_in+sample+'_clean-Read1.fastq '+directory_in+sample+'_clean-Read2.fastq | samtools view -b -h -o '+directory_out+sample+'_no.bam | bam polishbam --HD "@HD VN:1.3 SO:coordinate" --in '+directory_out+sample+'_no.bam --out '+directory_out+sample+'.bam'
subprocess.call(cmd,shell=True)
print('paired reads mapped to '+sample+'.fasta')

# BWA mem unpaired reads and @HD tag
cmd = 'bwa mem '+directory_out+sample+'.fasta '+directory_in+sample+'_clean-Read12-single.fastq | samtools view -b -h -o '+directory_out+sample+'_no_up.bam | bam polishbam --HD "@HD VN:1.3 SO:coordinate" --in '+directory_out+sample+'_no_up.bam --out '+directory_out+sample+'_up.bam'
subprocess.call(cmd,shell=True)
print('unpaired reads mapped to '+sample+'.fasta')


# merge BAM files
cmd = 'samtools merge -f '+directory_out+sample+'_all.bam '+directory_out+sample+'.bam '+directory_out+sample+'_up.bam'
subprocess.call(cmd,shell=True)
print('BAMs merged')

# sort and index BAM files
cmd = 'samtools sort '+directory_out+sample+'_all.bam -o '+directory_out+sample+'_all_sorted.bam'
subprocess.call(cmd,shell=True)
cmd = 'samtools index '+directory_out+sample+'_all_sorted.bam'
subprocess.call(cmd,shell=True)
print('BAM indexed and sorted')

# remove duplicates
cmd = 'picard MarkDuplicates I='+directory_out+sample+'_all_sorted.bam O='+directory_out+sample+'_all_sorted_deduplicated.bam M='+directory_out+sample+'_marked_dup_metrics.txt REMOVE_DUPLICATES=true'
subprocess.call(cmd,shell=True)
print('reads deduplicated for sample '+sample)

# calculate coverage
cmd = 'samtools depth '+directory_out+sample+'_all_sorted_deduplicated.bam > '+directory_out+sample+'.cov'
print('coverage calculated for sample '+sample)

# define function to replace nth position of sequence with N
def n2N(sqnc, pstn):
	sqnc = list(sqnc)
	sqnc[int(pstn)-1] = "N"
	return "".join(sqnc)

# process coverage
with open(directory_out+sample+'.cov', "r") as covfile:
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
with open(directory_out+sample+'_trimmed.fasta', "w") as outfile:
	SeqIO.write(list(sequences.values()), outfile, "fasta")
print('trimmed seqs written to '+sample+'_trimmed.fasta')

# remove unnecessary files
#cmd = "find ../coverage -type f ! -name '*.fasta' -delete"
#subprocess.call(cmd, shell=True)
#print('tidied up.')