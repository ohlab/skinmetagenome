#!/usr/bin/python
#inferNumberReads_fromPatho_tsv.py
#Author: Allyson Byrd
#Purpose: Infer number of reads from the patho .tsv file, genomes which fall below a certain percentage are grouped into "Other"
 
#Inputs: 1) input.tsv: tsv output file from Pathoscope
	#2) outFile.txt: output file of counts for each genome calculated from tsv file 
	#3) cutoff: i.e. 1 or 5 genomes which fall below this percentage are grouped into a category called other 

#Output: outFile.txt: output file of counts for each genome calculated from tsv file 
	#Format: genome		counts

#Usage: python inferNumberReads_fromPatho_tsv.py input.tsv outFile.txt cutoff

import sys 

#inputs
patho_tsv = sys.argv[1] #tsv file from Pathoscope
output = sys.argv[2] #output file of counts calculated from tsv file
#example: 1 or 5 
if len(sys.argv) > 3: 
	cutoff = int(sys.argv[3])*.01 #genomes which fall below this percentage i.e. (1% or 5%) are grouped into a category called other
else:
	cutoff = 0

total_above_cutoff = 0
other = 0

subject_dic = {}

#processing the tsv file
FILE = open(patho_tsv ,'r')

for line in FILE:

	line = line.split('\t')
	if line[0] == 'Total Number of Aligned Reads:':
		all_reads = int(line[1])
	elif line[0] == 'Genome':
		pass
	else:
		genome = line[0]
		final_guess = float(line[1])
				
		if final_guess > cutoff:
			subject_dic[genome] = final_guess
			total_above_cutoff += final_guess
		#genomes with a File Guess below the cutoff are grouped as other 
		else:
			other += final_guess
		

FILE.close()

total = other + total_above_cutoff 

FILE = open(output, 'w')

#adjusting the percentage of each genome by the sum of 'Best Hit'
#This is necessary because 'Best Hit' does not add to one
#'Best Hit' does not add up to one because reads which align to many things are automatically discarded from the Pathoscope output 

for genome in subject_dic:

	adjusted_percent = subject_dic[genome]/total
	num_reads = round(adjusted_percent * all_reads)
	FILE.write('%s\t%s\n' % (genome, num_reads) )

if other > 0:
	FILE.write('Other\t%s\n' %(round(other/total * all_reads)))
	
FILE.close()


