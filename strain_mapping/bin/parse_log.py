#!/usr/bin/python
#parse_log.py
#Author: Allyson Byrd
#Purpose: Parse log from python comparative_genomics.py ref.fa ORG bin/ genome/ &> log.txt to more user friendly format
	#Can easily see the % alignment between the reference and other genomes
	#Useful to see if a particular genome is an outlier and needs to be removed 
#Inputs:  1) log.txt from python comparative_genomics.py ref.fa ORG bin/ genome/ &> log.txt
			#Format:
				# 2417915 of 2560265 nt (94.4%) reference positions are in aligned regions
				# 142350 of 2560265 nt (5.6%) reference positions are in unaligned regions
				# ...
				# mkdir: cannot create directory `genome/new_names//masked': File exists
				# Strain	Propionibacterium_acnes_KPA171202.fa	Propionibacterium_acnes_KPA171202.fa
				# Total length	2560265
				# Masked	223995	8.74889903975
				# Unmasked	2336270	91.2511009603
				# Reference genome: Propionibacterium_acnes_KPA171202.fa	 2560265 bp
				# Left out strains: none
				# perl /data/byrdal/strain_tracking_pipeline//bin_final/SNV_table_2.pl genome/new_names//Propionibacterium_acnes_KPA171202.fa genome/new_names/Propionibacterium_acnes_SK137.fa > Propionibacterium_acnes_KPA171202.fa_Propionibacterium_acnes_SK137.fa
				# ..
		 #2) output log_parsed.txt 
 
#Output: log_parsed.txt: More easily read version of the input log file 
		
#Usage: python parse_log.py log.txt log_parsed.txt

import sys

#inputs
input = sys.argv[1] #logfile
output = sys.argv[2]

FILE = open(input, 'r')

alignment_values = []
genomes = []

for line in FILE:
	
	if "in aligned regions" in line:
		percent = float(line.split('%')[0].split('(')[-1])
		alignment_values.append(percent)
	elif line.startswith("Masked"):
		masked_line = line
	elif line.startswith("Unmasked"):
		unmasked_line = line
	elif line.startswith("Reference"):
		reference_line = line
	elif line.startswith("Left out"):
		left_out_line = line
	elif line.startswith("perl"):
		line = line.split(' ')
		genome_name= line[3].split('/')[-1]
		genomes.append(genome_name)
	else:
		pass

FILE.close()

#associating genomes with their values 
genome_value_dic = {}
assert len(genomes) == len(alignment_values)

for i in xrange(len(genomes)):

	genome_value_dic[genomes[i]] = alignment_values[i]
	i += 1 


FILE_OUT = open(output, 'w')

FILE_OUT.write(reference_line)
FILE_OUT.write(masked_line)
FILE_OUT.write(unmasked_line)
FILE_OUT.write('\n')
FILE_OUT.write(left_out_line)
FILE_OUT.write('\n')

sorted_list = sorted(genome_value_dic, key=lambda x : genome_value_dic[x],reverse=True)

for key in sorted_list:
    FILE_OUT.write("%s\t%s\n" % (key, genome_value_dic[key]))

FILE_OUT.close()

		

