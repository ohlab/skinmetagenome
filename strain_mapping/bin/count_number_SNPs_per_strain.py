#!/usr/bin/python
#count_number_SNPs_per_strain.py
#Author: Allyson Byrd
#Purpose: Creates a file that contains the number of SNPs per strain

#Input: 1) unique_SNPs.txt: file which lists every unique SNP per strain
			#File format: genome	genome_piece	SNP_position	SNP_value
			#i.e. Propionibacterium_acnes_J139.fa	taxid:679194|Propionibacterium_acnes_J139	538283	T
		#2) SNP_strain_counts.txt: file which lists the number of unique SNPs per strain
	
#Output: SNP_strain_counts.txt: file which lists the number of unique SNPs per strain

#Usage: python count_number_SNPs_per_strain.py unique_SNPs.txt SNP_strain_counts.txt

import sys

#Inputs
input = sys.argv[1] #file which lists every unique SNP per strain
output = sys.argv[2] #file which lists the number of unique SNPs per strain

strain_dic = {}

FILE = open(input, 'r')

for line in FILE:
	strain = line.split('\t')[1]
	if strain_dic.has_key(strain):
		strain_dic[strain] += 1
	else:
		strain_dic[strain] = 1

FILE.close()


FILE_OUT = open(output, 'w')

for strain in strain_dic:

	FILE_OUT.write('%s\t%s\n' % (strain, strain_dic[strain]) )


FILE_OUT.close()
