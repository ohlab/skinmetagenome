#!/usr/bin/python
#parse_bowtie2_output.py
#Author: Allyson Byrd
#Purpose: Parse a bowtie2 samfile and identify reads mapping to SNPs informative to a strain
 
#Inputs: 1) samFile: Bowtie2 alignment file from aligning reads against several genomes from different strains of the same species
	#2) outFile: indicates the number of unique SNPs found per strain, counts are normalized based on the total number of unique SNPs per strain 
	#3) unique_SNPs.tx: file which list SNPs unique to a strain
		#Format: genome	genome_piece	SNP_position	SNP_value
		#i.e. Propionibacterium_acnes_J139.fa	taxid:679194|Propionibacterium_acnes_J139	538283	T
	#4) SNP_strain_counts.txt: number of unique SNPs found per strain
		#Format: strain		number unique SNPs identified in the genome
	
#Output: outFile: indicates the number of unique SNPs found per strain, counts are normalized based on the total number of unique SNPs per strain 
	#Format: strain		number of unique SNPs found		normalized number SNPs found

#Usage: python parse_bowtie2_output.py samFile outFile unique_SNPs.txt SNP_strain_counts.txt

##Misc Notes##
#2 "unique" SNPs to different strains in the same read would prevent the read from mapping because it can't be a perfect match to different strains 
#To account for the fact that some genomes have a higher number of unique SNPs, the number of SNPs found is divided by the total number of unique SNPs which exists for that particular strain or ST.

#Explanation of the sam flags
#http://picard.sourceforge.net/explain-flags.html

import sys

sam_file = sys.argv[1]
outfile = sys.argv[2]
SNP_file = sys.argv[3]	
SNP_counts_file = sys.argv[4]

#holds letter value of SNP, when want to verify the read actually has the correct letter 
SNP_value_dic = {} #key: strain:location	value:letter 
strain_SNPs_dic = {} #key: strain	value: list of SNP locations

#returns whether a test number is between two other numbers 
def inside(testnum, beginRange, endRange):
   return beginRange <= testnum <= endRange

#returns the complement nucleotide    
def reverse(letter):
	if letter == 'A':
		return 'T'
	elif letter == 'T':
		return 'A'
	elif letter == 'G':
		return 'C'
	else:
		return 'G'


#Reading in SNPs unique to a strain
#going to always use the long name of the strain, not the short one, because the long one is how it shows up in the sam file 
num_snps = 0
SNP_FILE = open(SNP_file, 'r')
for line in SNP_FILE:
	#print line
	num_snps += 1
	line = line.strip().split('\t')
	strain = line[1]
	location = int(line[2])
	value = line[3]
	SNP_value_key = '%s__%s' % (strain, location) # use __ no _ as key because some names have _
	SNP_value_dic[SNP_value_key] = value
	
	if not strain_SNPs_dic.has_key(strain):
		strain_SNPs_dic[strain] = []
	strain_SNPs_dic[strain].append(location)
		

SNP_FILE.close

print "The number of input SNPs is %s" % num_snps

#sorting the dictionary so the coordinates are in ascending order 
#this is important when looping through all the SNPs so we know we can break when we pass a certain point
for strain in strain_SNPs_dic:
	strain_SNPs_dic[strain].sort()


#storing the total number of SNPs per each strain for normalization purposes 
strain_count_dic = {}

FILE = open(SNP_counts_file, 'r')
for line in FILE:
	line = line.strip().split('\t')
	strain = line[0]
	count = int(line[1])
	strain_count_dic[strain] = count
FILE.close()

strain_hits_dic = {} #key: strain	#value: # times unique SNP found 

#parsing the bowtie2 sam file 
SAM_FILE = open(sam_file, 'r')

for line in SAM_FILE:
	#bypassing the reference headers 
	if line.startswith('@'):
		pass
	else:
		line2 = line
		line = line.strip().split('\t')
		#name of read that aligned 
		read = line[0]
		#sum of thrown flags 
		value = int(line[1])
		#reference sequence the read aligned to  
		rname = line[2]
		
		#if value & 4 == 0:
		if rname != '*':
			#bowtie2 records start coordinate starting from 1 
			#nucmer also records from 1 so the values should coincide 
			start = int(line[3])
			read_seq = line[9]
			end = start + len(read_seq) - 1
			#now i need to loop go through all the dictionaries and see if we find something 
			#looping through strain_SNPs_dic
			strain = rname
			
			if strain_SNPs_dic.has_key(strain):
				SNP_list = strain_SNPs_dic[strain]
				for SNP in SNP_list:
					if inside(SNP, start, end):
						#need to verify that the correct letter is found 
						#key in SNP needs to be equal to the value of the letter in the read string, can be determined by looking at the coordinate of the SNP and the start of the read 				
						SNP_value_key = '%s__%s' % (strain, SNP)
						actual_letter = SNP_value_dic[SNP_value_key]
					
						read_SNP_index = SNP - start 
						read_letter = read_seq[read_SNP_index]
#checking that the value we have for the SNP matches the value we find in the read 
# 						if actual_letter == read_letter or reverse(actual_letter) == read_letter:				
# 							print "The letters DO match %s" % value					
# 						else:
# 							 print "The letters don't match %s" % value
						if not strain_hits_dic.has_key(strain):
 							strain_hits_dic[strain] = 1
 						else:
 							strain_hits_dic[strain] += 1

						break
					elif SNP > end:
						break
					else:
						pass

SAM_FILE.close()

#writing out the results, each of the SNP counts is normalized based on the total number of unique SNPs found for a particular strain
FILE_OUT = open(outfile, 'w')

FILE_OUT.write("

for strain in strain_hits_dic:
	norm_count = float(strain_hits_dic[strain])/strain_count_dic[strain]
	FILE_OUT.write("%s\t%s\t%s\n" % (strain, strain_hits_dic[strain], norm_count) )
	
FILE_OUT.close()
	
	
	