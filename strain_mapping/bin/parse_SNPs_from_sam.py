#!/usr/bin/python
#parse_SNPs_from_sam.py
#Author: Allyson Byrd
#Purpose: Parse a bowtie2 samfile and identify reads mapping to SNPs informative to a strain
 
#Inputs: 1) samFile: Bowtie2 alignment file from aligning reads against several genomes from different strains of the same species
	#2) outFile1: indicates the number of unique SNPs found per strain, counts are normalized based on the total number of unique SNPs per strain 
	#3) outFile2: lists the specific SNPs which were found in each file 
	#4) unique_SNPs.tx: file which list SNPs unique to a strain
		#Format: genome	genome_piece	SNP_position	SNP_value
		#i.e. Propionibacterium_acnes_J139.fa	taxid:679194|Propionibacterium_acnes_J139	538283	T
	#5) SNP_strain_counts.txt: number of unique SNPs found per strain
		#Format: strain		number unique SNPs identified in the genome
	
#Output: 1) outFile1: indicates the number of unique SNPs found per strain, counts are normalized based on the total number of unique SNPs per strain 
			#Format: strain		total_number_SNPs_found		normalized_total_number_SNPs_found		norm_num_to_100		num_unique_SNPs_found	num_unique_SNPs_in_genome	%_of_those_unique_SNPs_found	
		#2) outFile2: lists the specific SNPs which were found in each file
			#Format: strain		location_in_strain		#_times_that_SNP_seen

#Usage: python parse_SNPs_from_sam.py samFile outFile1 outFile2 unique_SNPs.txt SNP_strain_counts.txt

##Misc Notes##
#2 "unique" SNPs to different strains in the same read would prevent the read from mapping because it can't be a perfect match to different strains 
#To account for the fact that some genomes have a higher number of unique SNPs, the number of SNPs found is divided by the total number of unique SNPs which exists for that particular strain or ST.
#this file of found SNPs is useful for comparing the SNPs found in 2 different samples

#Explanation of the sam flags
#http://picard.sourceforge.net/explain-flags.html

import sys

#Inputs
sam_file = sys.argv[1]
outfile = sys.argv[2]
outfile2 = sys.argv[3]
SNP_file = sys.argv[4]	
SNP_counts_file = sys.argv[5]

#holds letter value of SNP, when want to verify the read actually has the correct letter 
SNP_value_dic = {} #key: strain:location	value:letter 
strain_SNPs_dic = {} #key: strain	value: list of SNP locations'
found_SNPs_dic = {} #key: SNP identifier (strain:location) value: number of times the SNP was found in the file 
unique_SNPs_strain_dic = {} #key: strain value: number of unique SNPs found 


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

print "The number of SNPs in out_shared is %s" % num_snps

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
 							
 						SNP_id = '%s:%s' % (strain,SNP)
 						if not found_SNPs_dic.has_key(SNP_id):
 							found_SNPs_dic[SNP_id] = 1
 							if not unique_SNPs_strain_dic.has_key(strain):
 								unique_SNPs_strain_dic[strain] = 0
 							unique_SNPs_strain_dic[strain] += 1
 						else:
 							found_SNPs_dic[SNP_id] += 1
						break
					elif SNP > end:
						break
					else:
						pass

SAM_FILE.close()

#processing the SNPs found for individual strains 
totalForNorm = 0
norm_count_dic = {}
for strain in strain_hits_dic:
	norm_count_num_SNPs = float(strain_hits_dic[strain])/strain_count_dic[strain]
	totalForNorm += norm_count_num_SNPs
	norm_count_dic[strain] = norm_count_num_SNPs

sorted_list = sorted(norm_count_dic, key=lambda x : norm_count_dic[x],reverse=True)

FILE_OUT = open(outfile, 'w')

FILE_OUT.write("Strain\ttotal_number_SNPs_found\tnormalized_total_number_SNPs_found\tnorm_num_to_100\tnum_unique_SNPs_found\tnum_unique_SNPs_in_genome\t%_of_those_unique_SNPs_found\n")	

for strain in sorted_list:
	norm_count = float(strain_hits_dic[strain])/strain_count_dic[strain]
	norm_to_100 = norm_count*100.0/totalForNorm
	percent_found = float(unique_SNPs_strain_dic[strain])/strain_count_dic[strain]
	FILE_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (strain, strain_hits_dic[strain], norm_count,norm_to_100,unique_SNPs_strain_dic[strain],strain_count_dic[strain],percent_found) )
FILE_OUT.close()

FILE_OUT2 = open(outfile2, 'w')

FILE_OUT2.write("strain\tlocation_in_strain\t#_times_that_SNP_seen\n")

for SNP_id in found_SNPs_dic:
	strain = SNP_id.split(':')[1]
	location = SNP_id.split(':')[2]
	FILE_OUT2.write('taxid:%s\t%s\t%s\n' % (strain, location, found_SNPs_dic[SNP_id]) )

FILE_OUT2.close()

	
	