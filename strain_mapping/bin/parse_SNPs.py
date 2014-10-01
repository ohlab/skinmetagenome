#!/usr/bin/python
#parse_SNPs.py
#Author: Allyson Byrd
#Purpose: Parse SNP output from SNV_table_2.pl and identify SNPs in the core unique to a particular strain 
 
#Inputs: 1) reference_SNPs_detailed_sorted.txt: file of SNPs identified in the core region of each genome versus the reference 
			#File format: SNP_position_in_reference		SNP_value_in_nonref_genome		nonref_genome	nonref_genome_piece		SNP_position_in_nonref_genome		SNP_values_reference
			#i.e. 22	C	Propionibacterium_acnes_J139.fa	taxid:679194|Propionibacterium_acnes_J139	98312	A
		#2) unique_SNPs.txt: name of output file that will contain the SNPs unique to a strain
		#3) SNPs_summary.txt : name of output file that will contain the number of strains which had a SNP at each particular coordinate 
		#4) path/to/reference.fa: path to reference genome 
		#5) num_genomes: number of genomes compared

#Output: 1) unique_SNPs.txt: File which list SNPs unique to a strain
			#File format: genome	genome_piece	SNP_position	SNP_value
			#i.e. Propionibacterium_acnes_J139.fa	taxid:679194|Propionibacterium_acnes_J139	538283	T
		#2) SNPs_summary.txt: File which lists number of strains which had a SNP at each particular coordinate 
			#File format: SNP:Value	#Strains
			#i.e. 438376::T		1
			
#Usage: python parse_SNPs.py reference_SNPs_detailed_sorted.txt unique_SNPs.txt SNPs_summary.txt path_to_reference.fa num_genomes

import sys

#Inputs
input = sys.argv[1] #file of SNPs identified in the core region of each genome versus the reference
output1 = sys.argv[2] # name of output file that will contain the SNPs unique to a strain
output2 = sys.argv[3] #name of output file that will contain the number of strains which had a SNP at each particular coordinate 
reference_path = sys.argv[4] #location of the reference genome 
num_strains = int(sys.argv[5]) #number of genomes compared

reference = reference_path.split('/')[-1].replace("_masked", "")
#opening the reference file to find the name of the piece by reading the first line
FILE = open(reference_path, 'r')

reference_piece = FILE.readline().strip().replace('>','')

#Purpose: count number of instances of each thing in a list
#Input: list
#Output: dictionary of counts 
def count_things(list):
	dic = {}
	for something in list:
		if dic.has_key(something):
			dic[something] += 1
		else:
			dic[something] = 1
	return dic

SNP_strain_dic = {}

SN_detail_dic = {} 

ref_values_dic = {}

#parsing the valid SNPs found from SNV_table_2.pl
FILE = open(input, 'r')

for line in FILE:

	line = line.strip().split('\t')
	location = line[0]
	letter = line[1]
	strain = line[2]
	piece = line[3]
	other_location = line[4]
	ref_letter = line[5]
	ref_values_dic[location] = ref_letter

	key = '%s::%s::%s' % (location, letter, strain)
	SN_detail_dic[key] = '%s\t%s' % (piece, other_location)
	
	SNP_key = '%s::%s' % (location, letter)
	
	if not SNP_strain_dic.has_key(SNP_key):
		SNP_strain_dic[SNP_key] = []

	SNP_strain_dic[SNP_key].append(strain)
	
FILE.close()

#in this file the SNPs are completely stripped of their identity related to the reference 
unique_SNPs = open(output1, 'w')

out_FILE = open(output2,'w')
out_FILE.write('SNP::Value\t#Strains\n')

for SNP in SNP_strain_dic:

	number_strains = len(SNP_strain_dic[SNP])

	#writing out useful statistics for every SNP
	out_FILE.write('%s\t%s\n' % (SNP, number_strains) )

	
	#means the SNP is unique to that strain and can potentially be used to identify it
	if number_strains == 1:
		strain = SNP_strain_dic[SNP][0]
		#print "this SNP is unique to a particular strain"
		unique_SNP_key = '%s::%s' % (SNP, strain)
		detail = SN_detail_dic[unique_SNP_key]
		letter = SNP.split('::')[1]
		#will write a file of SNPs unique to each strain Strain	Piece	Location	Letter
		unique_SNPs.write('%s\t%s\t%s\n' % (strain, detail,letter))
	#this means this SNP is unique to the reference 
	if number_strains == num_strains -1:
		strain = reference 
		piece = reference_piece
		ref_location = SNP.split('::')[0]
		letter = ref_values_dic[location] 
		unique_SNPs.write('%s\t%s\t%s\t%s\n' % (strain, piece, ref_location, letter) )

out_FILE.close()
unique_SNPs.close()



