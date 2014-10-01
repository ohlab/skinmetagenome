#!/usr/bin/python
#make_core_SNP_string.py
#Author: Allyson Byrd
#Purpose: To create a multifasta file from the SNPs identified in the core 
		#useful for creating a phylogenetic tree, fasta headers formatted for phyml

#Inputs: 1) reference.fa: name of the reference genome
	#2) SNPs.txt: reference_SNPs_detailed_sorted.txt from comparative_genomics.py SNP
		#Format: SNP	query_value	strain	query_seq	location_query	ref_value
		#i.e. 22	C	ATCC_11828	taxid|1091045|Propionibacterium_acnes_ATCC_11828	994433	A
	#3) out.fa: name of the output fasta file, where each genome is represented by a string of values at each SNP location
	
#Output: out.fa: name of the output fasta file, where each genome is represented by a string of values at each SNP location
	#headers are trimmed to less than 10 characters because that is the limit for phyml 

#Usage: python make_core_SNP_string.py reference.fa SNPs.txt out.fa

#Note: genome names need to be of the format genus_species_strain.fa
	#if the first 10 characters of a strain name are the same, the names should be manually edited in the 


import os, sys 

reference = sys.argv[1] #name of the reference genome 
SNPs = sys.argv[2]	#reference_SNPs_detailed_sorted.txt
out_fasta = sys.argv[3] #multifasta out file 

all_strains_dic = {}
ref_letter_dic = {}

short_name_dic = {}

found_dic = {} #counts number of SNPs found per strain

cline = 'cut -f 3 reference_SNPs_detailed_sorted.txt | sort | uniq > strains_list.txt'
os.system(cline)

FILE = open('strains_list.txt', 'r')

for line in FILE:
	strain = line.strip()
	all_strains_dic[strain] = []
	found_dic[strain] = 0
	
	#creating a short name for creating the fasta header by removing genus and species  
	genus = strain.split('_')[0]
	species = strain.split('_')[1]
	short_name = strain.replace(genus,'').replace(species,'').replace('.fa','').replace('_','')
	
	short_name_dic[strain] = short_name

FILE.close()

cline = 'rm strains_list.txt'
os.system(cline)

#adding the reference to the dictionaries 
all_strains_dic[reference] = []
found_dic[reference] = 0

genus = reference.split('_')[0]
species = reference.split('_')[1]
short_name = reference.replace(genus,'').replace(species,'').replace('.fa','').replace('_','')
short_name_dic[reference] = short_name


print "Fasta sequences will be created for %s strains" % len(all_strains_dic)

#Processing each of the SNPs
FILE = open(SNPs, 'r')

for line in FILE:
	line = line.strip().split('\t')
	location = line[0]
	query_value = line[1]
	strain_found = line[2]
	ref_value = line[5]
	
	found_dic[strain_found] += 1
	
	#going to automatically fill everything in with the reference value 
	#then change that value if its found in the file 
	#this ensures that every strain has a value for every SNP
	if not ref_letter_dic.has_key(location):
		for strain in all_strains_dic:
			all_strains_dic[strain].append(ref_value)
		ref_letter_dic[location] = ref_value
	
	all_strains_dic[strain_found][-1] = query_value 

FILE.close()

print "%s SNPs found in core regions" % (len(ref_letter_dic))


FILE_OUT = open(out_fasta, 'w')

for strain in all_strains_dic:

	if found_dic[strain] == 0 and strain != reference:
		print "%s has no SNP differences from the reference" % strain
	
	seq = ''.join(all_strains_dic[strain])
	
	assert len(ref_letter_dic) == len(seq)
	
	short_name = short_name_dic[strain]
	
	if len(short_name) > 10:
		shortened_name = short_name[:10]
		print "%s was shortended to %s" % (short_name, shortened_name)
		short_name = shortened_name
	
	FILE_OUT.write('>%s\n' % short_name)
	FILE_OUT.write('%s\n' % seq)

print "\tWarning: If several strains have the same shortened name, the names need to be manually changed in the fasta file" 

FILE_OUT.close()

	
		
	
