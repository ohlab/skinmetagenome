#!/usr/bin/python
#compare_SNPs.py
#Author: Allyson Byrd
#Purpose: Compare the SNPs found using 2 different references and return those SNPs which are identical and those SNPs which are different
 
#Inputs: 1) unique_SNPs_refA.txt: File which list SNPs unique to a strain when compared to reference A
			#File format: genome	genome_piece	SNP_position	SNP_value
			#i.e. Propionibacterium_acnes_J139.fa	taxid:679194|Propionibacterium_acnes_J139	538283	T
		#2) unique_SNPs_refB.txt: File which list SNPs unique to a strain when compared to reference B
			#File format: genome	genome_piece	SNP_position	SNP_value
			#i.e. Propionibacterium_acnes_J139.fa	taxid:679194|Propionibacterium_acnes_J139	538283	T
		#3) unique_SNPs_shared.txt: Lists the SNPs which were found in both files
			#File format: same as input files
		#4) unique_SNPs_diff.txt: Lists the SNPs which were not found in both input files 
			#File format: same as input files

#Output: 1) unique_SNPs_shared.txt: Lists the SNPs which were found in both files
		#2) unique_SNPs_diff.txt: Lists the SNPs which were found in both files

#Usage: python compare_SNPs.py unique_SNPs_refA.txt unique_SNPs_refB.txt unique_SNPs_shared.txt unique_SNPs_diff.txt

import sys

#Inputs
snp_file_1 = sys.argv[1]
snp_file_2 = sys.argv[2]

out_shared = sys.argv[3]
out_diff = sys.argv[4]

SNPs_dic = {} #key: line in input file #value: True if SNP is shared in both references False if SNP is unique to a reference
#True SNPs written to out_shared
#False SNPs written to out_diff

SNP1 = open(snp_file_1, 'r')

for line in SNP1:
	line = line.strip()
	#key is just the whole line 
	SNPs_dic[line] = False 

SNP1.close()

print "There are %s unique SNPs recorded in the first file" % len(SNPs_dic)

SNP2 = open(snp_file_2, 'r')
	
for line in SNP2:
	line = line.strip()
	
	#if the same line is found in both files, it means that SNP was found with both references
	if SNPs_dic.has_key(line):
		SNPs_dic[line] = True
	#if line was not found in the previous file, means SNP is unique to this reference
	else:
		SNPs_dic[line] = False
		
SNP2.close()

print "There are %s unique SNPs recorded in both files" % len(SNPs_dic)

#keep track of the number of SNPs shared or different between the different strains
shared_strain_dic = {}
diff_strain_dic = {}

num_shared = 0
num_diff = 0

SNPs_shared = open(out_shared, 'w')
SNPs_diff = open(out_diff, 'w')

for SNP in SNPs_dic:
	strain = SNP.split('\t')[0]
	
	#the SNP was found with both references 
	if SNPs_dic[SNP]:
		num_shared += 1
		SNPs_shared.write('%s\n' % SNP)
		if shared_strain_dic.has_key(strain):
			shared_strain_dic[strain] += 1
		else:
			 shared_strain_dic[strain] = 1
			 
	#the SNP was only found with one reference 
	else:
		num_diff += 1
		SNPs_diff.write('%s\n' % SNP)
		if diff_strain_dic.has_key(strain):
			diff_strain_dic[strain] += 1 
		else:
			diff_strain_dic[strain] = 1

SNPs_shared.close()
SNPs_diff.close()

print "The number of shared SNPs is %s" % num_shared
print "The number of different SNPs is %s" % num_diff

#outputs the number of SNPs shared and different for each strain 
print "Strain\tNum Shared\tNum Diff"
for strain in shared_strain_dic:
	if shared_strain_dic.has_key(strain):
		shared = shared_strain_dic[strain]
	else:
		shared = 0
	if diff_strain_dic.has_key(strain):
		diff = diff_strain_dic[strain]
	else:
		diff = 0
	
	print "%s\t%s\t%s" % (strain, shared, diff)
		
		
	
