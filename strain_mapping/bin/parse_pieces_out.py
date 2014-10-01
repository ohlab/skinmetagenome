#!/usr/bin/python
#parse_pieces_out.py
#Author: Allyson Byrd
#Purpose: Parse the counts output from the noncore pieces, so pieces from the same organism are summed together
	#Note: no normalization is performed  
#Inputs: 1) input.txt: file of counts for each piece, ideal if already normalized based on the length of the piece  
			#File format:piece 		norm_%		norm_counts		prenorm_counts
		#2) output.txt: summed counts for each organism 
			#File format:strain 		norm_%		norm_counts		prenorm_counts
		#3) pieces.fa_lengths_taxids: file which gives the length taxid and name for each piece
			#File format: piece		length		taxid		organism_name

#Output: output.txt: summed counts for each organism 
		#File format:strain 		norm_%		norm_counts		prenorm_counts
		
#Usage: python parse_pieces_out.py normalized_counts counts_per_species

import sys 

input = sys.argv[1] # file of counts for each piece, ideal if already normalized based on the length of the piece 
output = sys.argv[2] #file of counts for each species
taxid_lengths_file = sys.argv[3] #file which lists the taxid for each of the pieces 

#associating each piece to a strain using information in taxid_lengths_file
FILE = open(taxid_lengths_file, 'r')

piece_tax_name_dic = {}

for line in FILE:
	line = line.strip().split('\t')
	piece = line[0]
	taxid = line[2]
	name = line[3]
	piece_tax_name_dic[piece] = 'taxid:%s|%s' % (taxid, name)
	
FILE.close()

hits_per_dic = {}
hits_counts_dic = {}
hits_prenorm_dic = {}

#processing the counts files
FILE = open(input, 'r')

for line in FILE:
	if not line.startswith('Tax'):
		line = line.strip().split('\t')
		piece = line[0]
		norm_per = float(line[1])
		norm_count = float(line[2])
		prenorm_count = float(line[3])
		
		strain = piece_tax_name_dic[piece]
	
		if not hits_per_dic.has_key(strain):
			hits_per_dic[strain] = norm_per 
			hits_counts_dic[strain] = norm_count
			hits_prenorm_dic[strain] = prenorm_count
		else:
			hits_per_dic[strain] += norm_per
			hits_counts_dic[strain] += norm_count
			hits_prenorm_dic[strain] += prenorm_count

FILE.close()

#print ST_hits_dic

FILE_OUT = open(output, 'w')

FILE_OUT.write("Strain\tNormalized %\tNormalized Counts\tPrenormalized Counts\n")

#printing the strains in decreasing order of abundance

sorted_list = sorted(hits_per_dic, key=lambda x : hits_per_dic[x],reverse=True)

#print sorted_list

for strain in sorted_list:

	FILE_OUT.write("%s\t%s\t%s\t%s\n" % (strain, hits_per_dic[strain], hits_counts_dic[strain], hits_prenorm_dic[strain] ) )
	
FILE_OUT.close()