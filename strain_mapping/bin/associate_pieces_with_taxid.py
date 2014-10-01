#!/usr/bin/python
#associate_pieces_with_taxid.py
#Author: Allyson Byrd
#Purpose: Associate pieces with a taxid and name, so pieces from the same organism can be easily associated with one another 
 
#Inputs: 1) input.txt: file which lists the pieces that need to be associated with a taxid, lengths can also be included in this file 
			#File format: genome_piece	any_additional_columns_of_information
			#Genome_piece should be named so that the gi number is always: gi|gi_number|ref|... because the gi number is how the piece is associated to a taxid
	#2) output.txt: outfile file which associates each of the input pieces with a taxid and organism name 
		#File format: genome_piece	any_additional_columns_of_information_in_input.txt		taxid		organism_name
	#3) gi_taxid.dmp: file from NCBI taxonomy that associates each gi number with a taxid num
	#4) names.dmp: file from NCBI taxonomy that associates each taxid num with a name 

#Output: output.txt: outfile file which associates each of the input pieces with a taxid and organism name 
		#File format: genome_piece	any_additional_columns_of_information_in_input.txt		taxid		organism_name

#Usage: python associated_pieces_with_taxid.py input.txt output.txt gi_taxid.dmp names.dmp

import sys

#Inputs
input = sys.argv[1] #input file which lists pieces
output = sys.argv[2] #where output will be written
gi_taxid = sys.argv[3]
taxid_names = sys.argv[4]

piece_line_dic = {}
gi_taxid_dic = {}
taxid_name_dic = {}

#extracting the gis from the pieces listed in input.txt

FILE = open(input, 'r')
for line in FILE:
	gi = line.split('\t')[0].split('|')[1]
	
	if not gi_taxid_dic.has_key(gi):
		gi_taxid_dic[gi] = ''

FILE.close()

#associating each of the gis with their taxonomy name using files from NCBI taxonomy 

gis_found = 0 #once all the gis are associated with a taxid we can stop reading in the file 
gi_FILE = open(gi_taxid, 'r')

for line in gi_FILE:
	line = line.strip().split('\t')
	gi = line[0]
	taxid = line[1]
	
	if gi_taxid_dic.has_key(gi):
		gi_taxid_dic[gi] = taxid
		gis_found += 1
		
		taxid_name_dic[taxid] = ''
		if gis_found == len(gi_taxid_dic): #all gis have been associated 
			break

gi_FILE.close()

#print gi_taxid_dic

#now that all the gis are associated with a taxid, all taxids need to be associated with a name

taxids_found = 0 #once all the taxids are associated with a name can stop reading in the file
taxid_FILE = open(taxid_names, 'r')

for line in taxid_FILE:
	line = line.split('|')
	taxid = line[0].strip()
	name = line[1].strip().replace(' ','_')
	label = line[3].strip()
	
	if taxid_name_dic.has_key(taxid) and taxid_name_dic[taxid] == '' and label  == 'scientific name':
		taxid_name_dic[taxid] = name
		taxids_found += 1 
		
		if taxids_found == len(taxid_name_dic): #all taxids have been associated with names 
			break
taxid_FILE.close()

#print taxid_name_dic

#using the information extracted from the NCBI files to associate each piece with a taxid and name
FILE_out = open(output, 'w')
FILE_in = open(input, 'r')

for line in FILE_in:
	gi = line.split('\t')[0].split('|')[1]
	taxid = gi_taxid_dic[gi]
	name = taxid_name_dic[taxid]
	line = line.strip()
	FILE_out.write('%s\t%s\t%s\n' % (line, taxid, name) )

FILE_out.close()
FILE_in.close()
