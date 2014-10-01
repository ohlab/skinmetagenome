#!/usr/bin/python
#parse_aligned_file.py
#Author: Allyson Byrd
#Purpose: Simplify file of redundant alignment coordinates to non-repetitive regions
#Input: 1) alignments_redundant.txt: file of repetitive alignment coordinates
			#File format: genome	genome_piece	start	end	piece_length
			#i.e. Propionibacterium_acnes_KPA171202.fa	gi|50841496|ref|NC_006085.1|	1	4	2560265
			#File needs to be sorted based on the piece and starting coordinate 
			#To sort file: sort -n -k 3 unsorted.txt > sorted.txt
		#2) alignments_parsed.txt: output file of non-repetitive alignment coordinates 
	
#Output: alignments_parsed.txt: file of non-repetitive alignment coordinates 

#Usage: python parse_aligned_file.py alignments_redundant.txt alignments_parsed.txt

#Able to remove those that are 1 off from each other and condense them together 
#for example
#strain	piece	start	end	length
#HL067PA1	gi|422571318|ref|NZ_GL383141.1|	6226	7733	15834
#HL067PA1	gi|422571318|ref|NZ_GL383141.1|	7734	15834	15834
#becomes
#HL067PA1	gi|422571318|ref|NZ_GL383141.1|	6226	15834	15834

import sys

#returns whether a number is contained within a given range 
def inside(testnum, beginRange, endRange):
   return beginRange <= testnum <= endRange

#Inputs 
input = sys.argv[1] #file of repetitive alignment coordinates, file needs to be sorted based on the piece and starting coordinate 
output = sys.argv[2] #file of non-repetitive alignment coordinates

#Dictionaries that will be populated throughout
strain_aligned_dic  = {} #key: strain:piece		value: list of start and stop values [(1,5),(10:23),(26:30),...]
piece_length_dic = {}	#key: piece		value: length of piece

#Processing the input files 
FILE = open(input,'r')

for line in FILE:
	line = line.strip().split('\t')
	strain = line[0]
	piece = line[1]
	start = int(line[2])
	end = int(line[3])
	
	#saving the lengths of each piece in a dictionary 
	#knowing the lengths is important when reversing the alignments
	if len(line) == 5:
		length = int(line[4])
		piece_length_dic[piece] = length

	#this works only when the input file is sorted by piece and start coordinate 
	key = '%s::%s' % (strain, piece)
	
	if not strain_aligned_dic.has_key(key):
		strain_aligned_dic[key] = [(start,end)]
		
	else:	
		#the previous start in the list comparing to
		start_to_compare = strain_aligned_dic[key][len(strain_aligned_dic[key])-1][0]
		#the previous end comparing to 
		end_to_compare = strain_aligned_dic[key][len(strain_aligned_dic[key])-1][1]
		
		#if the start is not within the range being looked at before it's created as a new range
		if not inside(start, start_to_compare, end_to_compare):
			#this catches those that are consecutive numbers, 144 144, 145 146
			#if the new start is just one greater than the previous end, the new end replaces the previous end 
			if (start - 1) == end_to_compare:
				#this updating the last set of ranges add to the list 
				strain_aligned_dic[key][len(strain_aligned_dic[key])-1] = (start_to_compare, end)
			else:
				strain_aligned_dic[key].append((start,end))
		#if the start is within the previous range, we need to see if the end is also contained
		else:
			#if the end is not contained, the new end replaces the other end 
			if not inside(end,start_to_compare,end_to_compare):
				strain_aligned_dic[key][len(strain_aligned_dic[key])-1] = (start_to_compare, end)
			#if the end is contained, the new set of numbers is perfectly contained within the previous set 
			else:
				pass

FILE.close()

#separating each of the ranges for every piece as a new line and writing out the results
OUT_FILE = open(output, 'w')
for key in strain_aligned_dic:
	strain = key.split('::')[0]
	piece = key.split('::')[1]
	length = piece_length_dic[piece]
	
	for range in strain_aligned_dic[key]:
		OUT_FILE.write('%s\t%s\t%s\t%s\t%s\n' % (strain, piece, range[0], range[1], length))
		
OUT_FILE.close()
		   