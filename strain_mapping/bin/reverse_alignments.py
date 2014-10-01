#!/usr/bin/python
#reverse_alignments.py
#Author: Allyson Byrd
#Purpose: Reverse list of alignment coordinates using the length of the segment
	#If input is list of aligned coordinates, list of unaligned coordinates will be returned and vice versa 
#Input: 1) aligned.txt: file of nonrepetitive alignment coordinates (likely generated from parse_aligned_file.py)
			#File format: genome	genome_piece	start	end	piece_length
			#i.e. Propionibacterium_acnes_KPA171202.fa	gi|50841496|ref|NC_006085.1|	1	4	2560265
			#File needs to be sorted based on the piece and starting coordinate 
			#To sort file: sort -n -k 3 unsorted.txt > sorted.txt
		#2) unaligned.txt: output file of reversed coordinates
	
#Output: unaligned.txt: file of reversed coordinates 

#Usage: python reverse_alignments.py aligned.txt unaligned.txt 
	#or python reverse_alignments.py unaligned.txt aligned.txt 

#Misc: If the input is aligned sequences, sequences with complete alignements are recorded in a list at the bottom of the outfile

import sys

#Inputs
input = sys.argv[1] #file of condensed coordinates
output = sys.argv[2] #file of reversed coordinated, with list of those completely covered at end 

#filled in from input file
pieces_dic = {} #key: strain:piece	#values: list of piece length, then start and stop coordinates [100,(1,5),(7,19),(24,28),...]

#inferred from information in pieces_dic
others_unaligned_dic = {}	#key: strain #values: list of (piece:unalign_start:unalign_end:unalign_length)

FILE = open(input, 'r')

for line in FILE:

	line = line.strip().split('\t')
	strain = line[0]
	piece = line[1]
	start = int(line[2])
	end = int(line[3])
	if len(line) == 5:
		piece_length = int(line[4])

	key = '%s:%s' % (strain, piece)	
	if not pieces_dic.has_key(key):
		pieces_dic[key] = [piece_length]
					
	pieces_dic[key].append((start,end))
	
	if not others_unaligned_dic.has_key(strain):
		others_unaligned_dic[strain] = []

FILE.close()

previous = ''

#stores only those segments where the whole thing was covered in the input file 
#if the input is aligned sequences, sequences in complete_pieces_dic are those sequences which were fully aligned 
complete_pieces_dic = {} #key: piece

#works most of the time, except when a reference region aligns to the non reference more than once
for piece in pieces_dic:
	x = 0
	for something in pieces_dic[piece]:
		if x == 0:
			length = something
		else:
			start = something[0]
			end = something[1]
			if start == 1 and end == length:
				#the whole segment aligned
				complete_pieces_dic[piece] = ''
			else:
				if x == 1:
					if start == 1:
						#means the beginning of the sequence aligned
						unalign_start = end + 1 
						#print pieces_dic[piece]
						#verifying if there is more than one range listed for a particular piece
						#the unalign end is 1 minus the next start sequence in the range
						if len(pieces_dic[piece]) != 2:
							unalign_end = pieces_dic[piece][x+1][0] -1
						else:
							unalign_end = length
					else:
						#means the beginning of the sequence failed to align
						unalign_start = 1
						unalign_end = start -1
						#need to add something here to add the other end
						
				unalign_length = unalign_end - unalign_start + 1
				info = '%s:%s:%s:%s' % (piece, unalign_start, unalign_end,unalign_length)
				if not info == previous:
					others_unaligned_dic[strain].append('%s:%s:%s:%s' % (piece, unalign_start, unalign_end,unalign_length))
				previous = info
				
				if x == len(pieces_dic[piece]) - 1:
				#means this is the last region of overlap
					if end == length:
						#means the end of the sequence aligned 
						pass
					else:
						unalign_start = end + 1
						unalign_end = length
				else:
					unalign_start = end + 1 
					unalign_end = pieces_dic[piece][x+1][0] -1
			
				#bt2 alignment 1-based offset into the forward reference strand where leftmost character of the alignment occurs
				unalign_length = unalign_end - unalign_start + 1
				info = '%s:%s:%s:%s' % (piece, unalign_start, unalign_end,unalign_length)
				if not info == previous:
					others_unaligned_dic[strain].append('%s:%s:%s:%s' % (piece, unalign_start, unalign_end,unalign_length))
				previous = info
		x += 1

#writing out the results, segments with complete alignment are listed at the end  
OUT = open(output, 'w')

for strain in others_unaligned_dic:
	#print strain
	for range in others_unaligned_dic[strain]:
		#print range
		split_range = range.split(':')
		strain = split_range[0]
		piece = split_range[1]
		start = split_range[2]
		end = split_range[3]
		length = split_range[4]
		OUT.write('%s\t%s\t%s\t%s\t%s\n' % (strain, piece, start, end, length) )

for piece in complete_pieces_dic:
	OUT.write("COMPLETE\t%s\n" % piece)

OUT.close()