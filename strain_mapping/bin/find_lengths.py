#!/usr/bin/python
#find_lengths.py
#Author: Allyson Byrd
#Purpose: Calculates the length of each sequence in a fasta file
	#ignores spacer sequences in the artificially merged files 
#Inputs: 1) input.fa #input fasta sequence
		#2) input.fa_lengths #gives the length of each sequence in the fasta file 
		#File format: fasta header		length
 
#Output: input.fa_lengths #gives the length of each sequence in the fasta file 

#Usage: python find_lengths.py input.fa input.fa_lengths 

import sys

inputFasta = sys.argv[1]
output = sys.argv[2]
file = open(inputFasta, 'r')

lengths_dic = {}

for line in file:
	
	line = line.strip()
	
	if line.startswith(">"):
		name = line[1:]
		lengths_dic[name] = 0
	#don't want to count these lines because I use them as spacer sequences in the artificially merged files
	elif line == 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN':
		pass
	else:
		lengths_dic[name] += len(line)
					
	
file.close()

#print "Number names in dic ", len(lengths_dic)

lengths_list = []

lengths = open(output, 'w')

for seq in lengths_dic:

	lengths.write("%s\t%s\n" % (seq, lengths_dic[seq]) )
	
	lengths_list.append(lengths_dic[seq])

lengths.close()

#print stats about the lengths 

print inputFasta
print "\tNumber sequences:\t%s" % len(lengths_list)
print "\tShortest:\t%s" % min(lengths_list)
print "\tLongest:\t%s" % max(lengths_list)
avg = sum(lengths_list)/len(lengths_list)
print "\tAverage:\t%s" % avg
print "\tTotal length:\t%s" % sum(lengths_list)
