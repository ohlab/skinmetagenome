#!/usr/bin/python
#readCountNormalizer.py
#Author: Allyson Byrd
#Purpose: Normalize percentages by the total number of reads and the length of the sequence in the database
 
#Inputs: 1) counts.txt: input file of prenormalized read counts
			#File format: genome	count
			#i.e. taxid:679194|Propionibacterium_acnes_J139	3452
	#2) normalized_counts.txt: output file of normalized percentages per genome
		#File format: genome		normalized %		normalized counts		prenormalized counts
	#3) lengths.txt: file of the length of every genome in counts.txt
	
#Output: normalized_counts.txt: output file of normalized percentages per genome

#Usage: python readCountNormalizer.py counts.txt normalized_counts.txt lengths.txt

import sys, math

readCount = sys.argv[1] #input file of prenormalized read counts
normalizedCounts = sys.argv[2] #output file of post normalized read counts
#if a lengths file is provided the counts are normalized with regard to the length of the subject sequence
if len(sys.argv) > 3:
	lengthFile = sys.argv[3]
	norByLength =True
#if false (no length file is provided) just normalizes all the counts to 1, with no regard to length of the subject sequence
else:
	norByLength =False

file = open(lengthFile, 'r')

length_dic = {}

#saving the length of each sequence to a dictionary 
for line in file:
	line = line.strip().split('\t')
	length = float(line[1])
	name = line[0]

	length_dic[name] = length 

file.close()

#processing the file of read counts 
file = open(readCount, 'r')

counts_dic = {}
totalCounts = 0
raw_totalCounts = 0

otherFound = False

for line in file:
	line = line.strip().split('\t')
	count = float(line[1])
	name = line[0]
	raw_totalCounts += count
	
	#this works because other is always the last genome printed, this raw total Counts will always be the same 
	if name == 'Other':
		other = [count*100.0/raw_totalCounts, int(count)]
		otherFound = True 
		
	try:
		if not norByLength:
			counts_dic[name] = (float(count), float(count))
			totalCounts += float(count)
		else:
			counts_dic[name] = [float(count)/length_dic[name], int(count)]
			totalCounts += counts_dic[name][0] 
	except:
		pass
		#print name
	
	
file.close()

#print "The total number of counts ", totalCounts

sorted_list = sorted(counts_dic, key=lambda x : counts_dic[x][1],reverse=True)

file = open(normalizedCounts, 'w')

file.write("Taxon\tNormalized %\tNormalized Count\tPrenormalized count\n")

sum = 0

#taking into account the Other category that is not normalized by length
if otherFound: 
	non_other_percent = 100.0 - other[0]
else:
	non_other_percent = 100.0 

for name in sorted_list:

	normalized_per = float(counts_dic[name][0])/totalCounts * non_other_percent
	normalized_count = normalized_per * totalCounts/100
	file.write("%s\t%s\t%s\t%s\n" % (name, normalized_per, normalized_count, counts_dic[name][1]) )

	sum += float(counts_dic[name][0])/totalCounts

if otherFound:
	file.write("Other\t%s\t%s\t%s\n" % (other[0],other[1],other[1]) )
	
file.close()

#verifying that the sums should equal one 
#print 'Sum:\t', sum

