#normalize read counts by the total number of reads and the length of the sequence in the database
#depending on which length file can normalize whole species and individual contigs 
import sys, math

readCount = sys.argv[1] #input file of prenormalized read counts
normalizedCounts = sys.argv[2] #output file of post normalized read counts
#lengths of each of the sequences in the readCount File
#lengthFile = sys.argv[3]
if len(sys.argv) > 3:
	lengthFile = sys.argv[3]
else:
	lengthFile = 'all_lengths.txt'
#if false just normalizes all the counts to 1, with no regard to length of the subject sequence
if len(sys.argv) > 4:
# if len(sys.argv) > 3:
	coverages =False
else:
	coverages =True
	
norByLength =True 


file = open(lengthFile, 'r')

length_dic = {}

#identifies species that are already present in the database by going through the list of how long each of the species are in each of the databases
#in cases where multiple substrains are in a single database all those substrain lengths are summed together 
for line in file:
	lineWhole = line
	line = line.strip()
	line = line.split('\t')
	length = float(line[1])
	name = line[0]
	
	if length > 0:
		length_dic[name] = length 
	
	'''
	if length <= 0:
		length = 1
	name = line[0]
	if length_dic.has_key(name):
		#print "Already in database ", name
		#length_dic[name] += length
		pass
	else:
		length_dic[name] = length 
	'''

file.close()

file = open(readCount, 'r')

counts_dic = {}
totalCounts = 0

for line in file:
	line = line.strip()
	line = line.split('\t')
	count = float(line[1])
	name = line[0]
	if coverages:
		coverage1 = line[2]
		coverage2 = line[3]

	try:
		if not norByLength:
			if coverages:
				counts_dic[name] = (float(count), coverage1, coverage2, float(count))
			else:
				counts_dic[name] = (float(count), float(count))
			totalCounts += float(count)
		else:
			if coverages:
				counts_dic[name] = [float(count)/length_dic[name], coverage1, coverage2, int(count)]
			else:
				counts_dic[name] = [float(count)/length_dic[name], int(count)]
			totalCounts += counts_dic[name][0] 
	except:
		pass
		#print name
	
	
file.close()

print "The total number of counts ", totalCounts

sorted_list = sorted(counts_dic, key=lambda x : counts_dic[x][1],reverse=True)

file = open(normalizedCounts, 'w')

if coverages:
	file.write("Taxon\tNormalized %\tCount\t%Coverage\t%Coverage>5x\n")
else:
	file.write("Taxon\tNormalized %\tCount\n")

sum = 0

for name in sorted_list:
	if coverages:
		file.write("%s\t%s\t%s\t%s\t%s\n" % (name, float(counts_dic[name][0] * 100)/totalCounts, counts_dic[name][3], counts_dic[name][1], counts_dic[name][2]))
	else:
		file.write("%s\t%s\t%s\n" % (name, float(counts_dic[name][0] * 100)/totalCounts, counts_dic[name][1]))

	sum += float(counts_dic[name][0])/totalCounts

file.close()

#verifying that the sums should equal one 
print 'Sum:\t', sum

