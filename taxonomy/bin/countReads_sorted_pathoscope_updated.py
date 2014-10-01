#counts the number of reads which aligned to each of the sequences in the database file
#does the same thing as updatedFastq but doesn't write a new reads file
#for use with sam files allowing multiple hits per read 
#counts all the hits per read not just the primary hit 
import sys

samFile = sys.argv[1]
#outfile which lists subject sequence	# reads which aligned 
HitsStatsFile = sys.argv[2]

alignedReads = 0
unalignedReads = 0

#stores names of reads which had alignments
readHits_dic = {}
#stores reference sequences which reads were aligned too 
subjectHits_dic = {}

file = open(samFile, 'r')

#parsing the sam file to see which reads aligned and didn't align 
for line in file:
	#bypassing the reference headers 
	if line.startswith('@'):
		pass
	else:
		#name of read that aligned 
		read = line.split("\t")[0].strip()
		#sum of thrown flags 
		try:
			value = int(line.split("\t")[1].strip())
		except:
			print line
		#reference sequence the read aligned to  
		rname = line.split("\t")[2].strip()
		
		if value & 4 == 0:
			if not rname.startswith("*"):
				if not subjectHits_dic.has_key(rname):
					subjectHits_dic[rname] = 1
				else: 	
					subjectHits_dic[rname] += 1
			if readHits_dic.has_key(read):
				readHits_dic[read] += 1
			else:
				readHits_dic[read] = 1
		else:
			unalignedReads += 1
		
file.close()
				
#print 'The lenght of the dictionary is ', len(readHits_dic)

#counting how many reads had multiple alignments 
mappingMultipleTimes = 0 
mappingOnce = 0
for read in readHits_dic:
	if readHits_dic[read] > 1:
		mappingMultipleTimes += 1
	else: 
		mappingOnce += 1

#printing a file which contains the number of hits per each reference sequence which had > 1 hit 

sorted_list = sorted(subjectHits_dic, key=lambda x : subjectHits_dic[x],reverse=True)

#print sorted_list

f = open(HitsStatsFile, 'w')

for something in sorted_list:
	f.write('%s\t%s\n' % (something, subjectHits_dic[something]))

#for hit in subjectHits_dic:
#	f.write("%s\t%s\n" % (hit, subjectHits_dic[hit]))

f.close()

print "The sam file is %s" % samFile
print "###################Alignment Stats#########################"
print "Total Reads Mapped\t%s" % (alignedReads)
print "Mapped >1 Time\t%s" % (mappingMultipleTimes) 
print "Mapped 1 Time\t%s" % (mappingOnce)
print "Total Reads Unmapped\t%s" % (unalignedReads) 
