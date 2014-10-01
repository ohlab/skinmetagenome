#Purpose: analyze a sam file
#Returns: 1) text file of subject sequences and number of reads which aligned to each sequence
#only the top hits are recorded 
#2) updated fastq or fasta file, returns either reads which aligned or reads which didn't align
#which reads are returned is specified by the 5th command line argument
#1 for keeping aligned reads 2 for keeping unaligned reads
#this functionality replaces bowtie2's -un and -al options which are really time consuming

import sys

#sam file from bowtie
#sam files from other aligners may or may not work depending on the flags thrown
samFile = sys.argv[1]
#initial reads file used to generate the same file
readFile = sys.argv[2]
#name of the fasta file that will be created
newFasta = sys.argv[3]
#name of file containing the number of hits per each reference sequence 
HitsStatsFile = sys.argv[4]
#which database was being mapped against 
keepWhat = int(sys.argv[5])
#1 = keep aligned reads 2 = keep unaligned reads
	
if keepWhat == 1:
	keepAligned = True
else:
	keepAligned = False


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
		value = int(line.split("\t")[1].strip())
		#reference sequence the read aligned to  
		rname = line.split("\t")[2].strip()
		
		
		#8 flag unset means alignment is the primary alignment
		if value & 2**8 == 0:
			if not rname.startswith("*"):
				if not subjectHits_dic.has_key(rname):
					subjectHits_dic[rname] = 1
				else: 	
					subjectHits_dic[rname] += 1
			#4 bit flag not set, read had an alignment
			if value & 4 == 0:
				readHits_dic[read] = 1
				alignedReads += 1
			#4 bit flag set, read had no alignments 	
			else:
				unalignedReads += 1
		#nonprimary alignments 			
		else:
			#4 bit flag not set, read had an additional alignment
			if value & 4 == 0:
				readHits_dic[read] += 1
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

print "the length of hits dictionary is ", len(readHits_dic)

f = open(newFasta,'w')

readFound = False 

readsFile = open(readFile, 'r')

numReads = 0
linecount = 0

for line in readsFile:
	#need to adjust for fastq files, maybe make a way to automatically detect what comes after the @ in the fastq file
	if linecount == 0:
		#working with a fasta file
		if line.startswith(">"):
			counter = 2
			print "File is a fasta"
		#working with a fastq file
		else:
			counter = 4
			print "File is a fastq"
	 
	if linecount % counter == 0:
		numReads += 1
		line = line.strip()
		#cutting off the > or @
		readName = line[1:]
		readName = readName.split(" ")[0]
		#bowtie keeps these endings 
		#if line[-2:] == '/1' or line[-2:] == '/2':
		#	readName = line[1:-2]
		 
		if readHits_dic.has_key(readName):
			readFound = True 
			if keepAligned:
				f.write("%s\n" % line)
		else:
			readFound = False
			if not keepAligned:
				f.write("%s\n" % line)
			
	elif readFound == True:
		if keepAligned:
			f.write(line)

	else:
		if not keepAligned:
			f.write(line)
	linecount += 1
		
readsFile.close()
f.close()

#printing a file which contains the number of hits per each reference sequence which had > 1 hit 

f = open(HitsStatsFile, 'w')

for hit in subjectHits_dic:
	f.write("%s\t%s\n" % (hit, subjectHits_dic[hit]))

f.close()


print "###################Alignment Stats#########################"
print "Total Reads\t%s" % (numReads)
print "Total Reads Mapped\t%s\t%s" % (alignedReads, float(alignedReads*100)/numReads)
print "Mapped >1 Time\t%s\t%s" % (mappingMultipleTimes, float(mappingMultipleTimes*100)/numReads) 
print "Mapped 1 Time\t%s\t%s" % (mappingOnce, float(mappingOnce*100)/numReads)
print "Total Reads Unmapped\t%s\t%s" % (unalignedReads, float(unalignedReads*100)/numReads) 
