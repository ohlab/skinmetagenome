#script to combine 2 sam files, as is right now just combines two sam files, could have duplicate reads 
#to be used to combine the sam files created by the bacterial databases 
#takes three input arguments name 1: sam file 1 2: sam file 2 3: name of combined sam file 
#script that combines sam files created from multiple databases
#keeps track of which reads had hits in the multiple different databases 
#[Arch, Bact, Euky, Virl], T for read had a hit, F for the read didn't have a hit 

#this version combines all the sam files from all the different databases 

import sys, os


sam_out = sys.argv[1] #only writing the reads which aligned to something to avoid creating a massive file of no hits
#updating to take into account the newly merged bacterial database 
headers = sys.argv[2]
sam_arch = sys.argv[3]
sam_euky = sys.argv[4]
sam_virl = sys.argv[5]
sam_euky_2 = sys.argv[6]
sam_bac1 = sys.argv[7]
sam_bac2 = sys.argv[8]
if len(sys.argv) > 9:
	sam_bac3 = sys.argv[9]
if len(sys.argv) > 10:
	sam_bac4 = sys.argv[10]

##going through the arch sam file 

read_dic = {}

out = open(sam_out, 'w')
headers = open(headers, 'w')

arch = open(sam_arch, 'r')

print "Processing Archae"
for line in arch:
	if line.startswith('@'):
		headers.write(line)
	else:
		read = line.split("\t")[0].strip()
		if read not in read_dic:
			read_dic[read] = [False, False, False, False]
			
		value = int(line.split("\t")[1].strip())
		if value & 4 == 0:
			read_dic[read][0] = True 
			out.write(line)
arch.close()

print "Processing Bacteria 1"
bac1 = open(sam_bac1, 'r')

for line in bac1:
	if line.startswith('@'):
		headers.write(line)
	else:
		read = line.split("\t")[0].strip()
		if read not in read_dic:
			read_dic[read] = [False, False, False, False]
			
		value = int(line.split("\t")[1].strip())
		if value & 4 == 0:
			read_dic[read][1] = True 
			out.write(line)
bac1.close()

print "Processing Bacteria 2"			
bac2 = open(sam_bac2, 'r')

for line in bac2:
	if line.startswith('@'):
		headers.write(line)
	else:
		read = line.split("\t")[0].strip()
		if read not in read_dic:
			read_dic[read] = [False, False, False, False]
			
		value = int(line.split("\t")[1].strip())
		if value & 4 == 0:
			read_dic[read][1] = True
			out.write(line) 
bac2.close()

if len(sys.argv) >9:
	print "Processing Bacteria 3"			
	bac3 = open(sam_bac3, 'r')

	for line in bac3:
		if line.startswith('@'):
			headers.write(line)
		else:
			read = line.split("\t")[0].strip()
			if read not in read_dic:
				read_dic[read] = [False, False, False, False]
			
			value = int(line.split("\t")[1].strip())
			if value & 4 == 0:
				read_dic[read][1] = True
				out.write(line) 
	bac3.close()

if len(sys.argv) >10:
	print "Processing Bacteria 3"			
	bac4 = open(sam_bac4, 'r')

	for line in bac4:
		if line.startswith('@'):
			headers.write(line)
		else:
			read = line.split("\t")[0].strip()
			if read not in read_dic:
				read_dic[read] = [False, False, False, False]
			
			value = int(line.split("\t")[1].strip())
			if value & 4 == 0:
				read_dic[read][1] = True
				out.write(line) 
	bac4.close()
	
	

print "Processing Eukaryotes"			
euky = open(sam_euky, 'r')

for line in euky:
	if line.startswith('@'):
		headers.write(line)
	else:
		read = line.split("\t")[0].strip()
		if read not in read_dic:
			read_dic[read] = [False, False, False, False]
			
		value = int(line.split("\t")[1].strip())
		if value & 4 == 0:
			read_dic[read][2] = True 
			out.write(line)
euky.close()

print "Processing Virus"			
virl = open(sam_virl, 'r')

for line in virl:
	if line.startswith('@'):
		headers.write(line)
	else:
		read = line.split("\t")[0].strip()
		if read not in read_dic:
			read_dic[read] = [False, False, False, False]
			
		value = int(line.split("\t")[1].strip())
		if value & 4 == 0:
			read_dic[read][3] = True 
			out.write(line)
virl.close()

print "Processing Eukaryotes_2"			
euky = open(sam_euky_2, 'r')

for line in euky:
	if line.startswith('@'):
		headers.write(line)
	else:
		read = line.split("\t")[0].strip()
		if read not in read_dic:
			read_dic[read] = [False, False, False, False]
			
		value = int(line.split("\t")[1].strip())
		if value & 4 == 0:
			read_dic[read][2] = True 
			out.write(line)
euky.close()

headers.close()
out.close()
			
#counting up the different possibilities

A = 0
B = 0
E = 0
V = 0
AB = 0
AE = 0
AV = 0
BE = 0
BV = 0
EV = 0
ABE = 0
ABV = 0
AEV = 0
BEV = 0 
ABEV = 0  
none = 0 
				
for read in read_dic:

	arch = read_dic[read][0]
	bac = read_dic[read][1]
	euky = read_dic[read][2]
	viral = read_dic[read][3]
	
	if not arch and not bac and not euky and not viral:
		none += 1
	elif arch and not bac and not euky and not viral:
		A += 1
	elif not arch and bac and not euky and not viral:
		B += 1
	elif not arch and not bac and euky and not viral:
		E += 1
	elif not arch and not bac and not euky and viral:
		V += 1
	elif arch and bac and not euky and not viral:
		AB += 1
	elif arch and not bac and euky and not viral:
		AE += 1
	elif arch and not bac and not euky and viral:
		AV += 1
	elif not arch and bac and euky and not viral:
		BE += 1 
	elif not arch and bac and not euky and viral:
		BV += 1
	elif not arch and not bac and euky and viral:
		EV += 1
	elif arch and bac and euky and not viral:
		ABE += 1
	elif arch and bac and not euky and viral:
		ABV += 1
	elif arch and not bac and euky and viral:
		AEV += 1
	elif not arch and bac and euky and viral:
		#print read
		BEV += 1 
	else:
		ABEV += 1 

total = A + B + E + V + AB + AE + AV + BE + BV + EV + ABE + ABV + AEV + BEV + ABEV + none 

print "A:\t", A
print "B:\t", B
print "E:\t", E
print "V:\t", V
print "AB:\t", AB
print "AE:\t", AE
print "AV:\t", AV
print "BE:\t", BE
print "BV:\t", BV
print "EV:\t", EV
print "ABE:\t", ABE
print "ABV:\t", ABV
print "AEV:\t", AEV
print "BEV:\t", BEV
print "ABEV:\t", ABEV
print "Hit:\t", total - none 
print "none:\t", none
print "Total:\t", total 
		
print "Number reads in dictionary:\t", len(read_dic)

