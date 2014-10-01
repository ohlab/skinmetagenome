#to see what the overlap of reads is between kingdoms in the updated sam file 
#[Arch, Bact, Euky, Virl], T for read had a hit, F for the read didn't have a hit 

import sys 

sam = sys.argv[1]
lineages = sys.argv[2]

lineages = {}

FILE = open(lineages, 'r')

for line in FILE:
	line = line.strip()
	list = line.split('\t') 
	if list[8] != 'NULL':
		key = list[8]
	else:
		key = list[7]
	
	lineages[key] = int(list[0])
	#print key
	#print list
	#sys.exit()

FILE.close()

#[superkingdom 0, kingdom 1, phylum 2, class 3, order 4, family 5, genus 6, species 7, subspecies 8]

read_dic = {}

FILE = open(sam, 'r')

for line in FILE:
	if line.startswith('@'):
		pass
	else:
		line = line.split('\t')
		read = line[0].strip()
		taxonID = line[2].split('|')[0].replace('taxid:','')
		
		try:
			kingdom = lineages[taxonID]
		except:
			print taxonID
			kingdom = -1
			
		if read not in read_dic:
			read_dic[read] = [False, False, False, False]
			
		#[Arch 0, Bact 1, Euky 2, Virl 3]	
		value = int(line[1])
		if value & 4 == 0:
			#this is some test to see what kingdom or superkingdom it is 
			#Bacteria
			if kingdom == 2:
				read_dic[read][1] = True
			#eukary
			elif kingdom == 2759:
				read_dic[read][2] = True
			#virus
			elif kingdom == 10239:
				read_dic[read][3] = True
			#Archea
			elif kingdom == 2157:
			#else:
			#	assert lineages[taxonID] == 2157
				read_dic[read][0] = True
			else:
				pass
				#print taxonID
FILE.close()

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
