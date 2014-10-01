##using files created from sam tools and genome coverage bed
##first use sam tools to convert sam file to bam file
##sort the bam file
##use sorted bam file in genome coverage bed using default option
##genome coverage bed computes a histogram of coverage
##what percentage of the genome has 0 coverage what percentage has 5x coverage etc.
#automatically globs all the related coverages.txt files 
##Purpose: to assign what % of the genomes has any coverage and what % has coverage over the value specified as the 3rd command line argument
##Outputs readname	% coverage	% coverage over threshold

import sys, glob

countsFile = sys.argv[1] #file of counts for each of the reads
coverageInFile = sys.argv[2]
outFile = sys.argv[3] #file writing to 
coverage = int(sys.argv[4]) #coverage threshold for second number, i.e. 5 

sam_path_list = glob.glob('%s' % coverageInFile)

sub_dic = {}

for path in sam_path_list:
	
	#controls which files it analyzes based on which db was used 
	print path
	#parsing the genomecov files 
	file = open(path, 'r')
	for line in file:
		line = line.strip()
		line = line.split('\t')
		name = line[0]
		depth = int(line[1])
		fraction = float(line[4])
		
		if depth == 0:
			sub_dic[name] = [(1.0 - fraction) * 100.0, 0.0]
			
		elif depth >= coverage:
			try:
				sub_dic[name][1] += (fraction * 100.0)
			except:
				sub_dic[name] = [100, 0.0]
				sub_dic[name][1] += (fraction * 100.0)

	file.close()
			
countsFile = open(countsFile, 'r')
outFile = open(outFile, 'w')

for line in countsFile:
	line = line.strip()
	name = line.split('\t')[0]
	
	outFile.write("%s\t%s\t%s\n" % (line, sub_dic[name][0], sub_dic[name][1]))

countsFile.close()
outFile.close()
