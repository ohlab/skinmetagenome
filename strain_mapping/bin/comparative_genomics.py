#!/usr/bin/python
#comparative_genomics.py
#Author: Allyson Byrd

#Purpose: Runs necessary scripts for comparative genomic analysis of multiple genomes

#Inputs: 1)reference.fa: name of the genome that will act as the reference
			#should be fully assembled and in one piece
		#2)task #CORE, SNP, or NONCORE, identifies what task you want to perform
			#CORE: Identify core region in reference genome
			#SNP: Identify SNPs in core regions
			#NONCORE: Identify noncore regions in the genomes 


#Usage: 
	#Identify core region in reference genome
		#python comparative_genomics.py reference.fa CORE script_dir genome_dir &> log.txt
			#log.txt can be processed by bin/parse_log.py for easier viewing 
		#Outputs:
			#reference_aligned_sorted_parsed.txt: coordinates of core regions in reference
			#reference_unaligned_sorted_parsed.txt: coordinates of noncore regions in reference 
			#genome_dir/core/reference.fa_extracted: reference genome where only core regions have been extracted
			#genome_dir/core/reference.fa_masked: reference genome where noncore regions have been masked
				
	#Identify core regions when leaving genomes out of the analysis 
		#python comparative_genomics.py reference.fa CORE script_dir genome_dir leave_out_1.fa leave_out_2.fa &> log.txt
		#Outputs: same as above
		
	#Identify SNPs in core regions
		#python comparative_genomics.py reference.fa SNP script_dir genome_dir
			#Uses a version of the reference genome where noncore regions have been masked (genome_dir/core/reference.fa_masked)
			#Uses the merged genomes because these are used in the strain tracking pipeline and we want the coordinates to match 
		#Outputs:
			#reference_SNPs_detailed_sorted.txt: detailed information for each SNP identified in the non-reference strains when compared to the reference strain core 
			#SNPs_summary.txt: number of strains which had a SNP at each of these reference sequence locations
			#unique_SNPs.txt: file which list SNPs unique to a strain
			#SNP_strain_counts.txt: file which lists the number of SNPs per strain
			#SNP_sequences.fa: fasta file, where each genome is represented by a string of values at each SNP location
				#this file can be used by phyml to make phylogenetic trees 
	
	#Identify noncore regions in the genomes 
		#python comparative_genomics.py reference.fa SNP script_dir genome_dir
			#Uses a version of the reference genome where noncore regions have been masked (genome_dir/core/reference.fa_masked)
		#Outputs:
			#other_aligned_parsed.txt: coordinates of core regions in non reference genomes
			#other_unaligned_parsed.txt: coordinates of noncore regions in non reference genomes 
			#all_unaligned_parsed.txt: coordinates of noncore regions in all genomes (non reference and reference)
			#genome_dir/noncore/genome.fa_extracted: for each genome a file where noncore regions have been extracted
			#genome_dir/noncore/genome.fa_masked: for each genome a file where core regions have been masked 
			#genome_dir/noncore/all_noncore_pieces.fa: all noncore pieces concatenated into one file 
			#genome_dir/noncore/all_noncore_pieces.fa_lengths_taxids: length and taxid of each of the pieces in the file 

#Notes: 
	#Requirements: 
		#mummer
		#fully assembled reference genome
	#Needs to be run in the output folder because various outputs will be written there
	#Many intermediate files are created and subsequently deleted 
			
	
#Required scripts that should be in script_dir 
	#associate_pieces_with_taxid.py: Associate pieces with a taxid and name, so pieces from the same organism can be easily associated with one another
	#count_number_SNPs_per_strain.py: Counts number of unique SNPs per strain
	#extract_regions.py: Writes out genomes where all coordinates in the input file have been extracted in one file and reverse coordinates are masked in another
	#find_lengths.py:Calculates the length of each sequence in a fasta file
	#make_core_SNP_string.py: To create a multifasta file from the SNPs identified in the core 
	#parse_aligned_file.py: Simplifies file of redundant alignment coordinates to non-repetitive regions
	#parse_SNPs.py: Parses SNP output from SNV_table_2.pl and identify SNPs in the core unique to a particular strain 
	#reverse_alignments.py: Reverses list of alignment coordinates using the length of the segment
	#SNV_table_2.pl: Runs nucmer comparisons 

#list of dictionaries populated throughout
	#ref_unaligned_dic = {}	#stores the coordinates in the reference strain that are unaligned
		#key: reference::reference::start::stop
		#value: list of every query strain which also had this exact region unaligned 
	#ref_SNP_dic = {}	#stores information for all SNPs found in relation to the reference 
		#key: location in the reference, i.e. 12312
		#value: list of every SNP in every sequence that is found at this location 
			#i.e. [genome::piece in genome::location in query::nucleotide value in query, C1::NC_018707.1::54625::A]
	#others_aligned_dic = {} #stores the coordinates of the aligned regions in the nonreference genomes 
		#key: genome
		#value: list of information about each aligned region for a strain 
			#i.e. [piece in genome::start::end::length of piece]
	#ref_SNP_value = {}	#stores the value of the SNP in the reference strain 
		#key: location in the reference, i.e. 12312
		#value: nucleotide value of SNP in the reference, i.e. A

import os,glob,sys

##input parameters 

reference = sys.argv[1] #name of sequence which will be used as the reference, this genome needs to be one contiguous piece 
task = sys.argv[2] #CORE, SNP, or NONCORE
script_dir = sys.argv[3] #directory where all the necessary accessory scripts are stored 
genome_dir = sys.argv[4] #directory where the genomes to be compared are located 

leave_out_strains = [] #remaining arguments are names of strains that should be excluded from the analysis 

i = 5
while i < len(sys.argv):
	leave_out_strains.append(sys.argv[i])
	i += 1

#validating a valid task was chosen 
task_options = ['CORE', 'SNP', 'NONCORE']
if task not in task_options:
	print "Please choose CORE, SNP, or NONCORE for the task"
	print "Exitting"
	sys.exit()

if task == 'CORE':
	reference_path = '%s/%s' % (genome_dir, reference)
else:
	reference_path = '%s/core/%s_masked' % (genome_dir, reference)

#finding the length of the reference sequence
reference_length = 0

FILE = open(reference_path , 'r')
FILE.readline()

for line in FILE:
	line = line.strip()
	reference_length += len(line)
	
FILE.close()

print "Reference genome: %s\t %s bp" % (reference,reference_length)

if task == 'CORE' or task == 'NONCORE':
	path_list = glob.glob('%s/*' % genome_dir)
else:
	assert task == 'SNP'
	path_list = glob.glob('%s/merged/*' % genome_dir)

#creating all the necessary dictionaries, see top for specifics on keys and values 
ref_unaligned_dic = {}	#stores the coordinates in the reference strain that are unaligned
ref_SNP_dic = {}	#stores the SNPs in the reference strain 
others_aligned_dic = {} #store the coordinates of the aligned regions in the reference strain
ref_SNP_value = {} #stores the value of the SNP in the reference strain 

if len(leave_out_strains) > 0:
	print "Left out strains: %s" % (','.join([str(x) for x in leave_out_strains]) )
else:
	print "Left out strains: none"

num_genomes = 0 #keeps track of the number of genomes that have been analyzed 

#looping through all the genomes we are interested in comparing 
for genome_path in path_list:

	genome = genome_path.split('/')[-1]
	strain = genome
	
	is_directory = os.path.isdir(genome_path) #indicates if the file path is actually indicates a directory, if so we want to skip it 
	
	#we don't want to analyze the reference or those genomes list in the leave_out list 
	if genome == reference or genome in leave_out_strains or is_directory:
		pass
	else:
		
		num_genomes += 1
		#running the mummer comparison
		#this first step successfully calculates SNPs relative to reference 
		
		output_name = '%s_%s' % (reference, genome)
		
		#running Sean's script, which automatically runs mummer and calculates SNPs and unaligned regions in the reference
		#to accurately calculate SNPs and unaligned regions, the reference sequence needs to be continuous 
		cline = 'perl %s/SNV_table_2.pl %s %s > %s' % (script_dir,reference_path, genome_path, output_name)
		print cline
		os.system(cline)
		
		#parsing the output from SNV_table_2.pl which is called reference_fasta_query_fasta
		output_1 = open(output_name, 'r')
		
		#started_unaligned = False
		previous_ref = -1
		
		#skipping the first line because it is header 
		output_1.readline()
		for line in output_1:
			
			line = line.strip().split('\t')
			location_ref = int(line[0])		#location of the SNP in the reference sequence 
			ref_value = line[1]				#nucleotide value of the SNP in the reference sequence 
			query_value = line[2]			#nucleotide value of the SNP in the query sequence 
			
			#Here we are parsing all the SNP data 
			#if the SNP is identified as "is_Ok"
			if line[6] == 'YES':
				location_query = int(line[3])	#location of the SNP in the query sequence 
				ref_piece = line[8]	#name of the reference sequence 
	
				query_seq = line[9]	#name of the query sequence 
				validSNP = True
		
				ref_key = location_ref	
				#using the location of the SNP in the reference to store the value of the SNP in the reference 
				ref_SNP_value[location_ref] = ref_value
		
				#now we going to store more detailed information about each of the SNPs and where they occur in the query strains 
				query_key = '%s::%s::%s::%s' % (strain, query_seq, location_query, query_value)
				if not ref_SNP_dic.has_key(ref_key):
					ref_SNP_dic[ref_key] = [query_key]
				else:
					ref_SNP_dic[ref_key].append(query_key)
			else:
				validSNP = False
			
			#This is important in the beginning so we mcan ask any part of the reference genome which isnt core 
			#recording information about all the unaligned regions is the reference 
			#because of the way the unaligned regions are indicated we have to decipher the start and stop ranges of the unaligned regions 
			if line[6] == 'no' and line[7] == 'unaligned':
				#we are just starting to analyze the unaligned regions 
				if previous_ref == -1:
					start = location_ref
					previous_ref = location_ref
				#the region of unalignment is continuing
				elif location_ref - previous_ref == 1:
					previous_ref = location_ref
				#we've identified the entire unaligned segment
				else:
					#start = start and stop previous_ref
					#print start
					#print previous_ref
					stop = previous_ref
					key = '%s::%s::%s::%s' % (reference, ref_piece, start, stop)
					if not ref_unaligned_dic.has_key(key):
						ref_unaligned_dic[key] = [strain]
					else:
						ref_unaligned_dic[key].append(strain)
					ref_unaligned_dic
					start = location_ref
					previous_ref = location_ref
					
		output_1.close()
		cline = 'rm %s' % output_name
		os.system(cline)
		
		if task == 'NONCORE':
			#running the opposite comparison to find the aligned regions in the nonreference  
			output_name = '%s_%s' % (genome, reference)
			cline = 'perl %s/SNV_table_2.pl %s %s > %s' % (script_dir, genome_path, reference_path, output_name)
			print cline 
			os.system(cline)
		
			#parsing the output of 
			#nucmer -p tmp sequence1 sequence2
			#show-coords -rTl tmp.delta > tmp.table 
			output_1 = open('tmp.table', 'r')
		
			others_aligned_dic[strain] = []
			pieces_dic = {} #don't think this is actually being used to store anything 
			#where pieces_dic[key][0] = length of the piece, everything else in the list is coordinates of alignment
			#example file 
			#[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[TAGS]
			#8	243139	144	243238	243132	243095	99.43	243139	2560265	gi|459815410|ref|NZ_AOGR01000001.1|	gi|50841496|ref|NC_006085.1|
			#because this is a reverse nucmer
			#the reference is the genomes we are looping through
			#the query is the genome we picked initially as the reference 
			for line in output_1:
				line = line.split('\t')
				if len(line) == 11:
					start = int(line[0]) #start of the alignment in the reference
					end = int(line[1])	#end of the alignment in the reference  
					piece_length = int(line[7])	
					reference_length = int(line[8])
				
					piece = line[9]	#name of the reference 
					if not pieces_dic.has_key(piece):
						pieces_dic[piece] = [piece_length]
					
					pieces_dic[piece].append((start,end))
				
					others_aligned_dic[strain].append('%s::%s::%s::%s' % (piece, start, end,piece_length) )
				
		
			#storing information for each genome piece in a dictionary, but this doesn't tell us which pieces are missing		
			output_1.close()
			
			cline = 'rm %s' % output_name
			os.system(cline)

#sys.exit()
#removing intermediate mummer files
cline = 'rm REP.* tmp.*'
os.system(cline)

##after having done all the pairwise alignments, now we want to write out the dictionaries 
if task == 'CORE':

	cline = 'mkdir %s/core' % genome_dir
	os.system(cline)
	
	REF_UNALIGN_OUT = open('reference_unaligned.txt', 'w')

	for key in ref_unaligned_dic:
		reference = key.split('::')[0]
		ref_piece = key.split('::')[1]
		start = key.split('::')[2]
		end = key.split('::')[3]
		REF_UNALIGN_OUT.write('%s\t%s\t%s\t%s\t%s\n' % (reference,ref_piece,start,end,reference_length) )

	REF_UNALIGN_OUT.close()

	os.system('sort -n -k 3 reference_unaligned.txt > reference_unaligned_sorted.txt')

	cline = "python %s/parse_aligned_file.py reference_unaligned_sorted.txt reference_unaligned_sorted_parsed.txt" % script_dir
	os.system(cline)

	cline = "python %s/reverse_alignments.py reference_unaligned_sorted_parsed.txt reference_aligned_sorted_parsed.txt" % script_dir 
	os.system(cline)

	cline = "python %s/extract_regions.py reference_aligned_sorted_parsed.txt 25 %s %s/core" % (script_dir, genome_dir, genome_dir)
	os.system(cline)
	
	#removing unnecessary intermediate files
	cline = "rm reference_unaligned.txt reference_unaligned_sorted.txt"
	os.system(cline)

elif task == 'SNP': 

	print "The total number of SNPs in relation to the reference is %s" % len(ref_SNP_dic)

	#parsing the information in ref_SNP_dic to print out the individual SNP information for each strain
	REF_SNP_DETAIL_OUT = open('reference_SNPs_detailed.txt', 'w')

	for SNP in sorted(ref_SNP_dic.iterkeys()):

		for something in ref_SNP_dic[SNP]: 
			SNP_split = something.split('::')
			strain = SNP_split[0]
			query_seq = SNP_split[1]
			location_query = SNP_split[2]
			query_value = SNP_split[3]
		
			ref_value = ref_SNP_value[SNP]
	
			REF_SNP_DETAIL_OUT.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (SNP, query_value, strain, query_seq, location_query, ref_value) )


	REF_SNP_DETAIL_OUT.close()

	os.system('sort -n -k 1 -k 2 reference_SNPs_detailed.txt > reference_SNPs_detailed_sorted.txt')

	reference_path = '%s/merged/%s' % (genome_dir, reference)
	cline = "python %s/parse_SNPs.py reference_SNPs_detailed_sorted.txt unique_SNPs.txt SNPs_summary.txt %s %s" % (script_dir, reference_path, num_genomes)
	os.system(cline)
	
	cline = 'python %s/count_number_SNPs_per_strain.py unique_SNPs.txt SNP_strain_counts.txt'  % script_dir
	os.system(cline)
	
	cline = 'python %s/make_core_SNP_string.py %s reference_SNPs_detailed_sorted.txt SNP_sequences.fa' % (script_dir, reference)
	os.system(cline)
	
	cline = "rm reference_SNPs_detailed.txt"
	os.system(cline)


elif task == 'NONCORE':

	cline = 'mkdir %s/noncore' % genome_dir
	os.system(cline)
	
	OTHER_ALIGN_OUT = open('other_aligned.txt', 'w')

	for strain in others_aligned_dic:

		for range in others_aligned_dic[strain]:
			split_range = range.split('::')
			piece = split_range[0]
			start = split_range[1]
			end = split_range[2]
			length = split_range[3]
			OTHER_ALIGN_OUT.write('%s\t%s\t%s\t%s\t%s\n' % (strain, piece, start, end, length) )

	OTHER_ALIGN_OUT.close()
	
	#these alignment coordinates are already in order so don't need to be sorted 
	cline = 'python %s/parse_aligned_file.py other_aligned.txt other_aligned_parsed.txt' % script_dir
	os.system(cline)

	cline = 'python %s/reverse_alignments.py other_aligned_parsed.txt other_unaligned_parsed.txt' %script_dir 
	os.system(cline)
	
	#combining with the unaligned coordinates of the reference (reference_unaligned_sorted_parsed.txt) so we can also write out the noncore region of the reference 
	cline = 'cat reference_unaligned_sorted_parsed.txt other_unaligned_parsed.txt > all_unaligned_parsed.txt'
	os.system(cline)

	cline = 'python %s/extract_regions.py all_unaligned_parsed.txt 25 %s %s/noncore' % (script_dir, genome_dir, genome_dir)
	os.system(cline)
	
	cline = 'cat %s/noncore/*_extracted > %s/noncore/all_noncore_pieces.fa' % (genome_dir, genome_dir)
	os.system(cline)
	
	cline = 'python %s/find_lengths.py %s/noncore/all_noncore_pieces.fa %s/noncore/all_noncore_pieces.fa_lengths' % (script_dir, genome_dir, genome_dir)
	os.system(cline)
	
	cline = 'python %s/associate_pieces_with_taxid.py %s/noncore/all_noncore_pieces.fa_lengths %s/noncore/all_noncore_pieces.fa_lengths_taxids %s/NCBI_taxonomy/gi_taxid_nucl.dmp %s/NCBI_taxonomy/names.dmp' % (script_dir, genome_dir, genome_dir, script_dir, script_dir)
	os.system(cline)

	cline = 'rm %s/noncore/all_noncore_pieces.fa_lengths other_aligned.txt' % genome_dir
	os.system(cline)

