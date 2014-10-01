#!/usr/bin/python
#extract_regions.py
#Author: Allyson Byrd
#Purpose: Write out genomes where all coordinates in the input file have been extracted in one file and reverse coordinates are masked in another
 
#Inputs: 1) alignments.txt: file of alignment coordinates with pieces of complete alignment listed at the end
			#File format: genome	genome_piece	start	end	piece_length
			#i.e. Propionibacterium_acnes_KPA171202.fa	gi|50841496|ref|NC_006085.1|	1	4	2560265
	#2) length_cutoff #pieces shorter than this cutoff will still not be written to extracted file 
	#3) input_genomes_dir #directory of the input genomes
	#4) output_genomes_dir #directory where the extracted and masked genomes will be written 
	
#Output: For each genomes in the coordinates file, a file where coordinates have been extracted (output_genomes_dir/genome.fa_extracted) and another file were reverse coordinates are masked (output_genomes_dir/genome.fa_masked)

#Usage: python extract_regions.py alignments.txt 25 input_genomes_dir output_genomes_dir

#Clarifications:
#coordinates in the list are extracted 
#coordinates not in the list are masked 
#if a region is marked as COMPLETE,it means we do NOT want to write it out here, but we DO want to write out regions which are NOT found as COMPLETE and are lacking any coordinates because it means they didn't align at all and we want to write out the whole piece

import sys, os

input_coords = sys.argv[1] #input coordinates of alignment 
length_cutoff = int(sys.argv[2]) #pieces shorter than this cutoff will still not be written to extracted file 

input_genomes_dir = sys.argv[3] #directory of the input genomes
output_genomes_dir = sys.argv[4] #directory where the extracted and masked genomes will be written

strains_dic = {}	# key: strain	#values: list of strain:piece [(strain1:pieceA), (strain1:pieceB),...]
coordinates_dic = {}	#key: strain:piece	#values:list of start:end [(start1:end1),(start2:end2),...]

#keeps track of which pieces are completely aligned 
complete_dic = {} #key: strain:piece 

#parsing input coords 
coords_FILE = open(input_coords, 'r')

for line in coords_FILE:
	line = line.strip().split('\t')
	
	if line[0] != 'COMPLETE':
		strain = line[0]
		piece = line[1]
		start = int(line[2])
		end = int(line[3])
		length = int(line[4])
	
		if not strains_dic.has_key(strain):
			strains_dic[strain] = []
	
		key = '%s::%s' % (strain, piece)
		if not coordinates_dic.has_key(key):
			coordinates_dic[key] = []
			strains_dic[strain].append(key)
		
		coordinateString = '%s::%s' % (start, end)
	
		coordinates_dic[key].append(coordinateString)
	else:
		key = line[1] #format is strain::piece
		complete_dic[key] = ''

coords_FILE.close()

#Processing each of the genomes(strains) that were referenced in the alignment file 
for strain in strains_dic:
	numNs = 0
	#used to keep track of the pieces for extraction later 
	complete_piece = {}  #key:header #value: sequence associated with that fasta header
	
	file_name = strain
	genome_FILE = open('%s/%s' % (input_genomes_dir, file_name), 'r')

	if True:
		masked_genome = open('%s/%s_masked' % (output_genomes_dir, file_name), 'w')
		extract_fragments = open('%s/%s_extracted' % (output_genomes_dir, file_name), 'w')

		total_count = 0 #used to keep track of how many bases are in the complete genome
		for line in genome_FILE:
			line = line.strip()
			if line.startswith('>'):
				count = 0
				header = line
				piece = line.split('|')[3]
				whole_piece = line.split(' ')[0].replace('>','')
				masked_genome.write(line+"\n")
				complete_piece[header] = ''
				key = '%s::%s' % (strain, whole_piece)
				key2 = '%s::%s' % (strain, whole_piece)

			else:
				complete_piece[header] += line
				if coordinates_dic.has_key(key2):
					newline = ""
					for letter in line:
						count += 1
						total_count += 1
						replaceWithN = True
						#need to make sure each letter falls within one of the coordinate ranges
						#if the letter doesn't fall within the range it needs to be masked with an N
						for coord in coordinates_dic[key2]:
							start = int(coord.split("::")[0])
							end = int(coord.split("::")[1])
							if start <= count <= end:
								replaceWithN = False
								#we already found in so now don't need to continuing looping through the loop
								break
						if replaceWithN:
							newline += "N"
							numNs += 1
						else:
							newline += letter
					#print newline
					masked_genome.write(newline+"\n")
				#this piece was completely something, it should be completely masked 
				elif complete_dic.has_key(key2):
					#print "Complete piece"
					newline = ""
					for letter in line:
						count += 1
						total_count += 1
						newline += "N"
						numNs += 1
					masked_genome.write(newline+"\n")
				#this piece was not something, so it should NOT be masked, so can just write the original line
				#this piece does not have a list of coordinates, nor is it in the complete list  
				else:
					#print "Complete piece not masked"
					assert not complete_dic.has_key(key)
					masked_genome.write(line+"\n")
					total_count += len(line)
			
			#line_count += 1	
			#if line_count % 100 == 0:
			#	print line_count
		#print len(complete_piece)
		print "Strain\t%s\t%s" % (strain, file_name)
		print "Total length\t%s" % total_count
		masked_per = float(numNs * 100.0/total_count)
		print "Masked\t%s\t%s" % (numNs, masked_per)
		unmasked = total_count - numNs
		print "Unmasked\t%s\t%s" % (unmasked, float(unmasked * 100.0/total_count) )
	
		#creating the file of extracted sequences 
		num_complete_pieces = 0
		for header in complete_piece:
			#print header
			#piece = header.split('|')[3]
			piece = header.split(' ')[0].replace('>','')
			sequence = complete_piece[header]
			key = '%s::%s' % (strain, piece)
			#print key
			header = header.replace(' ','_')
			#print "Number of fragments is %s" % len(coordinates_dic[key])
			#we are interested in portions of this sequence 
			if coordinates_dic.has_key(key):
				for coord in coordinates_dic[key]:
					start = int(coord.split("::")[0])
					end = int(coord.split("::")[1])
					length = end - start 
					#start doing minus 1 to compensate for the fact mummer coordinates are reported starting at 1, not 0
					#end is normal because the end coordinate is not contained in a string split
					extracted_piece = sequence[start -1:end]
					if length > length_cutoff:
						extract_fragments.write('%s_%s_%s\n' % (header,start,end) )
						extract_fragments.write('%s\n' % extracted_piece)
			#this piece does not have a list of coordinates, nor is it in the complete list, therefore we are interested in the whole sequence 
			elif not complete_dic.has_key(key):
				extracted_piece = sequence
				extract_fragments.write('%s_complete\n' % header )
				extract_fragments.write('%s\n' % extracted_piece)
			#this piece was marked as COMPLETE in the input file, thus we aren't interested in its contents 
			else:
				#print key
				pass

		genome_FILE.close()
		masked_genome.close()
		extract_fragments.close()

				
	

	
	
