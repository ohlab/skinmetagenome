#!/usr/bin/python
#run_strain_stracking.py
#Authors: Allyson Byrd
#Purpose: Parses config file and generates strain tracking pipeline shell script 
#Input: config file
#Output: shell script 
#Usage: python run_strain_tracking.py config

import sys
import os

if len(sys.argv) > 1:
	config = sys.argv[1]
else:
	print "Error: No config file"
	print "Usage: python run_clinical_pathoscope.py config"
	sys.exit(1)

fastq_file = ''

bowtie2_bin = ''
script_bin = ''

output_dir = ''

whole_genome_bt2 = ''
whole_genome_len = ''

noncore_bt2 = ''
noncore_len = ''

SNP_file = ''
SNP_count_file = ''

#additional bowtie2 parameters 
parameter = '--very-sensitive --score-min L,-0.6,0.006'
bowtie2_k = 10 #will be replaced by number_genomes in config file
num_threads = 8

#Parsing the inputs in the config file
FILE_IN = open(config, 'r')
#fatal errors prevent the creation of the shell script 
fatal_errors = 0
for line in FILE_IN:
	if line.startswith('#'):
		pass
	if '=' in line:
		line = line.split('=')
		variable = line[0].strip()
		value = line[1].strip()
		
		#error checking of all user variables 
		if variable == '*filename':
			if value == '':
				print 'Fatal Warning: A *filename must be provided'
				fatal_errors += 1
			elif os.path.exists(value):
				fastq_file = value
			elif not os.path.exists(value):
				print 'Fatal Warning: %s does not exist' % value
				fatal_errors += 1
			else:
				pass
				
		elif variable == '*script_bin':
			if value == '':
				print 'Fatal Warning: The location of *script_bin must be provided'
				fatal_errors += 1
			elif os.path.exists(value):
				script_bin = value
			else:
				print 'Fatal Warning: The script_bin directory %s does not exist' % value
				fatal_errors += 1
				
		elif variable == '*whole_genome_bt2':
			if value == '':
				print 'Fatal Warning: The location of *whole_genome_bt2 must be provided'
				fatal_errors += 1
			elif os.path.exists('%s.1.bt2' % value):
 				whole_genome_bt2 = value
 			else:
 				print 'Fatal Warning: %s does not exist as a bowtie2 index' % value
 				fatal_errors += 1
				
		elif variable == '*whole_genome_lengths':
			if value == '':
				print 'Fatal Warning: The location of *whole_genome_lengths must be provided'
				fatal_errors += 1
			elif os.path.exists(value):
				whole_genome_len = value
			else:
				print 'Fatal Warning: %s does not exist' % value
				fatal_errors += 1
		
		elif variable == '*number_genomes':
			if value == '':
				print 'Fatal Warning: The *number_genomes must be provided'
				print 'Please list the number of genomes in your genome file'
				fatal_errors += 1
			else:
				try:
					value = int(value)
					if value <= 0:
						print "Fatal Warning: '*number_genomes' must be greater than zero, not %s" % value 
						fatal_errors += 1
					else:
						bowtie2_k = value
				except:
					print "Fatal Warning: '*number_genomes' must equal an integer, not %s" % value 
					fatal_errors += 1

		elif variable == '*noncore_bt2':
			if value == '':
				print 'Fatal Warning: The location of *noncore_bt2 must be provided'
				fatal_errors += 1
			elif os.path.exists(value):
				noncore_bt2 = value
			else:
				print 'Fatal Warning: %s does not exist' % value
				fatal_errors += 1
				
		elif variable == '*noncore_lengths_taxids':
			if value == '':
				print 'Fatal Warning: The location of *noncore_lengths_taxids must be provided'
				fatal_errors += 1
			elif os.path.exists(value):
				noncore_len = value
			else:
				print 'Fatal Warning: %s does not exist' % value
				fatal_errors += 1
		
		elif variable == '*SNPs':
			if value == '':
				print 'Fatal Warning: The location of *SNPs must be provided'
				fatal_errors += 1
			elif os.path.exists(value):
				SNP_file = value
			else:
				print 'Fatal Warning: %s does not exist' % value
				fatal_errors += 1
		
		elif variable == '*SNP_counts':
			if value == '':
				print 'Fatal Warning: The location of *SNPs_counts must be provided'
				fatal_errors += 1
			elif os.path.exists(value):
				SNP_count_file = value
			else:
				print 'Fatal Warning: %s does not exist' % value
				fatal_errors += 1

		elif variable == 'num_threads':
			if value == '':
				print "No value for 'num_threads' was provided"
				print "\tBy default the pipeline will be run with % s" % num_threads
			else:
				try:
					value = int(value)
					if value <= 0:
						print "Warning: 'num_threads' must be greater than zero, not %s" % value 
						print "\tBy default the pipeline will be run with num_threads = % s" % num_threads
					else:
						num_threads = value
				except:
					print "Warning: 'num_threads' must equal an integer, not %s" % value 
					print "\tBy default the pipeline will be run with num_threads = % s" % num_threads
						
		elif variable == 'output_directory':
			if value == '':
				cwd = os.getcwd()
				print "Warning: No valid output directory given"
				print "\tBy default results will be written to %s" % cwd
				output_dir = cwd
			elif os.path.exists(value):
				output_dir = value 
			else:
				cwd = os.getcwd()
				print 'Warning: The output bin directory %s does not exist' % value
				print "\tBy default results will be written to %s" % cwd
				output_dir = cwd

FILE_IN.close()

if fatal_errors > 0:
	print "Too many fatal errors: Exitting"
	sys.exit(1)
else:
	sample_name = fastq_file.split('/')[-1]
	sample_name = sample_name.replace('.fq','')
	sample_name = sample_name.replace('.fastq','')
	sample_name = sample_name.replace('.fa','')
	sample_name = sample_name.replace('.fasta','')

#writing the shell script 
print "\nWriting the shell file as %s/run_strain_tracking_%s.sh\n" % (output_dir,sample_name)

FILE_OUT = open('%s/run_strain_tracking_%s.sh' % (output_dir,sample_name), 'w')

FILE_OUT.write("echo 'Running strain tracking on sample %s'\n\n" % sample_name)

FILE_OUT.write("echo 'Running bowtie with complete strain sequences'\n")
sam = '%s/%s.sam' % (output_dir, sample_name)
FILE_OUT.write("bowtie2 -x %s -U %s -S %s -p %s %s -k %s\n" % (whole_genome_bt2, fastq_file, sam, num_threads, parameter, bowtie2_k) )

FILE_OUT.write("echo 'Counting SNPs'\n")
FILE_OUT.write("python %s/parse_SNPs_from_sam.py %s %s/%s_SNPs.txt %s/%s_SNPs_list.txt %s %s\n" % (script_bin, sam, output_dir, sample_name, output_dir, sample_name, SNP_file, SNP_count_file ) )

FILE_OUT.write("echo 'Running pathoscope'\n")
FILE_OUT.write("python %s/pathoscope/pathoscope2.py -t sam -e %s -f %s -outdir %s -noUpdatedAlignFile -thetaPrior 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000\n" % (script_bin, sample_name, sam, output_dir) )
FILE_OUT.write("rm %s\n" % sam)

FILE_OUT.write("echo 'Counting number of hits based of Pathoscope.tsv file'\n")
FILE_OUT.write("python %s/inferNumberReads_fromPatho_tsv.py %s/%s-sam-report.tsv %s/updated_%s.txt 1\n" % (script_bin, output_dir, sample_name, output_dir, sample_name ) )

FILE_OUT.write("echo 'Normalizing the read counts'\n")
FILE_OUT.write("python %s/readCountNormalizer.py %s/updated_%s.txt %s/updated_%s_norm.txt %s\n" % (script_bin, output_dir, sample_name, output_dir, sample_name, whole_genome_len) )
FILE_OUT.write("rm %s/updated_%s.txt\n\n" % (output_dir, sample_name) )

FILE_OUT.write("echo 'Running bowtie with noncore regions'\n")
sam = '%s/%s_pieces.sam' % (output_dir, sample_name)
FILE_OUT.write("bowtie2 -x %s -U %s -S %s -p %s %s -k %s\n" % (noncore_bt2, fastq_file, sam, num_threads, parameter, bowtie2_k) )

FILE_OUT.write("echo 'Running pathoscope on the pieces'\n")
FILE_OUT.write("python %s/pathoscope/pathoscope2.py -t sam -e %s_pieces -f %s -outdir %s -noUpdatedAlignFile\n" % (script_bin, sample_name, sam, output_dir) )
FILE_OUT.write("rm %s\n" % sam)

FILE_OUT.write("echo 'Counting number of hits based of Pathoscope.tsv file'\n")
FILE_OUT.write("python %s/inferNumberReads_fromPatho_tsv.py %s/%s_pieces-sam-report.tsv %s/updated_%s_pieces.txt 0\n" % (script_bin,output_dir, sample_name, output_dir, sample_name) )

FILE_OUT.write("echo 'Normalizing the counts based on the length of the segment'\n")
FILE_OUT.write("python %s/readCountNormalizer.py %s/updated_%s_pieces.txt  %s/updated_%s_pieces_norm.txt %s\n" % (script_bin, output_dir, sample_name, output_dir, sample_name, noncore_len) )

FILE_OUT.write("echo 'Parsing the data to have information for each strain'\n")
FILE_OUT.write("python %s/parse_pieces_out.py %s/updated_%s_pieces_norm.txt %s/updated_%s_pieces_norm_parsed.txt %s\n" % (script_bin,output_dir, sample_name, output_dir, sample_name, noncore_len) )
FILE_OUT.write("rm %s/updated_%s_pieces.txt %s/updated_%s_pieces_norm.txt\n\n" % (output_dir, sample_name, output_dir, sample_name) )

FILE_OUT.write("echo 'Finished!'\n")

FILE_OUT.close()

















