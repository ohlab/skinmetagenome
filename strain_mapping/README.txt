Everything necessary for comparative genomic analysis and strain tracking of metagenomic samples 

CONTENTS:

	A) Preparing genomes:
		•giving_genomes_useable_names.py:  
			Take genome files downloaded from NCBI and rename them to have useable names 

	B) Comparative genomics:
		Main script:
			•comparative_genomics.py:
				Runs necessary scripts for comparative genomic analysis of multiple genomes
		Accessory scripts called:	
			•associate_pieces_with_taxid.py: 
				Associate pieces with a taxid and name, so pieces from the same organism can be easily associated with one another
			•count_number_SNPs_per_strain.py: 
				Counts number of unique SNPs per strain
			•extract_regions.py: 	
				Writes out genomes where all coordinates in the input file have been extracted in one file and reverse coordinates are masked in another
			•find_lengths.py:
				Calculates the length of each sequence in a fasta file
			•make_core_SNP_string.py: 
				To create a multifasta file from the SNPs identified in the core 
			•parse_aligned_file.py: 
				Simplifies file of redundant alignment coordinates to non-repetitive regions
			•parse_SNPs.py: 
				Parses SNP output from SNV_table_2.pl and identify SNPs in the core unique to a particular strain 
			•reverse_alignments.py: 
				Reverses list of alignment coordinates using the length of the segment
			•SNV_table_2.pl: 
				Runs nucmer comparisons 
	
	C) Strain tracking pipeline:
		Main script:
			•run_strain_tacking.sh:
				Runs pipeline for strain tracking 
		Accessory scripts called:
			•inferNumberReads_fromPatho_tsv.py
				Infer number of reads from the patho .tsv file, genomes which fall below a certain percentage are grouped into "Other"
			•parse_SNPs_from_sam.py
				Parse a bowtie2 samfile and identify reads mapping to SNPs informative to a strain 
			•parse_pieces_out.py
				Parse the counts output from the noncore pieces, so pieces from the same organism are summed together
			•readCountNormalizer.py
				 Normalize percentages by the total number of reads and the length of the sequence in the database

INSTRUCTIONS:

	A) Preparing genomes for comparative genomic analysis:
		Requirements: python
		1. Download genomes you wish to compare from NCBI
			Notes: 
				-the fasta file headers need to be in the format: >gi|50841496|ref|NC_006085.1| Propionibacterium acnes KPA171202 chromosome, complete genome
				-all contigs from the same genome need to be in the same fasta file (the ncbi ftp has them split into separate files)
				
		2. Download 2 NCBI taxonomy files to the bin/NCBI_taxonomy folder
			wget ftp://anonymous@ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
			wget ftp://anonymous@ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
			
		3. Rename genomes downloaded from NCBI to something more useable 
			$python bin/giving_genomes_useable_names.py genome_org/ genome/ bin/NCBI_taxonomy/gi_taxid_nucl.dmp bin/NCBI_taxonomy/names.dmp T
			Notes: 
				-For example: genome_org/NC_006085.fna is renamed to genome/Propionibacterium_acnes_KPA171202.fa
				-If the final parameter equals T, merged versions of the genomes are generated in genome/merge. 
					Merged indicates that a genome with multiple fasta sequences has been artificially assembled into a single contiguous sequence with strings of Ns.
					Merged genomes are necessary for the Pathoscope portion of the strain tracking pipeline. 

	B) Comparative genomics analysis (for additional information please read notes in comparative_genomics.py):
		**All of these commands need to be executed in the desired output folder
		Requirements: mummer, python
			-both need to be in your $PATH
		1. Identify core region in reference genome
			$python ../bin/comparative_genomics.py reference.fa CORE script_dir genome_dir &> log.txt
				Notes:
					-reference.fa should be a genome that's fully assembled 
					-log.txt can be processed by bin/parse_log.py for easier viewing 
				 	$python ../bin/parse_log.py log.txt log_parsed.txt
			
		1B. Identify core regions when leaving genomes out of the analysis (necessary if you want to exclude a genome in genome_dir from the analysis) 
			$python comparative_genomics.py reference.fa CORE script_dir genome_dir leave_out_1.fa leave_out_2.fa &> log.txt
		
		2. Identify SNPs in core regions
			$python ../bin/comparative_genomics.py reference.fa SNP script_dir genome_dir
				Notes:
					-Uses version of the reference genome where noncore regions have been masked (genome_dir/core/reference.fa_masked)
					-Uses the merged genomes because these are used in the strain tracking pipeline and we want the coordinates to match 
		3. Identify noncore regions in the genomes 
			$python ../bin/comparative_genomics.py reference.fa SNP script_dir genome_dir
				Notes:
					-Uses version of the reference genome where noncore regions have been masked (genome_dir/core/reference.fa_masked)

	C)Strain tracking pipeline
		Requirements: python 2.7 or higher and bowtie2 
			-both need to be in your $PATH
		1. Process merged genomes for strain tracking
			-concatenate all merged genomes into one file
			$cat genome_dir/merged/* > genome_dir/merged/all_genomes.fa
			-find the length of each genome in the file, used for normalization
			$python bin/find_lengths.py genome_dir/merged/all_genomes.fa genome_dir/merged/all_genomes.fa_lengths
		
		2. Create bowtie2 indices of the merged genomes and noncore pieces
			$bowtie2-build genome_dir/all_genomes.fa genome_dir/all_genomes.fa
			$bowtie2-build genome_dir/noncore/all_noncore_pieces.fa genome_dir/noncore/all_noncore_pieces.fa
		
		3. Fill in config.txt 
		
		4. Generate the shell script 
			$python run_strain_tracking.py config.txt
		
		5. Run the generated shell script
			$sh run_strain_tracking_sample.sh
					 
			
module load mummer 

time python bin/giving_genomes_useable_names.py genome_org/ genome/ bin/NCBI_taxonomy/gi_taxid_nucl.dmp bin/NCBI_taxonomy/names.dmp T

#the following commands need to be executed in the output folder 
module load mummer
time python ../bin/comparative_genomics.py Propionibacterium_acnes_KPA171202.fa CORE ../bin/ ../genome Propionibacterium_acnes_J139.fa &> log.txt
time python ../bin/parse_log.py log.txt log_parsed.txt
time python ../bin/comparative_genomics.py Propionibacterium_acnes_KPA171202.fa SNP ../bin/ ../genome/ Propionibacterium_acnes_J139.fa
time python ../bin/comparative_genomics.py Propionibacterium_acnes_KPA171202.fa NONCORE ../bin/ ../genome/ Propionibacterium_acnes_J139.fa

 


