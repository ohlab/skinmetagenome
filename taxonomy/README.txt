This is a version of the taxonomic mapping that runs on a GNU/Linux parallel processing system. Bacterial, fungal, viral, and archaeal genomes were downloaded from NCBI, the Human Microbiome Project, Saccharomyces Genome Database, Fungal Genome Initiative, FungiDB. 


DISTRIBUTION:

This file				README.txt
Main shell script		classify.sh
						classify_singlet.sh
Helper scripts			./bin/pathoscope/	#version 1.0
						./bin/updatedfastq.py
						./bin/combine_SamsAll.py
						./bin/countReads_sorted_pathoscope_updated.py
						./bin/overlap_betweenKingdoms.py
						./bin/assignPercentCoverage.py
						./bin/readCountNormalizer.py 
Helper reference files	./bin/all_lengths.txt
						./bin/all_taxids_lineages.txt
Reference packages		./dbs <--https://drive.google.com/file/d/0Bx2MDtKknM1eUjJETkJZRFYzaEE/edit?usp=sharing
#please note, we are not updating this database package and so we recommend an expiration date of 1 year post publication. 


USAGE:

1. Download and install Pathoscope 1.0, Bowtie2, Python2.7, Samtools, Bedtools, Perl
2. Download and unpack the database reference files (current last update 3/24/2014), and make bowtie indices. This will take a long time upfront but only needs to be done once. Information on how to add genomes to the database is provided in ./dbs/README.txt. We encourage you to share your new databases :) 
 
https://drive.google.com/file/d/0Bx2MDtKknM1eUjJETkJZRFYzaEE/edit?usp=sharing
sh makebowtiedb.sh

3. If you are running a small number of files locally, use classify_singlets.sh
for f in *.fastq
do
sed 's/FILENAME/$f/g' classify_singlet.sh > $f.classify.sh;
nohup time sh $f.classify.sh > nohup.$f.log 2>&1&
done 

4. If you are using e.g., a GNU/Linux parallel processing system, use classify.sh. Make a file containing your fastq files called samples.txt, then run pipeline using
perl classify.sh samples.txt

5. Final files (updated_FILE_multipleSubtypes_coverages_normalized_sorted.txt) are a list of taxonomic classifications with relative abundances normalized by genome length, and an estimated genome coverage. 

6. Note, all reads are treated as singletons so you don't have to account for paired-end reads. 

 



#create a list of your fastq files. Here, all reads are treated as singletons so you do not have to account for pairs. 

ls ./fastq/*.fastq > filelist.txt

for f in `cat filelist.txt`
do 
sed 

