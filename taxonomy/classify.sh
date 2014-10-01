#Requires Bowtie2, Python2.7, Bowtie2, Samtools, Bedtools, and Pathoscope v1. 
#Designed for a GNU/Linux parallel processing system, e.g., here, 16 processors.

#classify.sh and classify_singlet.sh are the same except classify.sh is designed for submission of a large number of files for a PBS cluster system and classify_singlet.sh can be run on a desktop machine with adequate hard drive space for large SAM files and bowtie indices. 

#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Basename;

# This Perl script creates a .qsub file (a bash script) and submits it to the queue.
# This is a template based on the original files by Andi Broka.

# Get the file name of the run file from the command line
my ($list_of_sequencing_files_name) = @ARGV;

open(RUNS, $list_of_sequencing_files_name);

my @array_of_sequencing_files = <RUNS>;
close(RUNS);

# Run QA_stats for each fastq file (with directory) given in the run text file

foreach my $cur_file_name (@array_of_sequencing_files) {
    chomp($cur_file_name);
    my $cur_file_whole_name = $cur_file_name;
    my @name = split(/[.]/,$cur_file_name);
	my $cur_file = $name[0];
	
# Assign names to job and to the directory that will hold the .qsub and .log files for the job, and use this to name the .qsub file itself
my $source_dir = "./fastq";		#assign locations of files
my $db_dir = "./dbs";			#assign locations of bowtie indices
my $scrips_dir = "./bin";		#assign script locations
my $out_dir = "./output";		#designate an output directory
my $num_threads = "16"; 		#designate number of processors to use
my $k_param = "10"; 			#designate number of hits retrieved
my $qsub_dir = "./qsub";		#designate location to save scripts
my $log_dir = "./logs";			#designate location to save logs
my $jname = "$cur_file"; 
my $qsubf = "$cur_file\_metagenomics.sh";

# Open the .qsub file for writing and then start dumping in commands until 'EOF'
open (QSUBFILE, ">$qsub_dir/$qsubf");
print (QSUBFILE <<EOF);


#!/bin/bash
#
# this file is classify.sh
#
#PBS -N $cur_file\_meta
#  PBS -m be
#PBS -j eo
#PBS -e $log_dir/logfile_$cur_file

echo "=========================================================="
echo "Starting on       : \$(date)"
echo "Running on node   : \$PBS_O_HOST"
echo "Current job ID    : \$PBS_JOBID"
echo "Current job name  : \$PBS_JOBNAME"
echo "=========================================================="
echo " "

module load bedtools
module load samtools 
module load bowtie 
module load python 

echo "Removing human reads"
time bowtie2 -x $db_dir/hg19_rRNA -U  $source_dir/$cur_file_whole_name --very-sensitive -S $out_dir/$cur_file\_human\.sam -p $num_threads 

time python $scrips_dir/updatedfastq.py $out_dir/$cur_file\_human\.sam $source_dir/$cur_file_whole_name $out_dir/$cur_file\_NOThuman\.fastq $out_dir/$cur_file\_NOThuman\.txt 2

rm $out_dir/$cur_file\_NOThuman\.txt
rm $out_dir/$cur_file\_human\.sam

echo "Running bowtie on Fungi A "
time bowtie2 -x $db_dir/all_fungiA -U $out_dir/$cur_file\_NOThuman\.fastq --very-sensitive -k $k_param -S $out_dir/$cur_file\_fungiA.sam -p $num_threads 

echo "Running bowtie on Fungi B"
time bowtie2 -x $db_dir/all_fungiB -U $out_dir/$cur_file\_NOThuman\.fastq --very-sensitive -k $k_param -S $out_dir/$cur_file\_fungiB.sam -p $num_threads 

echo "Running bowtie on NCBI virus"
time bowtie2 -x $db_dir/virus.fa -U $out_dir/$cur_file\_NOThuman\.fastq --very-sensitive -k $k_param -S $out_dir/$cur_file\_virus.sam -p $num_threads 

echo "Running bowtie on arch"
time bowtie2 -x $db_dir/hmp_arch_merged.fa -U $out_dir/$cur_file\_NOThuman\.fastq --very-sensitive -k $k_param -S $out_dir/$cur_file\_arch.sam -p $num_threads 

echo "Running bowtie on NCBI/HMP bac 1"
time bowtie2 -x $db_dir/all_bacB_1B -U $out_dir/$cur_file\_NOThuman\.fastq --very-sensitive -k $k_param -S $out_dir/$cur_file\_bac1.sam -p $num_threads 

echo "Running bowtie on NCBI/HMP bac 2"
time bowtie2 -x $db_dir/all_bacB_2 -U $out_dir/$cur_file\_NOThuman\.fastq --very-sensitive -k $k_param -S $out_dir/$cur_file\_bac2.sam -p $num_threads 

echo "Running bowtie on NCBI/HMP bac 3"
time bowtie2 -x $db_dir/all_bacB_3B -U $out_dir/$cur_file\_NOThuman\.fastq --very-sensitive -k $k_param -S $out_dir/$cur_file\_bac3.sam -p $num_threads 

rm $out_dir/$cur_file\_NOThuman\.fastq

echo "Combining sam with multiple subtypes"
echo "Looking for overlap between reads in updated sam file"
python $scrips_dir/combine_SamsAll.py $out_dir/$cur_file\_multipleSubtypes_noHeaders.sam  $out_dir/$cur_file\_headers.txt $out_dir/$cur_file\_arch.sam $out_dir/$cur_file\_fungiA.sam $out_dir/$cur_file\_virus.sam $out_dir/$cur_file\_fungiB.sam  $out_dir/$cur_file\_bac1.sam $out_dir/$cur_file\_bac2.sam $out_dir/$cur_file\_bac3.sam

echo "Counting number of hits for each subject in the combined sam file"
time python $scrips_dir/countReads_sorted_pathoscope_updated.py $out_dir/$cur_file\_multipleSubtypes_noHeaders.sam $out_dir/$cur_file\_multipleSubtypes.txt

cat $out_dir/$cur_file\_headers.txt $out_dir/$cur_file\_multipleSubtypes_noHeaders.sam > $out_dir/$cur_file\_multipleSubtypes.sam 

echo "deleting all the individual sam files"
rm $out_dir/$cur_file\_headers.txt
rm $out_dir/$cur_file\_multipleSubtypes_noHeaders.sam
rm $out_dir/$cur_file\_fungiA.sam 
rm $out_dir/$cur_file\_virus.sam
rm $out_dir/$cur_file\_fungiB.sam
rm $out_dir/$cur_file\_arch.sam
rm $out_dir/$cur_file\_bac1.sam
rm $out_dir/$cur_file\_bac2.sam
rm $out_dir/$cur_file\_bac3.sam

echo "Running pathoscope"
time python ./bin/pathoscope/pathoscope.py -t sam -e $cur_file\_multipleSubtypes.sam -f $out_dir/$cur_file\_multipleSubtypes.sam -outdir $out_dir

rm $out_dir/$cur_file\_multipleSubtypes.sam

echo "Counting number of hits for each subject in the updated sam file"
time python $scrips_dir/countReads_sorted_pathoscope_updated.py $out_dir/updated_$cur_file\_multipleSubtypes.sam $out_dir/updated_$cur_file\_multipleSubtypes.txt

echo "Looking for overlap between reads in updated sam file"
time python $scrips_dir/overlap_betweenKingdoms.py $out_dir/updated_$cur_file\_multipleSubtypes.sam $scrips_dir/all_taxids_lineages.txt

echo "Running coverage calculation on updated sam file"

echo "Converting sam to bam"
time samtools view -bS $out_dir/updated_$cur_file\_multipleSubtypes.sam > $out_dir/updated_$cur_file\_multipleSubtypes.bam  

echo "Sorting bam"
time samtools sort $out_dir/updated_$cur_file\_multipleSubtypes.bam $out_dir/updated_$cur_file\_multipleSubtypes.sorted   

echo "Running genomeCoverageBed"
time genomeCoverageBed -ibam $out_dir/updated_$cur_file\_multipleSubtypes.sorted.bam > $out_dir/updated_$cur_file\_multipleSubtypes.sorted_coverages.txt

rm $out_dir/updated_$cur_file\_multipleSubtypes.sam
rm $out_dir/updated_$cur_file\_multipleSubtypes.bam
rm $out_dir/updated_$cur_file\_multipleSubtypes.sorted.bam

echo "Calculating the coverage for each of the subjects"

time python $scrips_dir/assignPercentCoverage.py $out_dir/updated_$cur_file\_multipleSubtypes.txt $out_dir/updated_$cur_file\_multipleSubtypes.sorted_coverages.txt $out_dir/updated_$cur_file\_multipleSubtypes_coverages.txt 5  

echo "Normalizing the read counts"
time python $scrips_dir/readCountNormalizer.py $out_dir/updated_$cur_file\_multipleSubtypes_coverages.txt $out_dir/updated_$cur_file\_multipleSubtypes_coverages_normalized.txt $scrips_dir/all_lengths.txt 

sort -k 4 -g -r $out_dir/updated_$cur_file\_multipleSubtypes_coverages_normalized.txt > $out_dir/updated_$cur_file\_multipleSubtypes_coverages_normalized_sorted.txt

rm $out_dir/updated_$cur_file\_multipleSubtypes_coverages_normalized.txt
rm $out_dir/updated_$cur_file\_multipleSubtypes.txt 
rm $out_dir/updated_$cur_file\_multipleSubtypes.sorted_coverages.txt
rm $out_dir/updated_$cur_file\_multipleSubtypes_coverages.txt

echo "--complete"

echo "=========================================================="
echo "Finished on       : \$(date)"
echo "=========================================================="

EOF

# Close the .qsub file and submit it to the queue
close (QSUBFILE);
system ("qsub -l nodes=1:c16 $qsub_dir/$qsubf");
}
