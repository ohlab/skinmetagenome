#Requires Bowtie2, python2.7 or greater, Samtools, Bedtools, and Pathoscope v1. 
#Designed for a GNU/Linux parallel processing system, e.g., here, 16 processors.
#classify.sh and classify_singlet.sh are the same except classify.sh is designed for submission of a large number of files for a PBS cluster system and classify_singlet.sh can be run on a desktop machine with adequate hard drive space for large SAM files and bowtie indices. 

#fastq files are contained in ./fastq
#databases are in ./dbs
#binaries are in your own defined directory


#make an output directory
mkdir output
echo "Removing human reads"
./bowtie2-2.2.3/bowtie2 -x ./dbs/hg19_rRNA -U  ./fastq/FILENAME --very-sensitive -S ./output/FILENAME\_human\.sam -p 16 

python2.7 ./bin/updatedfastq.py ./output/FILENAME\_human\.sam ./fastq/FILENAME ./output/FILENAME\_NOThuman\.fastq ./output/FILENAME\_NOThuman\.txt 2

rm ./output/FILENAME\_NOThuman\.txt
rm ./output/FILENAME\_human\.sam

echo "Running bowtie on Fungal databases"
./bowtie2-2.2.3/bowtie2 -x ./dbs/all_fungiA -U ./output/FILENAME\_NOThuman\.fastq --very-sensitive -k 10 -S ./output/FILENAME\_fungiA.sam -p 16 
./bowtie2-2.2.3/bowtie2 -x ./dbs/all_fungiB -U ./output/FILENAME\_NOThuman\.fastq --very-sensitive -k 10 -S ./output/FILENAME\_fungiB.sam -p 16 

echo "Running bowtie on Viral databases"
./bowtie2-2.2.3/bowtie2 -x ./dbs/virus.fa -U ./output/FILENAME\_NOThuman\.fastq --very-sensitive -k 10 -S ./output/FILENAME\_virus.sam -p 16 

echo "Running bowtie on Archaeal databases "
./bowtie2-2.2.3/bowtie2 -x ./dbs/hmp_arch_merged.fa -U ./output/FILENAME\_NOThuman\.fastq --very-sensitive -k 10 -S ./output/FILENAME\_arch.sam -p 16 

echo "Running bowtie on Bacterial databases"
./bowtie2-2.2.3/bowtie2 -x ./dbs/all_bacB_1B -U ./output/FILENAME\_NOThuman\.fastq --very-sensitive -k 10 -S ./output/FILENAME\_bac1.sam -p 16 
./bowtie2-2.2.3/bowtie2 -x ./dbs/all_bacB_2 -U ./output/FILENAME\_NOThuman\.fastq --very-sensitive -k 10 -S ./output/FILENAME\_bac2.sam -p 16 
./bowtie2-2.2.3/bowtie2 -x ./dbs/all_bacB_3B -U ./output/FILENAME\_NOThuman\.fastq --very-sensitive -k 10 -S ./output/FILENAME\_bac3.sam -p 16 

rm ./output/FILENAME\_NOThuman\.fastq

echo "Combining sam with multiple subtypes"
echo "Looking for overlap between reads in updated sam file"
python2.7 ./bin/combine_SamsAll.py ./output/FILENAME\_multipleSubtypes_noHeaders.sam  ./output/FILENAME\_headers.txt ./output/FILENAME\_arch.sam ./output/FILENAME\_fungiA.sam ./output/FILENAME\_virus.sam ./output/FILENAME\_fungiB.sam  ./output/FILENAME\_bac1.sam ./output/FILENAME\_bac2.sam ./output/FILENAME\_bac3.sam

echo "Counting number of hits for each subject in the combined sam file"
python2.7 ./bin/countReads_sorted_pathoscope_updated.py ./output/FILENAME\_multipleSubtypes_noHeaders.sam ./output/FILENAME\_multipleSubtypes.txt

cat ./output/FILENAME\_headers.txt ./output/FILENAME\_multipleSubtypes_noHeaders.sam > ./output/FILENAME\_multipleSubtypes.sam 

echo "deleting all the individual sam files"
rm ./output/FILENAME\_headers.txt
rm ./output/FILENAME\_multipleSubtypes_noHeaders.sam
rm ./output/FILENAME\_fungiA.sam 
rm ./output/FILENAME\_virus.sam
rm ./output/FILENAME\_fungiB.sam
rm ./output/FILENAME\_arch.sam
rm ./output/FILENAME\_bac1.sam
rm ./output/FILENAME\_bac2.sam
rm ./output/FILENAME\_bac3.sam

echo "Running pathoscope"
python2.7 ./bin/pathoscope/pathoscope.py -t sam -e FILENAME\_multipleSubtypes.sam -f ./output/FILENAME\_multipleSubtypes.sam -outdir ./output

rm ./output/FILENAME\_multipleSubtypes.sam

echo "Counting number of hits for each subject in the updated sam file"
python2.7 ./bin/countReads_sorted_pathoscope_updated.py ./output/updated_FILENAME\_multipleSubtypes.sam ./output/updated_FILENAME\_multipleSubtypes.txt

echo "Looking for overlap between reads in updated sam file"
python2.7 ./bin/overlap_betweenKingdoms.py ./output/updated_FILENAME\_multipleSubtypes.sam ./bin/all_taxids_lineages.txt

echo "Running coverage calculation on updated sam file"

echo "Converting sam to bam"
samtools view -bS ./output/updated_FILENAME\_multipleSubtypes.sam > ./output/updated_FILENAME\_multipleSubtypes.bam  

echo "Sorting bam"
./samtools/bin/samtools sort ./output/updated_FILENAME\_multipleSubtypes.bam ./output/updated_FILENAME\_multipleSubtypes.sorted   

echo "Running genomeCoverageBed"
./bedtools-2.17.0/bin/genomeCoverageBed -ibam ./output/updated_FILENAME\_multipleSubtypes.sorted.bam > ./output/updated_FILENAME\_multipleSubtypes.sorted_coverages.txt

rm ./output/updated_FILENAME\_multipleSubtypes.sam
rm ./output/updated_FILENAME\_multipleSubtypes.bam
rm ./output/updated_FILENAME\_multipleSubtypes.sorted.bam

echo "Calculating the coverage for each of the subjects"
python2.7 ./bin/assignPercentCoverage.py ./output/updated_FILENAME\_multipleSubtypes.txt ./output/updated_FILENAME\_multipleSubtypes.sorted_coverages.txt ./output/updated_FILENAME\_multipleSubtypes_coverages.txt 5  

echo "Normalizing the read counts"
python2.7 ./bin/readCountNormalizer.py ./output/updated_FILENAME\_multipleSubtypes_coverages.txt ./output/updated_FILENAME\_multipleSubtypes_coverages_normalized.txt ./bin/all_lengths.txt 

sort -k 4 -g -r ./output/updated_FILENAME\_multipleSubtypes_coverages_normalized.txt > ./output/updated_FILENAME\_multipleSubtypes_coverages_normalized_sorted.txt

rm ./output/updated_FILENAME\_multipleSubtypes_coverages_normalized.txt
rm ./output/updated_FILENAME\_multipleSubtypes.txt 
rm ./output/updated_FILENAME\_multipleSubtypes.sorted_coverages.txt
rm ./output/updated_FILENAME\_multipleSubtypes_coverages.txt

echo "--complete"
