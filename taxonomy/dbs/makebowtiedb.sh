#make bowtie indices for taxonomic mapping  
tar tvf dbs_fa.tar.gz;	#these are the original fasta files. 

cd dbs
for f in *.fa
do
base=`basename $f .fa`
./bin/bowtie2-2.2.3/bowtie2-build $f $base;
done
