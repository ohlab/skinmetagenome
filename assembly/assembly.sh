#Requires Velvet, Bowtie2, Python2.7, khmer, and seqtk
#Designed for a GNU/Linux parallel processing system, e.g., here, 16 processors and at least 24 GB of RAM. 
#binaries are in your own defined directory
#fastq files are in ./fastq

#this shell script expects a list of processed fastq files as input. 
ls ./fastq/*.fastq > filelist.txt

#uninterleave your files, and make an interleaved paired end file. 
for f in `cat filelist.txt`
do
java -jar ./bin/uninterleave.jar -fastq ./fastq/$f;
./bin/shuffleSequences_fastq.pl ./fastq/$f.mate1 ./fastq/$f.mate2 ./fastq/$f.interleaved.pe;
done

#make a velvet temp script. 33-69 is the kmer range interrogated. Assembly is followed by mapping reads back to contigs with Bowtie2. 
nano velvet.sh
for i in `seq 33 4 69`; 
do
TMPSCRIPT="FOO.velvet.$i.qsub"
echo "#\$ -S ./bin/bash" > $TMPSCRIPT
echo "#\$" >> $TMPSCRIPT
echo "cd FOLDER;" >> $TMPSCRIPT
echo "./velvet_1.2.10/velveth FOO_velvet.assembly.$i $i -fastq -short ./fastq/FOO.singletons -shortPaired ./fastq/FOO.interleaved.pe;" >> $TMPSCRIPT
echo "./velvet_1.2.10/velvetg FOO_velvet.assembly.$i -exp_cov auto -cov_cutoff auto -scaffolding no -read_trkg yes;" >> $TMPSCRIPT
echo "python2.7 ./khmer/sandbox/extract-long-sequences.py 300 FOO_velvet.assembly.$i/contigs.fa > FOO_velvet_assembly.$i.fa;" >> $TMPSCRIPT
echo "rm -r FOO_velvet.assembly.$i;" >> $TMPSCRIPT
echo "./bowtie2-2.2.3/bowtie2-build FOO_velvet_assembly.$i.fa FOO.$i.velvet_assembly;" >> $TMPSCRIPT
echo "./bowtie2-2.2.3/bowtie2 -p 16 -x FOO.$i.velvet_assembly -q -1 ./fastq/FOO.mate1 -2 ./fastq/FOO.mate2 -U ./fastq/FOO.singletons --sensitive --no-hd --no-unal -S FOO.$i.velvet.map > FOO.$i.velvet.stats.txt 2>&1;" >> $TMPSCRIPT
echo "rm FOO.$i.velvet.map;" >> $TMPSCRIPT
echo "rm FOO.$i.velvet_assembly*bt2;" >> $TMPSCRIPT
qsub -l nodes=1:c16 $TMPSCRIPT
done

#make a list of your interleaved files to process
ls ./fastq/*.interleaved.pe | sed 's/.interleaved.pe//g' > interleaved.txt

#run all the kmer ranges for each of your fastq files. 
for f in `cat interleaved.txt`
do
FOLDER="your folder here"
sed s/FOO/$f/g velvet.sh | sed s/$g > $f.velvet_running.sh
sh $f.velvet_running.sh
done


#once all the runs are complete, get some assembly statistics: 
for f in *.fa
do
python2.7 ./khmer/sandbox/assemstats3.py 300 $f >> all_assemstats.txt
done
mv all_assemstats.txt temp_all_assemstats.txt
grep .fa temp_all_assemstats.txt >  all_assemstats.txt

#now, get bowtie stats. 
grep -r "reads" *.stats.txt > totalreads_stats.txt
grep -r "were paired" *.stats.txt > bowtie.numpairedreads_stats.txt
grep -r "overall" *.stats.txt > overall.overallalign_stats.txt
grep -r "aligned concordantly exactly 1 time" *.stats.txt > concordantapair_align_stats.txt
grep -r "aligned discordantly" *.stats.txt > discordantpair_align_stats.txt

paste concordantapair_align_stats.txt discordantpair_align_stats.txt overall.overallalign_stats.txt bowtie.numpairedreads_stats.txt totalreads_stats.txt> all.txt
sed 's/ \{2,\}//g' all.txt > bowtiestats.txt

#Take files all_assemstats.txt and bowtiestats.txt. This will output a file "velvet_filestoget.txt" which is the "best" assembly based on your desired parameters for assembly optimization. Here, % alignment of reads back to the assembly is the most heavily weighted. 
R CMD BATCH ./bin/velvet_stats.R

#now, these are the files to get; make a shell file to move if you want:
sed 's/"//g' velvet_filestoget.txt | sed 's/^[A-Z]/mv &/g' | sed 's/.fa$/.fa selected/g' | sed 's/_assembly/*_assembly/'> velvet_filestoget.sh
mkdir selected
sh velvet_filestoget.sh

#now, let's get all aligned reads.
#get a list of your reads to use.
ls *.fa | sed 's/\.20_velvet_assembly//' | sed 's/_velvet_assembly//' | cut -f 1 -d '.' > ./selected/best.txt

for f in `cat best.txt`
do
TMPSCRIPT="$f.getunaligned.qsub"
echo "#\$ -S /bin/bash" > $TMPSCRIPT
echo "#\$" >> $TMPSCRIPT
echo "cd FOLDER;" >> $TMPSCRIPT
echo "./bowtie2-2.2.3/bowtie2-build $f*_velvet_assembly.*.fa $f.ref;" >> $TMPSCRIPT
echo "./bowtie2-2.2.3/bowtie2 -p 16 -x $f.ref -q -U ./fastq/$f.fastq --sensitive --no-hd --un $f.unaligned -S $f.velvet.map > $f.velvet.stats2.txt 2>&1;" >> $TMPSCRIPT
qsub -l nodes=1:c16 $TMPSCRIPT
done

##############################################################
#now, we combine unaligned reads. Depending on the size of your dataset, this can become very enormous. In our paper, we ran 2 iterative rounds of pooling unaligned reads, first combining unaligned reads that belong to an individual, then pooling all aligned reads. This is an example of how we run digital normalization with khmer to reduce the sample complexity prior to adaptive assembly. 

cat *.unaligned > all.unaligned

#if you wish to subsample your data prior to digital normalization: 
./seqtk-master/seqtk sample -s100 all.unaligned NUMBEROFREADS > all.subsampled.unaligned


nano velvet_khmer.sh
python2.7 ./khmer/scripts/normalize-by-median.py -C 20 -k 20 -N 4 -x 16e9 --savehash all.unaligned_ref.kh -R all.unaligned_1.report all.unaligned; 
python2.7 ./khmer/scripts/filter-abund.py -V all.unaligned_ref.kh all.unaligned.keep; 
python2.7 ./khmer/scripts/normalize-by-median.py -C 5 -k 20 -N 4 -x 16e9 all.unaligned.keep.abundfilt; 
python2.7 ./khmer/sandbox/strip-and-split-for-assembly.py all.unaligned.keep.abundfilt.keep;
python2.7 ./khmer/sandbox/readstats.py all.unaligned*.[e,p]e;
java -jar ./bin/uninterleave.jar -fastq all.unaligned

TMPSCRIPT="FOO.postkhmer_velvet.$i.qsub"
echo "#\$ -S /bin/bash" > $TMPSCRIPT
echo "#\$" >> $TMPSCRIPT
echo "cd FOLDER;" >> $TMPSCRIPT
for i in `seq 25 4 61`; 
do
echo "./velvet_1.2.10/velveth all.unaligned_postkhmer_velvet.assembly.$i $i -fasta -short all.unaligned.keep.abundfilt.keep.se -shortPaired all.unaligned.keep.abundfilt.keep.pe;" >> $TMPSCRIPT
echo "./velvet_1.2.10/velvetg all.unaligned_postkhmer_velvet.assembly.$i -exp_cov auto -cov_cutoff auto -scaffolding no -read_trkg yes;" >> $TMPSCRIPT
echo "python2.7 ./khmer/sandbox/extract-long-sequences.py 300 all.unaligned_postkhmer_velvet.assembly.$i/contigs.fa > all.unaligned_postkhmer_velvet_assembly.$i.fa;" >> $TMPSCRIPT
echo "rm -r all.unaligned_postkhmer_velvet.assembly.$i;" >> $TMPSCRIPT

echo "./bowtie2-2.2.3/bowtie2-build all.unaligned_postkhmer_velvet_assembly.$i.fa all.unaligned.$i.postkhmer_velvet_assembly;" >> $TMPSCRIPT
echo "./bowtie2-2.2.3/bowtie2 -p 16 -x all.unaligned.$i.postkhmer_velvet_assembly -q -1 all.unaligned.mate1 -2 all.unaligned.mate2 -U all.unaligned.singletons --sensitive --no-hd --no-unal -S all.unaligned.$i.postkhmer_velvet.map > all.unaligned.$i.postkhmer_velvet.stats.txt 2>&1;" >> $TMPSCRIPT
echo "rm FOO.$i.postkhmer_velvet.map;" >> $TMPSCRIPT
echo "rm FOO.$i.postkhmer_velvet_assembly*bt2;" >> $TMPSCRIPT
done
qsub -l nodes=1:g72 $TMPSCRIPT
#end

sh velvet_khmer.sh

#and repeat steps to get assembly statistics from above. 

