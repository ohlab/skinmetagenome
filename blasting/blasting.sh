#this generates the blastn and blastx searchables. fastacmd from legacy versions (pre Blast+) is used as is BLAST 2.2.26, and USEARCH7.0 or greater. 
#this script follows a workflow to blastn against reference genomes, then blastx sequences that don't match, and then gather the composite taxonomies. 
#binaries are in your own defined directory
#fastq files are in ./fastq

#download nr from the NCBI website and GI file
wget ftp.ncbi.nih.gov/blast/db/nr.00.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.01.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.02.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.03.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.04.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.05.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.06.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.07.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.08.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.09.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.10.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.11.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.12.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.13.tar.gz
wget ftp.ncbi.nih.gov/blast/db/nr.14.tar.gz
wget ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz

wget ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz

for f in *.tar.gz
do
tar xvf $f
done

#do some formatting
fastacmd -D 1 -d ./blastdb/nr.00 -o nr.00
fastacmd -D 1 -d ./blastdb/nr.01 -o nr.01
fastacmd -D 1 -d ./blastdb/nr.02 -o nr.02
fastacmd -D 1 -d ./blastdb/nr.03 -o nr.03
fastacmd -D 1 -d ./blastdb/nr.04 -o nr.04
fastacmd -D 1 -d ./blastdb/nr.05 -o nr.05
fastacmd -D 1 -d ./blastdb/nr.06 -o nr.06
fastacmd -D 1 -d ./blastdb/nr.07 -o nr.07
fastacmd -D 1 -d ./blastdb/nr.08 -o nr.08
fastacmd -D 1 -d ./blastdb/nr.09 -o nr.09
fastacmd -D 1 -d ./blastdb/nr.10 -o nr.10
fastacmd -D 1 -d ./blastdb/nr.11 -o nr.11
fastacmd -D 1 -d ./blastdb/nr.12 -o nr.12
fastacmd -D 1 -d ./blastdb/nr.13 -o nr.13
fastacmd -D 1 -d ./blastdb/nr.14 -o nr.14

ls nr.*[0-9]$ > nrlist.txt

#make a udb file
for f in `cat nrlist.txt`
do
TMPSCRIPT="usearch.$f.qsub"
echo "#\$ -S /bin/bash" > $TMPSCRIPT
echo "#\$" >> $TMPSCRIPT
echo "cd FOLDER;" >> $TMPSCRIPT
echo "./usearch7.0.1001_i86linux64 -makeudb_ublast $f -output $f.udb;" >> $TMPSCRIPT
qsub -l nodes=1 $TMPSCRIPT;
done
ls *.udb> nrudblist.txt


#if you need to convert your fastq files to fasta
for f in *.fastq
do 
perl ./bin/fastq2fa.pl $f > $f.converted.fa;
done


#blastn against reference genomes database
cat /taxonomy/dbs/*.fa > allgenomes.fa
ncbi-blast-2.2.26+/bin/makeblastdb -dbtype nucl -in /taxonomy/dbs/allgenomes.fa


##################
#first, blastn against reference databases
#blastn
for f in *.fa; 
do
TMPSCRIPT="$f.all.qsub"
echo "#\$ -S /bin/bash" > $TMPSCRIPT
echo "#\$" >> $TMPSCRIPT
echo "cd FOLDER;" >> $TMPSCRIPT
echo "./ncbi-blast-2.2.26+/bin/blastn -query $f -db /taxonomy/dbs/allgenomes.fa -outfmt 6 -num_threads 16 -evalue 1e-5 -max_target_seqs 10 -out $f.all.blast;" >> $TMPSCRIPT
qsub -l nodes=1:c16 $TMPSCRIPT
done

#get top hits
for f in *.fa; 
do
cat $f.all.blast | sort -k1,1 -k11g,11g | sort -u -k1,1| sort -r > $f.sortuniq.blast;
done

#now, let's retrieve all the sequences/genes/contigs that didn't map 
for f in *.fa
do
perl ./bin/notmatched.pl $f $f.sortuniq.blast
done

ls *.unmapped..txt> unmapped.txt


#blastx
TMPSCRIPT="usearch.qsub"
echo "#\$ -S /bin/bash" > $TMPSCRIPT
echo "#\$" >> $TMPSCRIPT
echo "cd FOLDER;" >> $TMPSCRIPT

for f in `cat ./blastdb/nrudblist.txt`
do
echo "./usearch7.0.1001_i86linux64 -ublast FOO --db ./blastdb/$f -blast6out FOO.$f.outfmt6.txt -evalue 0.01 -threads 4 -accel 0.5 -maxhits 2 --log FOO.$f.log;" >>$TMPSCRIPT
done

for f in `cat unmapped.txt`
do
sed s/FOO/$f/g usearch.qsub> $f.qsub
qsub -l nodes=1:g8 $f.qsub
done

#depending on your needs, you can then sort out and pull out only the top hit. 
TMPSCRIPT="uniq_blastx.swarm"
echo "#\$ -S /bin/bash" > $TMPSCRIPT
echo "#\$" >> $TMPSCRIPT
echo "cd FOLDER;" >> $TMPSCRIPT
for f in `cat unmapped.txt`
do
echo "cat $f.*.outfmt6.txt | sort -k1,1 -k11g,11g | sort -u -k1,1| sort -r > $f.sortuniq.blastx;" >> $TMPSCRIPT
done
qsub -l nodes=1:g8 uniq_blastx.swarm

#now, combine all blastn and blastx files to get a blastn or blastx file. 
for f in *.fa
do
cut -f 1,2 $f.sortuniq.blast > $f.blastnlist.txt
cut -f 1,2 $f.sortuniq.blastx > $f.blastxlist.txt
done

#get all the gis
cat *.outfmt6.txt | cut -f 2 | cut -f 2 -d '|' | sort | uniq > allgis.txt
perl ./bin/gitotaxid_2.pl allgis.txt ./blastdb/gi_taxid_prot.dmp #this requires >30GB RAM

#so now we have as our reference files to go into R allgis.txt.output.txt and blastlist.txt
#some processing of taxonomy strings
R CMD BATCH blastn_taxIDtotaxonomy.R
R CMD BATCH blastx_taxIDtotaxonomy.R

#now, combining all the taxonomies
for f in *.fa
do
perl ./bin/matchblastn_taxid_tofile.pl blastnlist_FINAL.txt $f*.sortuniq.blast
perl ./bin/matchblastx_taxid_tofile.pl blastxlist_FINAL.txt $f*.sortuniq.blastx
cat $f*.sortuniq.blast.output.txt $f*.sortuniq.blastx.output.txt > $f.blast.matched.txt
done

