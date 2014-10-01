#Requires Bowtie2, Python2.7, R, USEARCH7.0 or greater, and MCL
#Designed for a GNU/Linux parallel processing system, e.g., here, 16 processors and at least 72 GB of RAM. 
#binaries are in your own defined directory


#uses the gene catalog from "genecatalog" built from contigs > 1000bp and containing contigs <1000bp. 
allgenecalled.1000.centroids.fa
allgenecalled.1000.lengths.txt


####Determine gene abundances for each sample. 
#Make a bowtie index from the gene catalog
./bowtie2-2.2.3/bowtie2-build allgenecalled.1000.centroids.nolines.fa allgenecalled.1000.centroids.nolines

nano bowtie.qsub
#\$ -S /bin/bash
#\$
cd FOLDER
./bowtie2-2.2.3/bowtie2 -p 16 -x allgenecalled.1000.centroids.nolines -q -U FOO --sensitive -k 1 -S FOO.1000.sensk1.sam > FOO.1000.sensk1.stats.txt 2>&1;
python2.7 ./bin/countreads.py FOO.1000.sensk1.sam FOO.1000.origsamcalc.txt;
python2.7 ./bin/normalize_counts.py FOO.1000.origsamcalc.txt FOO.1000.origsamcalc.normalized.txt allgenecalled.1000.lengths.txt;
perl ./bin/mergefiles_1.pl FOO.1000.origsamcalc.txt.normalized.txt
#end

for f in `cat listoffiles.txt`
do
sed s/FOO/$f/g bowtie.qsub> $f.bowtie.qsub
qsub -l nodes=1:c16 $f.bowtie.qsub
done

#now combine all the files into one
cut -f 1 allgenecalled.1000.lengths.txt > allgenecalledlist.txt
cat header.txt > dataframe.k1.txt	#header.txt is a tab delimited file with a first column called "gene" then your original sample names
paste allgenecalledlist.txt *.1000.origsamcalc.normalized.txt.output.txt >> dataframe.k1.txt

#a little cleanup by removing genes that have no hits
perl ./bin/removezeros_3.pl dataframe.k1.txt


###final file of gene abundances for clustering is 
dataframe.k1.txt.nozero.txt

#To reduce compuational burden, we split our samples by site microenvironment prior to clustering and required that a gene be present in at least x% of samples. This is an example of clustering across all samples but you can split the data in any way. 

perl ./bin/filterbelow10.pl dataframe.k1.txt.nozero.txt 227	#227 is the number of samples that the gene must occur in. Exclude genes which occur in fewer than ~20% of samples (starting with 283). 
sed '1d' dataframe.k1.227.renamed.txt | sort | uniq >dataframe.k1.227.renamed.txt.temp.txt;
cat dataframe.k1.227.renamed.txt.header.txt dataframe.k1.227.renamed.txt.temp.txt > dataframe.k1.227.renamed.txt.uniq.txt;

#clustering is for 85% correlation across samples and inflation factor -I is 2. This step can be very memory intensive 
time mcxarray -data dataframe.k1.227.renamed.txt.uniq.txt -t 32 -skipr 1 -skipc 1 -o dataframe.k1.227.renamed.txt.uniq.txt.mci -write-tab dataframe.k1.227.renamed.txt.uniq.txt.tab --spearman -co 0.85 -tf 'abs(),add(-0.85)' --write-binary;
time mcl dataframe.k1.227.renamed.txt.uniq.txt.mci -I 2 --d -te 32;
time mcxdump -icl out.dataframe.k1.227.renamed.txt.uniq.txt.mci.I20 -tabr dataframe.k1.227.renamed.txt.uniq.txt.tab -o dump.out.dataframe.k1.227.mci.85.I20


#####output file is:
mv dump.out.dataframe.k1.227.mci.85.I20 dataframe.k1.227.85.I20

#now, with the clusters, let's get a consensus taxonomy. Considering clusters of >10 genes: 
#1. Get the top X number of clusters. If number of fields > 9, retrieve. Convert to list.
awk 'NF > 9' dataframe.k1.227.85.I20 > dataframe.k1.227.85.I20.list
tr '\t' '\n' < dataframe.k1.227.85.I20.list > dataframe.k1.227.85.I20.headers

#2. Then, match against allgenecalled.1000.centroids.nolines.fa
#3. If match, print fasta sequence to a new file. Create a file for each cluster.
time perl ./bin/retrieveforblast.pl dataframe.k1.227.85.I20.headers allgenecalled.1000.centroids.nolines.fa

#4. Then, blast each group against nr, then take the LCA for at least 50% of the genes.
nano usearch.qsub 
for f in `cat /nr/blastdb/nrudblist.txt`
do
TMPSCRIPT="$f.usearch.qsub"
echo "#\$ -S /bin/bash" > $TMPSCRIPT
echo "#\$" >> $TMPSCRIPT
echo "cd FOLDER;" >> $TMPSCRIPT
echo "./usearch7.0.1001_i86linux64 -ublast dataframe.k1.227.85.I20.headers.fasta --db /nr/blastdb/$f -blast6out dataframe.k1.227.85.I20.headers.fasta.$f.outfmt6.txt -evalue 0.01 -threads 4 -accel 0.5 -maxhits 2 --log dataframe.k1.227.85.I20.headers.fasta.$f.log;" >>$TMPSCRIPT
qsub -l nodes=1:g8 $TMPSCRIPT
done

cat dataframe.k1.227.85.I20.headers.fasta.nr.*.udb.outfmt6.txt > dataframe.k1.227.85.I20.headers.fasta.blastx
cat dataframe.k1.227.85.I20.headers.fasta.blastx | sort -k1,1 -k11g,11g  | sort -u -k1,1| sort -r | cut -f 1,2> dataframe.k1.227.85.I20.headers.fasta.sortuniq.blastx

#now, let's get all the gis to get some taxonomy information
cat *.outfmt6.txt | cut -f 2 | cut -f 2 -d '|' | sort | uniq > allgis.txt
perl /blasting/bin/gitotaxid_2.pl allgis.txt /blasting/blastdb/gi_taxid_prot.dmp #this requires >30GB RAM
R CMD BATCH blastx_taxIDtotaxonomy.R
perl /blasting/blastdb/matchblastx_taxid_tofile.pl blastxlist_FINAL.txt dump.dataframek1_90_I2.sortuniq.blastx

#now, let's match back to the clusters
perl matchtoblastx.pl dump.dataframek1_90_I2.sortuniq.blastx.output.txt dump.dataframek1_90_I2.list

#5. For each cluster, get the LCA as a species, genus, family, order, class, phylum, no rank. 
mkdir I2_90_dataframek1 
split -l 1 -d -a 4 dump.dataframek1_90_I2.list.taxa.txt I2_90_dataframek1/line.dataframek1_90_I2.
cd I2_90_dataframek1

#convert tabs to newlines in the output files
for f in line*
do 
tr '\t' '\n' < $f > temp.$f;
mv temp.$f $f
done

#so, for each file, we want to ask, what's the consensus taxonomy? Let's put a 50% concordance to species, genus, and family level. If can't get at least 50% at family level discard. Then group things that are the same at the species level (let other levels remain unique).
cat line* |sort | uniq > ../taxalist.txt
cat *taxalist.txt | sort | uniq > alltaxalist.txt

for f in line*[0-9]
do
perl ./bin/matchtax_tocatalog.pl $f
perl ./bin/count.pl $f.taxmatch.txt 0 > $f.species.txt
perl ./bin/count.pl $f.taxmatch.txt 1 > $f.genus.txt
perl ./bin/count.pl $f.taxmatch.txt 2 > $f.family.txt
done

R CMD BATCH taxonomy.R &

#the final output here is a list of the clusters with a taxonomy, size, concordance with LCA. 
species.txt
genus.txt
family.txt




