#This shell script expects a bam file aligned to HG19 or equivalent as input
#Multiple BAM files for a single sample should be merged using samtools like:
#  samtools merge [-nr] [-h inh.sam] <out.bam> <in1.bam> <in2.bam> [...]

#Requires: http://samtools.sourceforge.net/
#Requires: http://www.phrap.org/
#Recommended: TrimBWAstyle and PrinSeq (see below)

BAM=$1;                                   #input bam file
BINDIR="./bin"                            #where is bam2pair.pl and mask_trimmer.pl
CROSSMATCH="/bin/cross_match"   		  #recommend v 1.090518
ADAPTER="./bin/nextera_ends.fa"           #adapter sequences to trim
TRIMBWA="./bin/trimBWAstyle.pl"           #location of trim_BWAstyle.pl
PRINSEQ="./bin/prinseq-lite.pl"           #location of prinseq-lite.pl

echo "working on $BAM"
echo "... removing human alignable reads"
#Removing reads that:
# 1. Don't pass platform QC              (samtool flags)
# 2. Map uniquely to the human reference (samtool flags)
# 3. Map multiply to the human reference (grap commands)
#    Note: if your bam files do not contain the multiple mapping flags, remove
#          grep commands

samtools view -uF 0x0200 $BAM | samtools view -uf 0x0004 - | samtools view -f 0x0008 - | grep 'H0:i:0' | grep 'H1:i:0'| grep 'H2:i:0' | $BINDIR/bam2pairs.pl -logtemp logXXXX -out $BAM.no_human.paired.fasta 2> /dev/null

#The fasta file output in the previous step is sorted into read pairs and the
#quality scores are encoded in the defline.

echo "... trimming adapters"
#Mask adapter ends
#  Note that different versions of cross_match produces slightly different outputs
$CROSSMATCH $BAM.no_human.paired.fasta $ADAPTER -minmatch 10 -minscore 10 -screen > cm.log 2>&1

#Trim ends and convert to fastq
perl -w $BINDIR/mask_trimmer.pl -logtemp logXXXX -in $BAM.no_human.paired.fasta.screen -longest ACGT -out $BAM.no_human.trimmed.fastq

echo "... quality and length trimming"
#Quality trimming
#Final trimming and filtering is accomplished using two packages that you need to
#download, install and add to the path variables above:
# 1. http://wiki.bioinformatics.ucdavis.edu/index.php/TrimBWAstyle.pl
# 2. http://prinseq.sourceforge.net/
cat $BAM.no_human.trimmed.fastq | perl -w $TRIMBWA | perl -w $PRINSEQ -fastq stdin -out_format 1 -line_width 0 -out_good stdout -min_len 50 -out_bad null > $BAM.processed.fasta


#if you need to convert your fastq files to fasta
for f in *.fastq
do 
perl ./bin/fastq2fa.pl $f > $f.converted.fa;
done
