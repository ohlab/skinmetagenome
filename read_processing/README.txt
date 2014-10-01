This is a single-processor version of a metagenomic read-processing pipeline.

DISTRIBUTION:

This file              README.txt
Main shell script      read_proc.sh
Helper script          bin/bam2pairs.pl
Helper script          bin/mask_trimmer.pl
Nextera adapter fasta  bin/nextera_ends.fa

INSTALLATION:

1. Download and install Samtools, TrimBWAstyle.pl, PrinSeq and cross_match
2. Move the bam2pairs.pl and mask_trimmer.pl to a suitable location
3. Move the nextera_ends.fa file (or another adapter fasta file) to
   a suitable loction
4. Update the path variables in the head of the read_proc.sh script

RUN:

sh read_proc.sh input.bam

A number of log files and intermediate files are produced. These can safely
be compressed or deleted. Final output is in a file called:

input.bam.processed.fasta
