This is a version of the adaptive iterative assembly that runs on a GNU/Linux parallel processing system.

DISTRIBUTION:

This file              README.txt
Main shell script      assembling.sh
Helper script          bin/uninterleave.jar
Helper script          bin/velvet_stats.R


USAGE:

1. Download and install Velvet/Metavelvet, Bowtie2, Python2.7, khmer0.7 or greater, R, and seqtk
2. Follow workflow provided in assembly.sh

Final files are fasta files containing contigs for each assembly. 


NOTES:
Depending on the sample complexity and sequencing depth per sample, one may wish to use the open source digital normalization tool khmer on individual samples to reduce dataset size. For our samples, we found that assembly quality was very similar between khmer-normalized samples followed by Velvet assembly, and samples assembled directly using Velvet. 

We also recommend using a range of kmers for assembly, which will depend on your dataset. Typically we start with a range of 25-71. Reads are split into mate pairs and singletons, and then mate pairs are interleaved. Assembly quality is then evaluated with some basic quality metrics using some tools from the khmer package, and by mapping reads back to contigs. 
