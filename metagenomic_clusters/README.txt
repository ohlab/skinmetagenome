This is a version that clusters gene abundances across samples to produce metagenomic clusters, and runs on a GNU/Linux parallel processing system.

DISTRIBUTION:

This file				README.txt
Main shell script		clusters.sh
Helper scripts			bin/count.pl
						bin/countreads.py
						bin/filterbelow10.pl
						bin/gitotaxid_2.pl
						bin/matchtax_tocatalog.pl
						bin/mergefiles_1.pl
						bin/normalize_counts.py
						bin/removezeros_3.pl
						bin/retrieveforblast.pl
						bin/taxonomy.R

USAGE:

1. Download and install Bowtie2, Python2.7, R, USEARCH7.0 or greater, and MCL. Also requires /blasting
2. Follow workflow provided in clusters.sh

Final output here is a list of the clusters with a taxonomy, size, concordance with LCA, with large number of intermediate files with clustering information.  
species.txt
genus.txt
family.txt
