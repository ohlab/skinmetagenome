This is a version that performs blastn against reference genome collection or blastx against the NCBI nr database using blastn and USEARCH, respectively and runs on a GNU/Linux parallel processing system.

DISTRIBUTION:

This file				README.txt
Main shell script		blasting.sh
Helper scripts			bin/blastx_taxIDtotaxonomy.R
Helper scripts			bin/blastn_taxIDtotaxonomy.R
						bin/matchblastx_taxid_tofile.pl
						bin/matchblastn_taxid_tofile.pl
						bin/matchtoblastx.pl
						bin/notmatched.pl
						bin/matchtax.pl
Reference files			blastdb/alltaxalist_reformatted.txt
						blastdb/taxIDmatch.txt

USAGE:

1. Download and install R, Perl, USEARCH v7.0 or greater. Also requires a legacy blast that contains fastacmd as well as BLAST 2.2.26 or thereabouts. 
2. Follow workflow provided in blasting.sh. Final output is a *.matched.txt file that contains your sequence name, the taxid or gi, and the taxonomy from that hit. Or blastn or blastx results separately are generated as intermediate files. 

