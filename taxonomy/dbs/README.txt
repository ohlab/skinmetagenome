Database is hosted at https://drive.google.com/file/d/0Bx2MDtKknM1eUjJETkJZRFYzaEE/edit?usp=sharing

Disclaimer (September 2014): We make no commitment beyond an amorphous endeavor to actively update any associated databases. However, most can be expanded in a fairly straightforward way, if you have specific references in mind. 

To add additional reference genomes, append a fasta formatted file to the end of any of the fasta files or create your own (e.g., append a fungal genome to all_fungiA.fa). Join contigs with a spacer of 100 N's. 

e.g.,
>Cool_fungus
GGGGGACTGNNNNN.......NNNNACTGGGGG

In the all_lengths.txt file, append a line containing the appropriate NCBI taxid number and the genome length in bp for your genome (not including the Ns). 

e.g., 
taxid:1000000000|Cool_fungus	10000000


Then remake the bowtie indices.  
./bin/bowtie2-2.2.3/bowtie2-build newfasta.fa newfasta;
