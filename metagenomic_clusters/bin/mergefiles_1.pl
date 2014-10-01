#make files the same length

#!/usr/bin/perl -w
use strict;

	
foreach my $file (@ARGV) {
	open (FILE, $file);

	my $output_file = "$file.output.txt";
	open (OUTPUT, ">$output_file");

my %reffile = &build_seq_info_hash($file);
# print_hash(%reffile); 
# exit;

my $ref = "allgenecalled.1000.lengths.txt";
open (REF, $ref);
	
	foreach (<REF>) {
 	chomp;
		$_ =~ s/^\s+//;
		$_ =~ s/\s+$//;

	my ($gene, $length) = split(/\t/);

	my $matched = $reffile{$gene} || "0";
#  	print OUTPUT "$_\t$matched\n";
 	print OUTPUT "$matched\n";
	}
}

exit;


####################
sub build_seq_info_hash{
	my $file = shift;
	my %hash;

	open(IN, $file);
	foreach (<IN>) {
		chomp;
		$_ =~ s/^\s+//;
		$_ =~ s/\s+$//;
		$_ =~ s/"//g;
		
		my($gene,$norm,$counts) = split(/\t/);
		
	$hash{$gene} = $norm;

					}
	close(IN);
	return %hash;
}
 

#subroutine of general use; prints a hash out
sub print_hash{
  my(%hash) = @_;
  foreach my $key(sort keys %hash) {
    print "'$key'\t'$hash{$key}'\n";
  }
}
