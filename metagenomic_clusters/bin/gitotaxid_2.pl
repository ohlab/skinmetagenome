#!/usr/bin/perl -w
use strict;

my $file = $ARGV[0];
my $ref = $ARGV[1];
my %reffile = &build_seq_info_hash($ref);
# print_hash(%reffile); 
# exit;
	
	my $output_file = "$file.output.txt";

	open (OUTPUT, ">$output_file");
	open (FILE, $file);
			
	foreach (<FILE>) {
 	chomp;

	my ($gi) = split(/\t/);
	my $match = $reffile{$gi} || "1";

  	  print OUTPUT "$gi\t$match\n";
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
		
my ($gi, $taxid) = split(/\t/, $_);
		
	$hash{$gi} = $taxid;

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
