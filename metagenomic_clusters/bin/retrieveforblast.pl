#!/usr/bin/perl -w
use strict;

my $fasta = $ARGV[1];
my $ref = $ARGV[0];
my %reffile = &build_seq_info_hash($ref);
# print_hash(%reffile); 
# exit;

	my $output_file = "$ref.fasta";

	open (OUTPUT, ">$output_file");
	open (FASTA, $fasta);
			
	while (my $header = <FASTA>) {
		chomp($header);
		my $sequence = <FASTA>;
		chomp($sequence);

		my $header2 = $header;
		$header2 =~ s/>//g;

	my $match = $reffile{$header2} || "donget";

	if ($match eq "getit") {
 	print OUTPUT "$header\n$sequence\n";
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
		
	$hash{$_} = "getit";

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
