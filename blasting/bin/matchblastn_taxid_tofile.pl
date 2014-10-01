#!/usr/bin/perl -w
use strict;

my $ref = $ARGV[0];
my $file = $ARGV[1];
my %reffile = &build_seq_info_hash($ref);
# print_hash(%reffile); 
# exit;
	
	my $output_file = "$file.output.txt";

	open (OUTPUT, ">$output_file");
	open (FILE, $file);
			
	foreach (<FILE>) {
 	chomp;
		$_ =~ s/^\s+//;
		$_ =~ s/\s+$//;
		$_ =~ s/"//g;

my ($queryid, $subjectid, $perid, $alignlength,$mismatch,$nogap,$querystart,$queryend, $seqstart, $seqend, $eval, $bitscore) = split(/\t/);

	my ($taxid, $species) = split(/[|]/, $subjectid);
	$taxid =~ s/taxid://g;

	my $match = $reffile{$taxid} || $species;

 	  print OUTPUT "$queryid\t$subjectid\t$match\n";
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
		
		my($taxon,$designation,$full, $full2) = split(/\t/);
		
	$hash{$taxon} = $designation;

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
