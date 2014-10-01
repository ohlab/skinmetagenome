#!/usr/bin/perl -w
use strict;


my $file = $ARGV[0];
my $ref = $ARGV[1];
my $type = $ARGV[2];

my %reffile = &build_seq_info_hash($ref);
# print_hash(%reffile); 
# exit;
	
	my $output_file = "$file.unmapped.$type.txt";
	my $output_file2 = "$file.allmapped.$type.txt";

	open (OUTPUT, ">$output_file");
	open (OUTPUT2, ">$output_file2");
	open (FILE, $file);
			
	while (my $header = <FILE>) {
		chomp($header);
		my $sequence = <FILE>;
		chomp($sequence);

		my ($a) = split (/\t/, $header);
		 $a =~ s/^.//s; 
	my $match = $reffile{$a} || "unmapped";

	if ($match eq "unmapped") {

 	  print OUTPUT ">$a\n$sequence\n";
 	  
 	  }

	print OUTPUT2 "$a\t$match\n";

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
		
	my ($a,$b,$c,$d,$e,$f,$g,$h,$i,$j,$k,$l) = split (/\t/);
		
	$hash{$a} = "mapped"

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
