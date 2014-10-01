#!/usr/bin/perl -w
use strict;

my $ref = $ARGV[0];
my $list = $ARGV[1];
my %reffile = &build_seq_info_hash($ref);
print_hash(%reffile); 
exit;
	
	my $output_file = "$list.taxa.txt";

	open (OUTPUT, ">$output_file");
	open (LIST, $list);
			
	foreach (<LIST>) {
		chomp;

	my @array = split(/\t/); 
	my @idarray; 
	
	foreach my $id (@array) {
	my $match = $reffile{$id} || "nohit";
	
	push(@idarray, $match); 
	 	}
	 	
	print OUTPUT "@idarray\n"
	
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
		
		my($id, $gi) = split(/\t/);
		
	$hash{$id} = $gi;

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
