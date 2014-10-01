#!/usr/bin/perl -w
use strict;

my $file = $ARGV[0];
my $ref = "/blasting/blastdbs/alltaxalist_reformatted.txt";
my $rHoH = build_hashohash($ref); 
# print_hash(%reffile); 
# exit;
	
	my $output_file = "$file.taxmatch.txt";

	open (OUTPUT, ">$output_file");
	open (FILE, $file);
			
	foreach (<FILE>) {
 	chomp;
		$_ =~ s/^\s+//;
		$_ =~ s/\s+$//;
		$_ =~ s/"//g;

		my ($a) = split (/\t/);
		
		my $species = $rHoH->{$a}->{'species'} || "FAIL";
		my $genus = $rHoH->{$a}->{'genus'} || "FAIL";
		my $family = $rHoH->{$a}->{'family'} || "FAIL";

 	  print OUTPUT "$species\t$genus\t$family\n";
	}



exit;


####################
sub build_hashohash{
	my $file = shift;
		my %HoH = ();

	open(IN, $file);
	foreach (<IN>) {
		chomp;
		$_ =~ s/^\s+//;
		$_ =~ s/\s+$//;
		$_ =~ s/"//g;
		
	my($match,$kingdom,$phylum,$class,$order,$family,$genus,$species,$new_species,$new_fam,$new_genus) = split(/\t/);
		
	$HoH{ $match }{ 'species' } = $new_species;
    $HoH{ $match }{ 'genus' } = $new_genus;	
    $HoH{ $match }{ 'family' } = $new_fam;	

					}
	close(IN);
	return \%HoH;
}
 

#subroutine of general use; prints a hash out
sub print_hash{
  my(%hash) = @_;
  foreach my $key(sort keys %hash) {
    print "'$key'\t'$hash{$key}'\n";
  }
}
