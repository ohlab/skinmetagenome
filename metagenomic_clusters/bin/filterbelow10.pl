#!/usr/bin/perl -w
use strict;
# use List::Util qw(sum);

my $file = $ARGV[0];
my $cutoff = $ARGV[1];
my $output_file = "$file.filt.$cutoff.txt";

	open (OUTPUT, ">$output_file");
	open (FILE, $file);

    print OUTPUT scalar <FILE>;	# Print the first line
	
	foreach (<FILE>) {
 	chomp;

		my @array = split (/\s/);
		my $identifier = shift @array;
# print "@array\n";
		
		my @zeros = grep (/^0$/, @array);
# print "@zeros\n";
		my $numzeros = @zeros;
# print "$numzeros\n";

# 		if ($numzeros < 29) { 
		if ($numzeros < $cutoff) { 
		print OUTPUT "$_\n";
		}
	}

exit;

