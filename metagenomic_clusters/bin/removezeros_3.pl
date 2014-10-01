#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);

my $file = $ARGV[0];
my $output_file = "$file.nozero.txt";

	open (OUTPUT, ">$output_file");
	open (FILE, $file);

    print OUTPUT scalar <FILE>;	# Print the first line
	
	foreach (<FILE>) {
 	chomp;

		my @array = split (/\t/);
		my $identifier = shift @array;
 		my $sum = 0;
 		$sum = sum @array;

		if ($sum != 0) {
		print OUTPUT "$_\n";
	}
	}

exit;

