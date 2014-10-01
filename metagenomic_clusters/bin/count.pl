#!/usr/bin/perl -w
#
# Script to count the number of times a tag appears in a file
#
use strict;

# User variables - adjust these as necessary
my $infile = $ARGV[0];
my $column_to_count = $ARGV[1];



# Structure to hold/count the tags
my %tag_count;

open(IN, $infile);
while(<IN>) {
  chomp;
  
		$_ =~ s/^\s+//;
		$_ =~ s/\s+$//;  
		
  #print "$_\n";
  my(@fields) = split("\t", $_);
  my $tag = $fields[$column_to_count];
  #print "$tag\n";
  
  if(!$tag_count{$tag}) {
    $tag_count{$tag} = 1;
  }
  else { # means there was an existing tag
    $tag_count{$tag} = $tag_count{$tag} + 1;
  }
}

foreach my $tag(sort keys %tag_count) {
   {
    print "$tag\t$tag_count{$tag}\n";
  }
}
