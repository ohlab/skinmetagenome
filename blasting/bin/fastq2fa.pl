#convert fastq to fasta files

#!/usr/bin/perl -w
use strict;

open(FQ,$ARGV[0]);
my $c=0;
while (my $d1=<FQ>)
{
    if ($d1=~/^\s*\n/){next}
    $c++;
    my $s=<FQ>;
    my $d2=<FQ>;
    my $q=<FQ>;
    #error check
    if ($d1 !~/^\@/ || $d2!~/^\+/ || $s !~/^[A-Z]+$/i)
    {
	die "Error $c\n".join("\n",$d1,$s,$d2,$q)."\n";
    }
    $d1=~s/^\@/>/;
    print $d1.$s;
}
