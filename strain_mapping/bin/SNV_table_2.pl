#!/usr/bin/perl -w
#Sean Conlan

#Purpose: Runs the nucmer comparisons and processes the output to identify SNPs and unaligned regions
#Usage: perl SNV_table_2.pl genome1 genome2

use warnings;
use strict;
use FindBin;                 # where was script installed?
use lib $FindBin::Bin;      # use that dir for libs, too
#use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;
#use lib "/home/conlans/bin/SPC_util";
#use mymath;

#SETTINGS
my $show_unaligned=1; #set to 0 to suppress reporting reference positions that are not aligned
                      #set to 1 to write a row for every reference position that is unaligned

#check options and executable locations
if (@ARGV<2){usage()};
my ($db,$fas)=@ARGV;
my $out=`which nucmer`;
if ($out=~/no\snucmer/){die "are you aure you have nucmer installed?\n"};
if (! -e $db || ! -e $fas){usage()};

#reference needs to be repeat masked
#print "#----- repeat masking $db -----\n";
system("nucmer --maxmatch --nosimplify --prefix=REP $db $db 2>/dev/null");
system("show-coords -rT REP.delta > REP.table");
my $rep=get_repeat_hash("REP.table");

#get reference length
my $seqio_object = Bio::SeqIO->new(-file => $db);
my $ref   = $seqio_object->next_seq;
my $ref_size=$ref->length();
my $ref_seq=$ref->seq;

system("nucmer -p tmp $db $fas 2>/dev/null");
my $snps=`show-snps -ClrT -x 10 tmp.delta`;
system("show-snps -ClrT -x 10 tmp.delta > tmp.snps"); #save a copy for review
system("show-coords -rTl tmp.delta > tmp.table");

#Keep track of what part is aligned
my $coord=`show-coords -rT tmp.delta`;
my $align_mask="u" x $ref_size;
if (length($align_mask) ne length($ref_seq)){die "Error generating alignment mask"};
my $aligned_count=0;
my @coord=split(/\n/,$coord);
my $ln=0;
foreach my $c (@coord)
{
    $ln++;
    if ($ln<4 || $c eq "NUCMER" || $c eq "" || $c=~/^\[S1\].+\[TAGS\]$/){next};
    my ($s1,$e1,@rest)=split(/\t/,$c);
    if ($e1<$s1){die "Error pardsing coords\n"}; #expect reference orientation
    for (my $q=$s1;$q<=$e1;$q++)
    {
	if (substr($align_mask,$q-1,1) eq "u"){$aligned_count++};
	substr($align_mask,$q-1,1,"A");
    }
}
print STDERR "$aligned_count of $ref_size nt (".sprintf("%.1f", (100*($aligned_count/$ref_size)))."%) reference positions are in aligned regions\n";
print STDERR ($ref_size-$aligned_count)." of $ref_size nt (".sprintf("%.1f", (100*(($ref_size-$aligned_count)/$ref_size)))."%) reference positions are in unaligned regions\n";

#Step through SNVs and write table
print join("\t","ref_p","ref","qry","qry_p","ref_ctx...*..........","qry_ctx...*..........","is_OK","explanation","tag1","tag2")."\n";
my @s=split(/\n/,$snps);
for (my $i=4;$i<@s;$i++)
{
    my ($p1,$s1,$s2,$p2,$buff,$dist,$len_r,$len_q,$ctx_r,$ctx_q,$frm,@tags)=split(/\s+/,$s[$i]);

    if (substr($align_mask,$p1-1,1) ne "A"){die "Aaack! alignment mask does not agree with variant call\n"};
    #apply filters

    my %reject;
    if ($s1 eq "." || $s2 eq "."){$reject{'indel'}++};                                   #don't trust indels (454)
    if ($ctx_q=~/([^ACTG\.])/){$reject{'ambig-'.$1}++};                                  #don't trust SNV near ambiguous bases
    #good at finding SNVs very close together
    #if ($buff<=2){$reject{'near_snv'}++};                                                #don't trust clusters of SNVs
    if ($dist<=20){$reject{'near_end'}++};                                               #don't trust SNV near sequence ends
    if ($ctx_r=~/(AAAAA|TTTTT)/ || $ctx_q=~/(AAAAA|TTTTT)/){$reject{'homopolymer'}++};   #don't trust SNV near homopolymers (454)
    
    #searches in repeat hash to see if falls in one of the region identified as repeats
    my $rep_flag=0;                                                                      #don't trust SNV in repeat regions
    foreach my $rs (sort keys %{$rep})
    {
	if ($p1>=$rs && $p1<=$rep->{$rs} && $rep_flag != 1){$rep_flag=1;$reject{'in_repeat'}++};
    }

    my $is_OK="YES";
    if (scalar(keys(%reject)) != 0){$is_OK="no"};

    print join("\t",$p1,$s1,$s2,$p2,$ctx_r,$ctx_q,$is_OK,join(",",sort(keys(%reject))),$tags[scalar(@tags)-2],$tags[scalar(@tags)-1])."\n";
}

#Add unaligned positions in the reference
if ($show_unaligned==1)
{
    for (my $p=1;$p<=$ref_size;$p++)
    {
	if (substr($align_mask,$p-1,1) eq "A"){next};
	my $cl="----------"; #left context
	my $cr="----------"; #right context
	#lines below commented out to speed things up a bit
	if ($p>10){$cl=uc(substr($ref_seq,($p-11),10))};
	if ($p<$ref_size-10){$cr=uc(substr($ref_seq,$p,10))};
	print join("\t",$p,uc(substr($ref_seq,($p-1),1)),"-","-",$cl.uc(substr($ref_seq,($p-1),1)).$cr,"---------------------","no","unaligned")."\n";
    }
}

#########################################
# SUBROUTINES

sub usage
{
    #usage
    print STDERR "get_filtered_SNVs.pl reference.fasta genome.fasta\n";
    exit(0);
}

sub get_repeat_hash
{
    #reads output from show-coords -rT pKpQIL_REP.delta > pKpQIL_REP.table
    my %dat;
    my $r=$_[0];
    open(INF,"<$r") or die "can not open repeats\n";
    while (my $newlin=<INF>)
    {
	chomp($newlin);
	my @v=split(/\t/,$newlin);
	if (scalar(@v)<9 || $v[0]!~/^\d+$/ || $v[1]!~/^\d+$/ || $v[0] == $v[2]){next};
	$dat{min($v[0],$v[1])}=max($v[0],$v[1]);
	$dat{min($v[2],$v[3])}=max($v[2],$v[3]);
    }
    close(INF);
    #print Dumper(\%dat);exit(0);
    return(\%dat);
}

#moved out of mymath for portability
sub min
  {
    my @num=@_;
    my $minimum=$num[0];
    for (my $p=0;$p<@num;$p++)
      {
        if ($num[$p]<$minimum){$minimum=$num[$p]};
      }
    return($minimum);
  }

sub max
  {
    my @num=@_;
    my $maximum=$num[0];
    for (my $p=0;$p<@num;$p++)
      {
        if ($num[$p]>$maximum){$maximum=$num[$p]};
      }
    return($maximum);
  }
