#!/usr/bin/perl -w
use warnings;
use strict;
use FindBin;                # where was script installed?
use lib $FindBin::Bin;      # use that dir for libs, too
use Getopt::Long;
use Data::Dumper;
use File::Temp qw/ tempfile tempdir /;

my $ifh;
my $ofh;
my $opt = {format => 'fasta', logtemp => 'logXXXX', loglabel => $0.":".$$ };
if (-t STDIN && ! @ARGV){usage()};
GetOptions ($opt, 'in=s', 'format=s', 'out=s', 'nosingles', 'help', 'logtemp=s', 'loglabel=s' );
if ($opt->{help}){usage()};
my ($log, $logfile) = tempfile($opt->{logtemp}, UNLINK => 0 );

#Figure out whether we are using files or pipes for input
if ($opt->{in})
{open($ifh,"<".$opt->{in}) or die "can not open ".$opt->{in}."\n"}
elsif (! -t STDIN)
{$ifh=\*STDIN}
else
{usage()};
#now $ifh, is reading from a file or STDIN as needed 

#Output
if ($opt->{out})
{
    open($ofh,">".$opt->{out}) or die "can not open ".$opt->{out}."\n";
}
else
{$ofh=\*STDOUT};


#Parse stream
print $log join("\t",time(),$opt->{loglabel},"INFO","START",1)."\n";
my %dat; #store pair data
my $written=0;
my $l=0;
while (my $newlin = <$ifh>)
{
    chomp($newlin);
    if ($newlin=~/^\@/)
    {
	#header
	next;
    }
    my ($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual)=split(/\t/,$newlin);
    $l++;
    my $key1=$qname;
    my $bits=dec2bin($flag);
    my $key2;
    if ($bits=~/...01....../){$key2="1"}
    elsif ($bits=~/...10....../){$key2="2"}
    else
    {
	#unpaired
	$key2="1";
    }
    $dat{$key1}{$key2}{'s'}=$seq;
    $dat{$key1}{$key2}{'q'}=$qual;
    if (scalar(keys(%{$dat{$key1}}))==2)
    {
	#pair complete, write and delete from hash
	if ($opt->{format} eq 'table')
	{
	    print $ofh "$key1";
	    if (exists($dat{$key1}{'1'})){print $ofh "\t".join("\t",$dat{$key1}{'1'}{'s'},$dat{$key1}{'1'}{'q'})};
	    if (exists($dat{$key1}{'2'})){print $ofh "\t".join("\t",$dat{$key1}{'2'}{'s'},$dat{$key1}{'2'}{'q'})};
	    print $ofh "\n";
	}
	elsif ($opt->{format} eq 'fasta')
	{
	    if (exists($dat{$key1}{'1'})){print $ofh ">".$key1."/1 QUAL=".$dat{$key1}{'1'}{'q'}."\n".$dat{$key1}{'1'}{'s'}."\n"};
	    if (exists($dat{$key1}{'2'})){print $ofh ">".$key1."/2 QUAL=".$dat{$key1}{'2'}{'q'}."\n".$dat{$key1}{'2'}{'s'}."\n"};
	}
	elsif ($opt->{format} eq 'fastq')
	{
	    if (exists($dat{$key1}{'1'})){print $ofh "@".$key1."/1\n".$dat{$key1}{'1'}{'s'}."\n"."+".$key1."/1\n".$dat{$key1}{'1'}{'q'}."\n"};
	    if (exists($dat{$key1}{'2'})){print $ofh "@".$key1."/2\n".$dat{$key1}{'2'}{'s'}."\n"."+".$key1."/2\n".$dat{$key1}{'2'}{'q'}."\n"};
	}

	$written++;
	delete($dat{$key1});
    }
    if ($l%100000 == 0)
    {
	print STDERR "Reading $l... keeping track of ".scalar(keys(%dat))." pairs and $written pairs have been written\n";
    }
}
#clear out data hash;
my $singles=0;
foreach my $key1 (keys %dat)
{
    my $single_flag=0;
    if (scalar(keys(%{$dat{$key1}}))==1){$single_flag=1};
    if ($single_flag==1 && $opt->{nosingle}){next};
    
    if ($opt->{format} eq 'table')
    {
	print $ofh "$key1";
	if (exists($dat{$key1}{'1'})){print $ofh "\t".join("\t",$dat{$key1}{'1'}{'s'},$dat{$key1}{'1'}{'q'})};
	if (exists($dat{$key1}{'2'})){print $ofh "\t".join("\t",$dat{$key1}{'2'}{'s'},$dat{$key1}{'2'}{'q'})};
	print $ofh "\n";
    }
    elsif ($opt->{format} eq 'fasta')
    {
	if (exists($dat{$key1}{'1'})){print $ofh ">".$key1."/1 QUAL=".$dat{$key1}{'1'}{'q'}."\n".$dat{$key1}{'1'}{'s'}."\n"};
	if (exists($dat{$key1}{'2'})){print $ofh ">".$key1."/2 QUAL=".$dat{$key1}{'2'}{'q'}."\n".$dat{$key1}{'2'}{'s'}."\n"};
    }
    if ($single_flag==0){$written++} else {$singles++};
    delete($dat{$key1});
}

print STDERR "Finished $l... keeping track of ".scalar(keys(%dat))." pairs and $written pairs and $singles singles have been written\n";
print $log join("\t",time(),$opt->{loglabel},"INFO","BAM_ROWS",$l)."\n";
print $log join("\t",time(),$opt->{loglabel},"INFO","PAIRS",$written)."\n";
print $log join("\t",time(),$opt->{loglabel},"INFO","SINGLES",$singles)."\n";
print $log join("\t",time(),$opt->{loglabel},"INFO","END",1)."\n";
close($ifh);
close($log);
##################

sub usage
{
    print "bam2pairs - converts the stream of data from samtools \'view\' into a fasta for cross_match or table for bowtie\n";
    print "usage: perl -w bam2pairs.pl -in bamview.txt -out out.fas [-format table|fasta][-nosingles]\n";
    print "usage: samtools view -F 0x0200 file.bam | perl -w bam2pairs.pl -out out.fas -format fasta\n";
    print "note : mates are paired sequentially in fasta with quality string in defline\n";
    print "note : logging: -logtemp logXXXX -loglabel MyScript\n";
    exit(0);
}

sub dec2bin
{
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    my $pad="0" x (11-length($str));
    return($pad.$str);
}
