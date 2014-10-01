#!/usr/bin/perl -w
use warnings;
use strict;
use FindBin;                # where was script installed?
use lib $FindBin::Bin;      # use that dir for libs, too
use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;
use File::Temp qw/ tempfile tempdir /;

my $format="fastq";

my $ifh;
my $ofh;
my $opt = {logtemp => 'logXXXX', loglabel => $0.":".$$ , matelengthdiff => 10};
if (-t STDIN && ! @ARGV){usage()};
GetOptions ($opt, 'in=s', 'out=s', 'longest=s', 'help', 'logtemp=s', 'loglabel=s', 'matelengthdiff=i' );
if ($opt->{help}){usage()};
my ($log, $logfile) = tempfile($opt->{logtemp}, UNLINK => 0 );

my $seqin;
#Figure out whether we are using files or pipes for input
if ($opt->{in})
{
    $seqin = Bio::SeqIO->new(
				-file     => $opt->{in},
				-format => 'fasta',
				);
}
elsif (! -t STDIN)
{
    $seqin = Bio::SeqIO->new(
				-fh     => \*STDIN,
				-format => 'fasta',
				);
}
else
{usage()};
#now $seqin, is reading from a file or STDIN as needed 

#Output
if ($opt->{out})
{
    open($ofh,">".$opt->{out}) or die "can not open ".$opt->{out}."\n";
}
else
{$ofh=\*STDOUT};

#Parse stream
print $log join("\t",time(),$opt->{loglabel},"INFO","START",1)."\n";
my %count; #hold accounting
my %read; #hold the read pair as it is found
my $n=0;
while (my $sq=$seqin->next_seq)
{
    $n++;
    my $s=$sq->seq;
    my $id=$sq->display_id;
    my $mate=1; 
    if ($id=~/\_1$/){$mate=1}
    elsif ($id=~/\_2$/){$mate=2}
    elsif ($id=~/\/1$/){$mate=1}
    elsif ($id=~/\/2$/){$mate=2}
    else {die "Unable to parse mate information from $id\n"};

    $id=~s/\_[12]$//; #trim off mate designator
    $id=~s/\/[12]$//; #trim off mate designator
    my $desc=$sq->desc;
    my $qual;
    if ($desc=~/QUAL\=(\S+)/){$qual=$1};

    #Trimming by keeping longest unmatched region
    my $match_coor=0;   #calculate for trimming qual
    my $match_length=0; #calculate for trimming qual
    if (exists($opt->{longest}))
    {
        my $match="";
        while ($s=~/([$opt->{longest}]+)/gi)
        {
            if (length($1)>length($match))
            {
                $match=$1;
                $match_length=length($match);
                $match_coor=pos($s)-$match_length;
            }
        }
	if ($match eq ""){next}; #skip completely masked sequences
        $s=$match;
	$qual=substr($qual,$match_coor,$match_length);
    }
    if (! $qual || length($qual) ne length($s))
    {
	print STDERR "Error while trimming quality score for $n sequence $id !\n";
	print STDERR "$s\n";
	print STDERR "$qual\n";
	die "fatal error\n";;
    }

    #data for this sequence stored in: $id,$s,$mate,$qual

    my @ml=sort(keys(%read));
    if ( scalar(keys(%read))==0 )
    {
	#store it and move on
	$read{$mate}{'name'}=$id;
	$read{$mate}{'seq'}=$s;
	$read{$mate}{'qua'}=$qual;
    }
    elsif ( scalar(keys(%read))>0 && $id eq $read{$ml[0]}{'name'})
    {
	#We have a pair of mates
	$read{$mate}{'name'}=$id;
	$read{$mate}{'seq'}=$s;
	$read{$mate}{'qua'}=$qual;
	my $line=format_read2(\%read);
	if ($line){print $ofh $line."\n"};
	undef(%read);
    }
    elsif ( scalar(keys(%read))>0 && $id ne $read{$ml[0]}{'name'})
    {
	#Current read does not match stored read, print single then store
	my $line=format_read2(\%read);
	if ($line){print $ofh $line."\n"};
	undef(%read);
	$read{$mate}{'name'}=$id;
	$read{$mate}{'seq'}=$s;
	$read{$mate}{'qua'}=$qual;
    }
    else {die "got lost while pairing reads\n"};
}

if (scalar(keys(%read))>0)
{
    #There is a final read pair or singlet in memory
    my $line=format_read2(\%read);
    if ($line){print $ofh $line."\n"};
    undef(%read);
}

close($ofh);
foreach my $k (sort keys %count)
{
    print $log join("\t",time(),$opt->{loglabel},"INFO",$k,$count{$k})."\n";    
}
print $log join("\t",time(),$opt->{loglabel},"INFO","END",1)."\n";

##################

sub usage
{
    print "mask_trimmer - trims crossmatch masked reads\n";
    print "usage: perl -w mask_trimmer -in in.fas -out out.table [-longest]\n";
    print "Expects pair sorted input fasta like:\n";
    print ' >C038EACXX:5:2108:03914:64159_1 QUAL=@@@DDDFD8FDHHE?ECAIIIIIIEGA?C9FBB?F?E...'."\n";
    print ' CCCCTGAGCTAGCCATGCTCTGACAGTCTCAGTTGCACACACGAGCCAGCAGAGGGCTGTCTCTCATACACGTC...'."\n";
    print ' >C038EACXX:5:2108:03914:64159_2 QUAL=#####################################...'."\n";
    print ' TCTCCCCGCCTCCCGACGCCCTTCGGAGATGTCTATAAGAGACAGCCCCTGCGCTAGCCATGCTCTGAAAGTCT...'."\n";
    print "-matelengthdiff specifies the maximum difference (default=10) in mate length\n";
    print "-longest ACTG extracts (trims) the longest stretch of sequence\n";
    print "note : logging: -logtemp logXXXX -loglabel MyScript\n";
    exit(0);
}

sub format_read2
{
    #a simpler read formatter, the original was trying to do a bunch of extra stuff
    my $read=$_[0];

    if (exists($read->{1}) && exists($read->{2}))
    {
	if ($read->{1}->{'name'} ne $read->{2}->{'name'}){die "Aack! reads with different names have been paired\n"};
	#Filter on length differences
	my $delta_length=abs(length($read->{1}->{'seq'})-length($read->{2}->{'seq'}));
	if ($delta_length>$opt->{matelengthdiff}){$count{'MATELENGTH_UNEQUAL_SKIP-'.$opt->{matelengthdiff}}++;return(0)};
	if ($format eq "fastq")
	{
	    my $str;
	    $str= "@".$read->{1}->{'name'}."/1\n".$read->{1}->{'seq'}."\n";
	    $str.="+".$read->{1}->{'name'}."/1\n".$read->{1}->{'qua'}."\n";
	    $str.="@".$read->{1}->{'name'}."/2\n".$read->{2}->{'seq'}."\n";
	    $str.="+".$read->{1}->{'name'}."/2\n".$read->{2}->{'qua'};
	    $count{'MATES'}++;
	    $count{'TOTAL_PASS_READS'}+=2;
	    return($str);
	}
    }
    elsif (exists($read->{1}))
    {
	if ($format eq "fastq")
	{
	    my $str;
	    $str= "@".$read->{1}->{'name'}."/1\n".$read->{1}->{'seq'}."\n";
	    $str.="+".$read->{1}->{'name'}."/1\n".$read->{1}->{'qua'};
	    $count{'SINGLE'}++;
	    $count{'TOTAL_PASS_READS'}++;
	    return($str);
	}
    }
    elsif (exists($read->{2}))
    {
	if ($format eq "fastq")
	{
	    my $str;
	    $str= "@".$read->{2}->{'name'}."/2\n".$read->{2}->{'seq'}."\n";
	    $str.="+".$read->{2}->{'name'}."/2\n".$read->{2}->{'qua'};
	    $count{'SINGLE'}++;
	    $count{'TOTAL_PASS_READS'}++;
	    return($str);
	}
    }
    else {return(0)}; #shouldn't happen...

}
