=pod 

=head1 Quality
Used for producing quality score for reads.

FASTQ quality format:
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
 ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
 ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
 .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
 !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
 |                         |    |        |                              |                     |
33                        59   64       73                            104                   126
 0........................26...31.......40                                
                          -5....0........9.............................40 
                                0........9.............................40 
                                   3.....9.............................40 
 0........................26...31........41                               

  S - Sanger        Phred+33,  raw reads typically (0, 40)
  X - Solexa        Solexa+64, raw reads typically (-5, 40)
  I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
  J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
	  with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
	  (Note: See discussion above).
  L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

*above fastq information from wikipedia http://en.wikipedia.org/wiki/FASTQ_format

IN Psim, TWO TYPES OF QUALITY ARE employed: START WITH \"!\" AND START WITH \"\@\".

Quality score:
Illumina:	quality of reads
				first tenth--normal(int($b/10+0.5),$QualityExp-3,$QualityVar*10,$QualityExp."-",0);
				last nine tenth--normal($b-int($b/10+0.5),$QualityExp,$QualityVar,"40-",0);
			quality of site
				first ten sites--Q=0.8*x+Qmean-5
				other sites--Q=-a*x^2+Qmean+3.5
				a=(3.5*n-41)/(n*(n+1)*(2*n+1)/6-385)
454:		quality of reads
				first tenth--normal(int($b/10+0.5),$QualityExp-3,$QualityVar*10,$QualityExp."-",0);
				others--normal($b-int($b/10+0.5),$QualityExp,$QualityVar,"40-",0);
			quality of site
				Q=-ax^2+Qmean+3.5
				a=21/(n+1)/(2*n+1)
SOLiD:		quality of reads
				first tenth--normal(int($b/10+0.5),$QualityExp-3,$QualityVar*10,$QualityExp."-",0);
				others--normal($b-int($b/10+0.5),$QualityExp,$QualityVar,"40-",0);
			quality of site
				Q=-ax^2+Qmean+3.5
				a=21/(n+1)/(2*n+1)

example:
for illumina(fixed read length):
	&Quality($Command,\@quality,$readleng,$Qmean,$Qsd,$readsnum,$Qtype);
for solid(fixed read length):
	&Quality($Command,\@quality,$readleng,$Qmean,$Qsd,$readsnum);
for 454(variable read length):
	&Quality($Command,\@quality,\@readlength,$Qmean,$Qsd);

=cut

package Quality;
require Exporter;
use warnings;
use strict;
use normal;

our @ISA=qw(Exporter);
our @EXPORT=qw(Quality);
our @VERSION=1.0;

sub Quality
{
	my ($Type,$quality,$length,$Qmean,$Qsd,$readnum,$Qtype)=@_;
	if($Type=~/roche/)
	{
		$readnum=scalar(@$length);
	}
	my (@Qfront,@Qbehind,@Qall);
	my $tenth=int($readnum/10+0.5);
#		print "TESTQUALITY $tenth,$Qmean,$Qsd,1\n";
	@Qfront=normal($tenth,$Qmean-3,$Qsd*10,$Qmean."-","1");

	@Qbehind=normal(($readnum-$tenth),$Qmean,$Qsd,"40-","1");
	push @Qall,@Qfront;
	push @Qall,@Qbehind;

	if($Type=~/illumina/i)
	{
		my $a=aillumina($length);
		my $qtype;
		if($Qtype=~/!/)
		{
			$qtype=33;
		}
		else
		{
			$qtype=64;
		}
		foreach my $Q(@Qall)
		{
			my $qreads="";
			for my $i(1..10)
			{
				my $q=int(0.8*$i+$Q-5);
				$q=40 if($q>40);
				$q=0 if($q<0);
				$q+=$qtype;
				$qreads.=chr($q);
			}
			for my $j(11..$length)
			{
				my $q=int(-$a*$j*$j+$Q+3.5);
				$q=40 if($q>40);
				$q=0 if($q<0);
				$q+=$qtype;
				$qreads.=chr($q);
			}
			push @$quality,$qreads;
		}
	}	
	elsif($Type=~/solid/i)
	{
		my $a=asolid454($length);
		foreach my $Q(@Qall)
		{
			my $qreads="";
			for my $i(1..$length)
			{
				my $q=int(-$a*$i*$i+$Q+3.5);
				$q=40 if($q>40);
				$q=0 if($q<0);
				$qreads.="$q ";
			}
			push @$quality,$qreads;
		}
	}
	else
	{
		for my $k(0..$#$length)
		{
			my $qreads="";
			my $a=asolid454($$length[$k]);
			for my $i(1..$$length[$k])
			{
				my $q=int(-$a*$i*$i+$Qall[$k]+3.5);
				$q=40 if($q>40);
				$q=0 if($q<0);
				$qreads.="$q ";
			}
			push @$quality,$qreads;
		}
	}
}
1;

sub aillumina
{
#	(3.5*n-41)/(n*(n+1)*(2*n+1)/6-385)
	my $readlength=shift;
	my $a=(3.5*$readlength-41)/($readlength*($readlength+1)*(2*$readlength+1)/6-385);
	return $a;
}
sub asolid454
{
	my $readlength=shift;
	my $a=21/($readlength+1)/(2*$readlength+1);
	return $a;
}
