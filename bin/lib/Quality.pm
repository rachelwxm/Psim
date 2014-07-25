=head1 Quality

Copyright 2014, Rachel Wu. All Rights Resered.  
This program is free software. You may copy or redistribute 
it under the same terms as Perl itself.  

=head2 Version

=head4 Author: Wu Xiaomeng(Rachel) rachelwu123@gmail.com

=head4 Version: 1.5

=head4 Update: 2014-04-24

=head2 USAGE

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

IN Psim, TWO TYPES OF QUALITY ARE employed: START WITH "!" AND START WITH "@".

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

=head2 PARAMETERS

	sequencing platform type
	quality array(array reference)
	fixed read length(string) or read length array(array reference)
	average quality score
	quality score standard deviation
	reads num[for Illumina and SOLiD]
	FASTQ quality format(! or @)[for Illumina]

=head2 EXAMPLE

	#for illumina(fixed read length):
	&Quality($Command,\@quality,$readleng,$Qmean,$Qsd,$readsnum,$Qtype);
	#for solid(fixed read length):
	&Quality($Command,\@quality,$readleng,$Qmean,$Qsd,$readsnum);
	#for 454(variable read length):
	&Quality($Command,\@quality,\@readlength,$Qmean,$Qsd);

=cut

package Quality;
require Exporter;
use warnings;
use strict;
use normal;

our @ISA=qw(Exporter);
our @EXPORT=qw(Quality);
our @VERSION=1.5;

sub Quality
{
	my %hash=@_;
#	my ($Type,$quality,$length,$Qmean,$Qsd,$readnum,$Qtype)=@_;
	if($hash{cmd}=~/roche/)
	{
		$hash{num}=scalar(@{$hash{readleng}});
	}
	my (@Qfront,@Qbehind,@Qall);
	my $tenth=int($hash{num}/10+0.5);
#$print "$tenth\t$hash{mean}\t$hash{sd}\n";
#	print "\"size\"\,$tenth,\"mean\",$hash{mean}-3,\"sd\",$hash{sd}*10,\"limit\",$hash{mean}.\"-\",\"type\",\"1\",\n";
	@Qfront=normal("size",$tenth,"mean",$hash{mean}-3,"sd",$hash{sd}*10,"limit",$hash{mean}."-0+","type","1");
#print "@Qfront\n";

	@Qbehind=normal("size",($hash{num}-$tenth),"mean",$hash{mean},"sd",$hash{sd},"limit","40-0+","type","1");
	push @Qall,@Qfront;
	push @Qall,@Qbehind;

	if($hash{cmd}=~/illumina/i)
	{
		my $a=aillumina($hash{readleng});
		my $qtype;
		if($hash{type}=~/!/)
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
			for my $j(11..$hash{readleng})
			{
				my $q=int(-$a*$j*$j+$Q+3.5);
				$q=40 if($q>40);
				$q=0 if($q<0);
				$q+=$qtype;
				$qreads.=chr($q);
			}
			push @{$hash{array}},$qreads;
		}
	}	
	elsif($hash{cmd}=~/solid/i)
	{
		my $a=asolid454($hash{readleng});
		foreach my $Q(@Qall)
		{
			my $qreads="";
			for my $i(1..$hash{readleng})
			{
				my $q=int(-$a*$i*$i+$Q+3.5);
				$q=40 if($q>40);
				$q=0 if($q<0);
				$qreads.="$q ";
			}
			push @{$hash{array}},$qreads;
		}
	}
	else
	{
		for my $k(0..$#{$hash{readleng}})
		{
			my $qreads="";
			my $a=asolid454(${$hash{readleng}}[$k]);
			for my $i(1..${$hash{readleng}}[$k])
			{
				my $q=int(-$a*$i*$i+$Qall[$k]+3.5);
				$q=40 if($q>40);
				$q=0 if($q<0);
				$qreads.="$q ";
			}
			push @{$hash{array}},$qreads;
		}
	}
}
1;

sub aillumina
{
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
