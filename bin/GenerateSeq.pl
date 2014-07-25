#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Math::Random qw(random_poisson);
#use Random qw(random_poisson);

my ($n,$l,$auto,$pre);
my $help;
GetOptions(
	"n:s"	=>	\$n,
	"l:s"	=>	\$l,
	"pre:s"	=>	\$pre,
	"auto"	=>	\$auto,
	"h"		=>	\$help,
	"help"	=> 	\$help,
);
die
"USAGE:    Used to generate randomly sequence(s) with given length(s) or average legnth
Author:   Rachel Wu(rachelwu123\@gmail.com)
Version:  1.0
Updata:   2014-04-17

perl $0 [-n] [-l] [-pre] [-auto] > OUTFILE
    -n  NUM                number of sequence(default=1)
    -l  NUM or NUM,NUM...  length of each sequence, separate different values by \"\,\"
    -pre CHAR              prefix of your sequence(s)
    -auto                  randomly generate \"\-n\" sequences for Poisson distribution with average length of \"\-l\"

Example:
    #generate a random reference file length of 10000 bp
    perl $0 -l 10000 > MyNewRefFile.fa
    #generate two reference sequences with the length of 10kb and 15kb
    perl $0 -n 2 -l 10000,15000 -pre example > MyNewRefFile.fa 
    #generate several(assume 3) references with average length of 10kb
    perl $0 -n 3 -l 10000 -auto -pre example > My NewRefFile.fa
\n" if(!$l || $help);

$n||=1;
$pre||="by Rachel";
my @Length=split /\,/,$l;

if($auto)
{
	my $lamada=shift @Length;
	undef @Length;
	for(1..$n)
	{
		my $l=random_poisson(1,$lamada);
		push @Length,$l;
	}
}
my $numbertmp=scalar(@Length);

die "
Number of length NOT MATCH your \-n parameter!\n" if($numbertmp!=$n);

my %hash=
(
 "0"	=>	"A",
 "1"	=>	"T",
 "2"	=>	"G",
 "3"	=>	"C",
 );

for my $i(1..$n)
{
	my $seq="";
	for(1..$Length[$i-1])
	{
		$seq.=$hash{int(rand(4))};
	}
	print ">sequence$i"."of$pre, length\=$Length[$i-1] bp\n";
	print $seq,"\n";
}
