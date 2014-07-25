#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my %Error;
$Error{"A"}{"0"}="T";
$Error{"A"}{"1"}="G";
$Error{"A"}{"2"}="C";
$Error{"a"}{"0"}="T";
$Error{"a"}{"1"}="G";
$Error{"a"}{"2"}="C";
$Error{"T"}{"0"}="A";
$Error{"T"}{"1"}="G";
$Error{"T"}{"2"}="C";
$Error{"t"}{"0"}="A";
$Error{"t"}{"1"}="G";
$Error{"t"}{"2"}="C";
$Error{"G"}{"0"}="A";
$Error{"G"}{"1"}="T";
$Error{"G"}{"2"}="C";
$Error{"g"}{"0"}="A";
$Error{"g"}{"1"}="T";
$Error{"g"}{"2"}="C";
$Error{"C"}{"0"}="A";
$Error{"C"}{"1"}="T";
$Error{"C"}{"2"}="G";
$Error{"c"}{"0"}="A";
$Error{"c"}{"1"}="T";
$Error{"c"}{"2"}="G";

my ($ref,$percent,$help,$rate);
GetOptions
(
	"f:s"	=>	\$ref,
	"p:s"	=>	\$percent,
	"h"		=>	\$help,
	"r:s"	=>	\$rate,
);

die 
"USAGE:    Used to randomly generate snp infomation with given reference sequence and SNP rate
Author:   Rachel Wu(rachelwu123\@gmail.com)
Version:  v1.0
Update:   2014-04-17

perl $0 [-f] [-p] [-r] [-h HELP]
    -f	reference sequence file
    -p	percent of SNP sites
    -r	rate of difference SNP types(eg \"0.7,0.3\" means 70% SNP sites have
        one SNP possiblity and 30% SNP sites have two SNP possibilities, 
        \"0.5,0.3,0.2\" means 50% SNP sites have just one SNP possiblity,
        30% SNP sites have two SNP possibilities and 20% SNP sites have
        three SNP possiblities)

EXAMPLE
    perl $0 -f ../example/example.fa -p 0.01 -r 0.7,0.3 > SNPinfo.txt
\n" if($help || !$ref || !$percent || !$rate);

my $sequence="";
my $name="";
open FILE,"<$ref";
while(<FILE>)
{
	if(/^>/)
	{
		if($sequence ne "")
		{
			&SNPGENERATE;
		}
		$name=$_;
		$name=~s/^>//;
		chomp $name;
		$sequence="";
	}
	else
	{
		chomp;
		$sequence.=$_;
	}
}
&SNPGENERATE;
close FILE;

sub SNPGENERATE
{
	my $RefLength=length($sequence);
	print "TEST: length of ref seuqence is $RefLength\n";
	my $SNPCount=int($percent*$RefLength+0.5);
	my %SNPInfo;	#key-site	value-base and possiblity
	
	my $Type=$rate=~tr/,/,/;
	if($Type==0)	#one SNP possiblity
	{
		for(1..$SNPCount)
		{
			my $Site=int(rand($RefLength));
			$Site=int(rand($RefLength)) while(exists $SNPInfo{$Site});

			my $Ini=substr($sequence,$Site,1);
			my $a=int(rand(3));
			while(!exists $Error{$Ini}{$a})
			{
				$Site=int(rand($RefLength));
				$Ini=substr($sequence,$Site,1);
			}
			my $SNPRate=(int(rand(9))+1)/10;
	$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\t$SNPRate";
#		$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\t$SNPRate\t$Ini";
		}
	}
	elsif($Type==1)	#two SNP possibilities
	{
		my @rates=split /\,/,$rate;
		my $single=int($SNPCount*$rates[0]+0.5);	#single SNP number
		for(1..$single)
		{
			my $Site=int(rand($RefLength));
			$Site=int(rand($RefLength)) while(exists $SNPInfo{$Site});

			my $Ini=substr($sequence,$Site,1);
			my $a=int(rand(3));
			while(!exists $Error{$Ini}{$a})
			{
				$Site=int(rand($RefLength));
				$Ini=substr($sequence,$Site,1);
			}
			my $SNPRate=(int(rand(9))+1)/10;
#	$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\t$SNPRate";
			$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\t$SNPRate\t$Ini";
		}
		for(($single+1)..$SNPCount)
		{
			my $Site=int(rand($RefLength));
			$Site=int(rand($RefLength)) while(exists $SNPInfo{$Site});

			my $Ini=substr($sequence,$Site,1);
			my $a=int(rand(3));
			my $b=int(rand(3));
			$b=int(rand(3)) while($b==$a);

			while(!exists $Error{$Ini}{$a})
			{
				$Site=int(rand($RefLength));
				$Ini=substr($sequence,$Site,1);
			}
			my $SNPRate1=(int(rand(9))+1)/10;
			my $SNPRate2=(int(rand(9))+1)/10;
			$SNPRate2=(int(rand(9))+1)/10 while($SNPRate2+$SNPRate1>1);

#$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\,$Error{$Ini}{$b}\t$SNPRate1\,$SNPRate2";
			$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\,$Error{$Ini}{$b}\t$SNPRate1\,$SNPRate2\t$Ini";
		}
	}
	elsif($Type==2)	#three SNP possiblities
	{
		my @rates=split /\,/,$rate;
		my $single=int($SNPCount*$rates[0]+0.5);	#single SNP number
			my $double=int($SNPCount*($rates[0]+$rates[1])+0.5);	#number tile double SNP

			for(1..$single)
			{
				my $Site=int(rand($RefLength));
				$Site=int(rand($RefLength)) while(exists $SNPInfo{$Site});

				my $Ini=substr($sequence,$Site,1);
				my $a=int(rand(3));
				while(!exists $Error{$Ini}{$a})
				{
					$Site=int(rand($RefLength));
					$Ini=substr($sequence,$Site,1);
				}
				my $SNPRate=(int(rand(9))+1)/10;
#				$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\t$SNPRate\t$Ini";
$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\t$SNPRate";
			}
		for(($single+1)..$double)
		{
			my $Site=int(rand($RefLength));
			$Site=int(rand($RefLength)) while(exists $SNPInfo{$Site});

			my $Ini=substr($sequence,$Site,1);
			my $a=int(rand(3));
			my $b=int(rand(3));
			$b=int(rand(3)) while($b==$a);

			while(!exists $Error{$Ini}{$a})
			{
				$Site=int(rand($RefLength));
				$Ini=substr($sequence,$Site,1);
			}
			my $SNPRate1=(int(rand(9))+1)/10;
			my $SNPRate2=(int(rand(9))+1)/10;
			$SNPRate2=(int(rand(9))+1)/10 while($SNPRate2+$SNPRate1>1);
$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\,$Error{$Ini}{$b}\t$SNPRate1\,$SNPRate2";
#			$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\,$Error{$Ini}{$b}\t$SNPRate1\,$SNPRate2\t$Ini";
		}
		for(($double+1)..$SNPCount)
		{
			my $Site=int(rand($RefLength));
			$Site=int(rand($RefLength)) while(exists $SNPInfo{$Site});

			my $Ini=substr($sequence,$Site,1);
			my $a=int(rand(3));
			my $b=int(rand(3));
			$b=int(rand(3)) while($b==$a);
			my $c;
			if(($a==0 && $b==1) || ($a==1 && $b==0)){$c=2}
			elsif(($a==0 && $b==2) || ($a==2 && $b==0)){$c=1}
			else{$c=0}

			while(!exists $Error{$Ini}{$a})
			{
				$Site=int(rand($RefLength));
				$Ini=substr($sequence,$Site,1);
			}
			my $SNPRate1=(int(rand(7))+1)/10;
			my $SNPRate2=(int(rand(10-10*$SNPRate1))+1)/10;
			my $SNPRate3=(int(rand(11-10*$SNPRate1-10*$SNPRate2))+1)/10;
$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\,$Error{$Ini}{$b}\t$SNPRate1\,$SNPRate2";
#			$SNPInfo{($Site+1)}="$Error{$Ini}{$a}\,$Error{$Ini}{$b}\t$SNPRate1\,$SNPRate2\t$Ini";
		}
	}
	else
	{
		die "Wrong Type parameter!\n";
	}
	foreach my $key(keys %SNPInfo)
	{
		print "$name\t$key\t$SNPInfo{$key}\n";
	}
}


close FILE;
