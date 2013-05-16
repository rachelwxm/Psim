=pod

=head1 MismatchDoubleStrain
parameter: 
	\@StartLeng
	mismatch length range(3_20)
	library type
	double strain percent
	lost rate of single strain
	\@DoubleSite,\@SingleSite
output:
	fragment

Example:
	#Library=1
	&MismatchDoubleStrain($Library,\$RefSeq,\@StartLeng,$Range,$Prefix,$lostSingle,$Mutation,\@Fragment)
	#Library=2
	&MismatchDoubleStrain($Library,\$RefSeq,\@StartLeng,$Range,$Prefix,$DoublePercent,$Mutation,\@Fragment)

=cut

package MismatchDoubleStrain;
require Exporter;
use warnings;
use strict;
use ReverseComplement;

our @ISA=qw(Exporter);
our @EXPORT=qw(MismatchDoubleStrain);
our $VERSION=1.0;

sub MismatchDoubleStrain
{
	my ($Library,$RefSeq,$SL,$Range,$Prefix,$lost,$MutationRate,$fragment,$name,$mutationreport)=@_;
	my $NumFrag=scalar(@$SL);
	my ($ll,$hl)=split /\_/,$Range;
	my $number=int(1/$MutationRate);
# four possibilities: 
#	s2<s1<e2<e2
#	s2<s1<e1<e2
#	s1<s2<e1<2e
#	s1<s2<e2<e1

#the first possibility:
#-----s1*********************e1
#s1*******************e2-------
#push @$name," $i\-1\ start=$s\ end=$e length=
	for my $i(0..(int($NumFrag/4)))
	{
		my ($s1,$l1)=split /\_/,$$SL[$i];
		my ($e1,$s2,$e2);
		$e1=$s1+$l1-1;
		$s2=$s1-int(rand($hl-$ll+1)+$ll);
		$e2=$e1-int(rand($hl-$ll+1)+$ll);
		my $seq1=substr($$RefSeq,$s1,($e1-$s1+1));
		my $seq2=ReverseComplement(substr($$RefSeq,$s2,($e2-$s2+1)));
		if($Library==1)	#single strain library preparation
		{
			if(rand(1)<$lost)
			{
				$seq1="";
			}
			else
			{
				for my $j(($e2-$s1+1)..($e1-$s1))
				{
					substr($seq1,$j,1)=&Mutation($number,$Prefix."\t$i\t$j\t",substr($seq1,$j,1),$mutationreport);
				}
				push @$name," $i\-1\ start\=$s1 end\=$e1 length\=".($e1-$s1+1);
			}
			if(rand(1)<$lost)
			{
				$seq2="";
			}
			else
			{
				for my $k(($e2-$s1+1)..($e2-$s2))
				{
					substr($seq2,$k,1)=&Mutation($number,$Prefix."\t$i\t$k\t",substr($seq2,$k,1),$mutationreport);
				}
				push @$name," $i-2 start\=$s2 end=$e2 length=".($e2-$s2+1);
			}
		}
		else
		{
			$seq1=substr($seq1,0,($e2-$s1+1));
			push @$name," $i\-1\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
			if(rand(1)>$lost)
			{
				$seq2="";
			}
			else
			{
				$seq2=ReverseComplement($seq1);
				push @$name," $i\-2\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
			}
		}
		push @$fragment,$seq1;
		push @$fragment,$seq2;
	}

#the second possibility:
#-----s1**************e1-------
#s1**************************e2
	for my $ii((int($NumFrag/4)+1)..(int($NumFrag/2)))
	{
		my ($s1,$l1)=split /\_/,$$SL[$ii];
		my ($e1,$s2,$e2);
		$e1=$s1+$l1-1;
		$s2=$s1-int(rand($hl-$ll+1)+$ll);
		$e2=$e1+int(rand($hl-$ll+1)+$ll);
		my $seq1=substr($$RefSeq,$s1,($e1-$s1+1));
		my $seq2=ReverseComplement(substr($$RefSeq,$s2,($e2-$s2+1)));
		if($Library==1)	#single strain library Prefixparation
		{
			if(rand(1)<$lost)
			{
				$seq1="";
			}
			else
			{
				push @$name," $ii\-1\ start\=$s1 end\=$e1 length\=".($e1-$s1+1);
			}
			if(rand(1)<$lost)
			{
				$seq2="";
			}
			else
			{
				for my $j(0..($e2-$e1-1))
				{
					substr($seq2,$j,1)=&Mutation($number,$Prefix."\t$ii\t$j\t",substr($seq2,$j,1),$mutationreport);
				}
				for my $k(($e2-$s1+1)..($e2-$s2))
				{
					substr($seq2,$k,1)=&Mutation($number,$Prefix."\t$ii\t$k\t",substr($seq2,$k,1),$mutationreport);
				}
				push @$name," $ii\-2\ start\=$s2 end\=$e2 length\=".($e2-$s2+1);
			}
		
		}
		else
		{
			for my $j(0..($e2-$e1-1))
			{
				substr($seq2,$j,1)=&Mutation($number,$Prefix."\t$ii\t$j\t",substr($seq2,$j,1),$mutationreport);
				
			}
			if(rand(1)>$lost)
			{
				$seq1="";
			}
			else
			{
				$seq1=ReverseComplement($seq2);
				push @$name," $ii\-1\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
			}
			push @$name," $ii\-2\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
		}
		push @$fragment,$seq1;
		push @$fragment,$seq2;
	}
#the third possibility:
#s1*******************e1-------
#------s1********************e2
	for my $iii((int($NumFrag/2)+1)..(int($NumFrag/4*3)))
	{
		my ($s1,$l1)=split /\_/,$$SL[$iii];
		my ($e1,$s2,$e2);
		$e1=$s1+$l1-1;
		$s2=$s1-int(rand($hl-$ll+1)+$ll);
		$e2=$e1+int(rand($hl-$ll+1)+$ll);
		my $seq1=substr($$RefSeq,$s1,($e1-$s1+1));
		my $seq2=ReverseComplement(substr($$RefSeq,$s2,($e2-$s2+1)));
		if($Library==1)	#single strain library preparation
		{
			if(rand(1)<$lost)
			{
				$seq1="";
			}
			else
			{
				for my $j(0..($s2-$s1-1))
				{
					substr($seq1,$j,1)=&Mutation($number,$Prefix."\t$iii\t$j\t",substr($seq1,$j,1));
				}
				push @$name," $iii\-1\ start\=$s1 end\=$e1 length\=".($e1-$s1+1);
			}
			if(rand(1)<$lost)
			{
				$seq2="";
			}
			else
			{
				for my $k(0..($e2-$e1-1))
				{
					substr($seq2,$k,1)=&Mutation($number,$Prefix."\t$iii\t$k\t",substr($seq2,$k,1));
				}
				push @$name," $iii\-2\ start\=$s1 end\=$e1 length\=".($e2-$s2+1);
			}
		
		}
		else
		{
			for my $j(0..($s2-$s1-1))
			{
				substr($seq1,$j,1)=&Mutation($number,$Prefix."\t$iii\t$j\t",substr($seq1,$j,1),$mutationreport);
			}
			for my $k(0..($e2-$e1-1))
			{
				substr($seq2,$k,1)=&Mutation($number,$Prefix."\t$iii\t$k\t",substr($seq2,$k,1),$mutationreport);
			}
			$seq1.=ReverseComplement(substr($seq2,0,($e2-$e1)));
			push @$name," $iii\-1\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
			if(rand(1)>$lost)
			{
				$seq2="";
			}
			else
			{
				$seq2=ReverseComplement($seq1);
				push @$name," $iii\-2\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
			}
		}
		push @$fragment,$seq1;
		push @$fragment,$seq2;
	}
#the fourth possibility:
#s1**************************e1
#------s1***********e2---------
	for my $iv((int($NumFrag/4*3)+1)..($NumFrag-1))
	{
		my ($s1,$l1)=split /\_/,$$SL[$iv];
		my ($e1,$s2,$e2);
		$e1=$s1+$l1-1;
		$s2=$s1-int(rand($hl-$ll+1)+$ll);
		$e2=$e1+int(rand($hl-$ll+1)+$ll);
		my $seq1=substr($$RefSeq,$s1,($e1-$s1+1));
		my $seq2=ReverseComplement(substr($$RefSeq,$s2,($e2-$s2+1)));
		if($Library==1)	#single strain library preparation
		{
			if(rand(1)<$lost)
			{
				$seq1="";
			}
			else
			{
				for my $j(0..($s2-$s1-1))
				{
					substr($seq1,$j,1)=&Mutation($number,$Prefix."\t$iv\t$j\t",substr($seq1,$j,1),$mutationreport);
				}
				for my $j(($e2-$s1+1)..($e1-$s1))
				{
					substr($seq1,$j,1)=&Mutation($number,$Prefix."\t$iv\t$j\t",substr($seq1,$j,1),$mutationreport);
				}
				push @$name," $iv\-1\ start\=$s1 end\=$e1 length\=".($e1-$s1+1);
			}
			if(rand(1)<$lost)
			{
				$seq2="";
			}
			else
			{
				push @$name," $iv\-2\ start\=$s2 end\=$e2 length\=".($e2-$s2+1);
			}
		}
		else
		{
			for my $j(0..($s2-$s1-1))
			{
				substr($seq1,$j,1)=&Mutation($number,$Prefix."\t$iv\t$j\t",substr($seq1,$j,1),$mutationreport);
			}
			$seq1=substr($seq1,0,($e2-$s1+1));
			push @$name," $iv\-1\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
			if(rand(1)>$lost)
			{
				$seq2="";
			}
			else
			{
				$seq2=ReverseComplement($seq1);
				push @$name," $iv\-2\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
			}
		}
		push @$fragment,$seq1;
		push @$fragment,$seq2;
	}
}


sub Mutation
{
	my ($num,$prefix,$char,$report)=@_;
	my $re;
	if((!int(rand($num))) && ($char=~/C/i))
	{
		$re="A";
push @$report,"$prefix$char\t$re";
	}
	else
	{
		$re=$char;
	}
	return $re;
}
1;
