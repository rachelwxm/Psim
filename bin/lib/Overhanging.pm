=head1 Overhanging

Copyright 2014, Rachel Wu. All Rights Resered.  
This program is free software. You may copy or redistribute 
it under the same terms as Perl itself.  

=head2 Version

=head4 Author: Wu Xiaomeng(Rachel) rachelwu123@gmail.com

=head4 Version: 1.5

=head4 Update: 2014-04-21

=head2 USAGE

Used to generate overhanging and variant short DNA fragments


=head2 PARAMETERS

	library preparation type--1=single strain,2=double strain
	reference sequence(string reference)
	fragment start site and length array(array reference) 
	overhang length range of two strain
	sequence name
	strain lost rate of library type=1 or lost rate of one strain of library type=2
	mutation possibility along fragment site(or site region)(hash reference)
	    rate[$i]=$MutationSite{int(scalar(keys %MutationSite)*$i/$fragmentlength)}
	probabilities of different mutation types(hash reference)
	array that keep fragment sequences(array reference)
	array that keep information of each fragment(array reference)
	array that keep mutation information(array reference)
	efficiency of fill-in reaction

=head2 EXAMPLE

	#Library=1
	&Overhanging($Library,\$Sequence,\@StartLeng,$Mismatch,$SeqName,$lostSingle,$MutationSite,$MutationArray,\@Fragment,\@PaloName,\@MUTATION)
 	&Overhanging($Library,\$Sequence,\@StartLeng,$Mismatch,$SeqName,$lost,%MutationRate,%MutationArray,\@Fragment,\@PaloName,\@MUTATION,$Efficiency);
	#Library=2
	&Overhanging($Library,\$Sequence,\@StartLeng,$Mismatch,$SeqName,$DoublePercent,$MutationSite,$MutationArray,\@Fragment,\@PaloName,\@MUTATION)

=cut

package Overhanging;
require Exporter;
use warnings;
use strict;
use ReverseComplement;

our @ISA=qw(Exporter);
our @EXPORT=qw(Overhanging);
our $VERSION=1.5;
sub Overhanging
{
	my ($Library,$RefSeq,$SL,$Range,$Prefix,$lost,$MutationSite,$MutationArray,$fragment,$name,$mutationreport,$Efficiency)=@_;
	my $NumFrag=scalar(@$SL);
	my $usernum=scalar(keys %$MutationSite);
	my ($ll,$hl)=split /\_/,$Range;
	my $h;
#===========GENERATE FOUR TYPES OF OVERHANG==============
# four possibilities: 
#	s2<s1<e2<e2
#	s2<s1<e1<e2
#	s1<s2<e1<2e
#	s1<s2<e2<e1

#the first possibility:
#-----s1*********************e1
#s2*******************e2-------
	for my $i(0..(int($NumFrag/4)))
	{
		my ($s1,$l1)=split /\_/,$$SL[$i];
		next if($l1<5);
		my ($e1,$s2,$e2);
		$e1=$s1+$l1-1;
		$s2=$s1-int(rand($hl-$ll+1)+$ll);
		if($l1<(5+$ll))
		{
			$e2=$e1;
		}
		else
		{
			(($l1-5)>$hl)?($h=$hl):($h=($l1-5));
			$e2=$e1-int(rand($h-$ll+1)+$ll);
		}
		my $seq1=substr($$RefSeq,$s1,$l1);
		my $l2=$e2-$s2+1;
		my $seq2=&ReverseComplement(substr($$RefSeq,$s2,$l2));
		$l2=length($seq2);
		for my $j(0..($l1-1))
		{
			my $rate=$$MutationSite{int($usernum*$j/$l1)};
			if(rand(1)<$rate)
			{
				substr($seq1,$j,1)=&Mutation($MutationArray,$Prefix."\t$i-1\t$j\t",substr($seq1,$j,1),$mutationreport);
			}
		}
		for my $j(0..($l2-1))
		{
			my $rate=$$MutationSite{int($usernum*$j/$l2)};
			if(rand(1)<$rate)
			{
				substr($seq2,$j,1)=&Mutation($MutationArray,$Prefix."\t$i-2\t$j\t",substr($seq2,$j,1),$mutationreport);
			}
		}
		if($Library==1)	#single strain library preparation
		{
			if(rand(1)>$lost)
			{
				if(int(rand(2))==0)
				{
					$seq1=&ReverseComplement($seq1);
				}
				push @$fragment,$seq1;
				push @$name," $i\-1\ start\=$s1 end\=$e1 length\=$l1";
			}
			if(rand(1)>$lost)
			{
				if(int(rand(2))==0)
				{
					$seq2=&ReverseComplement($seq2);
				}
				push @$fragment,$seq2;
				push @$name," $i-2 start\=$s2 end=$e2 length=$l2";
			}
		}
		else	#double strain library preparation
		{
			if(int(rand(2))==0)
			{
				$seq1=&ReverseComplement($seq1);
			}
			$seq1=substr($seq1,0,($e2-$s1+1));
			push @$fragment,$seq1;
			push @$name," $i\-1\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
			if(rand(1)>$lost)
			{
				if(int(rand(2))==0)
				{
					$seq2=&ReverseComplement($seq2);
				}
				$seq2=substr($seq2,0,($e2-$s1));
				push @$name," $i\-2\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
				push @$fragment,$seq2;
			}
		}
	}

#the second possibility:
#-----s1**************e1-------
#s2**************************e2
	for my $ii((int($NumFrag/4)+1)..(int($NumFrag/2)))
	{
		my ($s1,$l1)=split /\_/,$$SL[$ii];
		next if($l1<5);
		my ($e1,$s2,$e2);
		$e1=$s1+$l1-1;
		$s2=$s1-int(rand($hl-$ll+1)+$ll);
		$e2=$e1+int(rand($hl-$ll+1)+$ll);
		my $seq1=substr($$RefSeq,$s1,($e1-$s1+1));
		my $seq2=ReverseComplement(substr($$RefSeq,$s2,($e2-$s2+1)));
		my $l2=length($seq2);
		for my $j(0..($l1-1))
		{
			my $rate=$$MutationSite{int($usernum*$j/$l1)};
			if(rand(1)<$rate)
			{
				substr($seq1,$j,1)=&Mutation($MutationArray,$Prefix."\t$ii-1\t$j\t",substr($seq1,$j,1),$mutationreport);
			}
		}
		for my $j(0..($l2-1))
		{
			my $rate=$$MutationSite{int($usernum*$j/$l2)};
			if(rand(1)<$rate)
			{
				substr($seq2,$j,1)=&Mutation($MutationArray,$Prefix."\t$ii-2\t$j\t",substr($seq2,$j,1),$mutationreport);
			}
		}
		if($Library==1)	#single strain library Prefixparation
		{
			if(rand(1)>$lost)
			{
				if(int(rand(2))==0)
				{
					$seq1=&ReverseComplement($seq1);
				}
				push @$fragment,$seq1;
				push @$name," $ii\-1\ start\=$s1 end\=$e1 length\=$l1";
			}
			if(rand(1)>$lost)
			{
				if(int(rand(2))==0)
				{
					$seq2=&ReverseComplement($seq2);
				}
				push @$fragment,$seq2;
				push @$name," $ii\-2\ start\=$s2 end\=$e2 length\=$l2";
			}
		}
		else
		{
			if(rand(1)<$Efficiency)
			{
				$seq1.=&ReverseComplement(substr($seq2,0,($e2-$e1)));
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq1=&ReverseComplement($seq1);
					}
					push @$name," $ii\-1\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
					push @$fragment,$seq1;
				}
				if(int(rand(2))==0)
				{
					$seq2=&ReverseComplement($seq2);
				}
				$seq2=substr($seq2,0,($e2-$s1+1));
				push @$name," $ii\-2\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
				push @$fragment,$seq2;
			}
			else
			{
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq1=&ReverseComplement($seq1);
					}
					push @$name," $ii\-1\ start\=$s1 end\=$e1 length\=$l1";
					push @$fragment,$seq1;
				}
				if(int(rand(2))==0)
				{
					$seq2=&ReverseComplement($seq2);
				}
				$seq2=substr($seq2,($e2-$e1),$l1);
				push @$name," $ii\-2\ start\=$s1 end\=$e1 length\=$l1";
				push @$fragment,$seq2;
			}
		}
	}
#the third possibility:
#s1*******************e1-------
#------s2********************e2
	for my $iii((int($NumFrag/2)+1)..(int($NumFrag/4*3)))
	{
		my ($s1,$l1)=split /\_/,$$SL[$iii];
		next if($l1<5);
		my ($e1,$s2,$e2);
		$e1=$s1+$l1-1;
		if($l1<(5+$ll))
		{
			$s2=$s1;
		}
		else
		{
			(($l1-5)>$hl)?($h=$hl):($h=($l1-5));
			$s2=$s1-int(rand($hl-$ll+1)+$ll);
		}
		$e2=$e1+int(rand($hl-$ll+1)+$ll);
		my $seq1=substr($$RefSeq,$s1,($e1-$s1+1));
		my $seq2=ReverseComplement(substr($$RefSeq,$s2,($e2-$s2+1)));
		my $l2=length($seq2);
		for my $j(0..($l1-1))
		{
			my $rate=$$MutationSite{int($usernum*$j/$l1)};
			if(rand(1)<$rate)
			{
				substr($seq1,$j,1)=&Mutation($MutationArray,$Prefix."\t$iii-1\t$j\t",substr($seq1,$j,1),$mutationreport);
			}
		}
		for my $j(0..($l2-1))
		{
			my $rate=$$MutationSite{int($usernum*$j/$l2)};
			if(rand(1)<$rate)
			{
				substr($seq2,$j,1)=&Mutation($MutationArray,$Prefix."\t$iii-2\t$j\t",substr($seq2,$j,1),$mutationreport);
			}
		}
		if($Library==1)	#single strain library preparation
		{
			if(rand(1)>$lost)
			{
				if(int(rand(2))==0)
				{
					$seq1=&ReverseComplement($seq1);
				}
				push @$fragment,$seq1;
				push @$name," $iii\-1\ start\=$s1 end\=$e1 length\=$l1";
			}
			if(rand(1)>$lost)
			{
				if(int(rand(2))==0)
				{
					$seq2=&ReverseComplement($seq2);
				}
				push @$fragment,$seq2;
				push @$name," $iii-2 start\=$s2 end=$e2 length=$l2";
			}
		}
		else
		{
			my $amp1=rand(1);
			my $amp2=rand(2);
			if($amp1<$Efficiency && $amp2<$Efficiency) #both amplify
			{
				$seq1.=&ReverseComplement(substr($seq2,0,($e2-$e1)));
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq1=&ReverseComplement($seq1);
					}
					push @$name," $iii\-1\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
					push @$fragment,$seq1;
				}
				$seq2.=&ReverseComplement(substr($seq1,0,($s1-$s2)));
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq2=&ReverseComplement($seq2);
					}
					push @$name," $iii\-2\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
					push @$fragment,$seq2;
				}
			}
			elsif($amp1<$Efficiency && $amp2>$Efficiency) #seq1 amplify, seq2 not amplify
			{
				$seq1=substr($seq1,($s2-$s1));
				$seq1.=&ReverseComplement(substr($seq2,0,($e2-$e1)));
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq1=&ReverseComplement($seq1);
					}
					push @$name," $iii\-1\ start\=$s2 end\=$e2 length\=$l2";
					push @$fragment,$seq1;
				}
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq2=&ReverseComplement($seq2);
					}
					push @$name," $iii\-2\ start\=$s2 end\=$e2 length\=$l2";
					push @$fragment,$seq2;
				}
			}
			elsif($amp1>$Efficiency && $amp2<$Efficiency) #seq2 amplify, seq1 not amplify
			{
				$seq2=substr($seq2,($e2-$e1));
				$seq2.=&ReverseComplement(substr($seq1,0,($s2-$s1)));
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq1=&ReverseComplement($seq1);
					}
					push @$name," $iii\-1\ start\=$s1 end\=$e1 length\=$l1";
					push @$fragment,$seq1;
				}
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq1=&ReverseComplement($seq1);
					}
					push @$name," $iii\-2\ start\=$s1 end\=$e1 length\=$l1";
					push @$fragment,$seq2;
				}
			}
			else	#neither amplify
			{
				$seq1=substr($seq1,($s2-$s1));
				$seq2=substr($seq2,($e2-$e1));
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq1=&ReverseComplement($seq1);
					}
					push @$name," $iii\-1\ start\=$s2 end\=$e1 length\=".($e1-$s2+1);
					push @$fragment,$seq1;
				}
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq2=&ReverseComplement($seq2);
					}
					push @$name," $iii\-2\ start\=$s2 end\=$e1 length\=".($e1-$s2+1);
					push @$fragment,$seq2;
				}
			}
		}
	}
#the fourth possibility:
#s1**************************e1
#------s2***********e2---------
	for my $iv((int($NumFrag/4*3)+1)..($NumFrag-1))
	{
		my ($s1,$l1)=split /\_/,$$SL[$iv];
		next if($l1<5);
		my ($e1,$s2,$e2);
		$e1=$s1+$l1-1;
		if($l1<(5+$ll))
		{
			$e2=$e1;
			$s2=$s1;
		}
		else
		{
			(($l1-5)>$hl)?($h=$hl):($h=($l1-5));
			$s2=$s1-int(rand($hl-$ll+1)+$ll);
			$e2=$e1-int(rand($h-$ll+1)+$ll);
			while($s2==$e2)
			{
				$e2=$e1-int(rand($h-$ll+1)+$ll);
			}
			($s2,$e2)=sort{$a<=>$b}($s2,$e2);
		}
		my $seq1=substr($$RefSeq,$s1,($e1-$s1+1));
		my $seq2=&ReverseComplement(substr($$RefSeq,$s2,($e2-$s2+1)));
		my $l2=length($seq2);
		for my $j(0..($l1-1))
		{
			my $rate=$$MutationSite{int($usernum*$j/$l1)};
			if(rand(1)<$rate)
			{
				substr($seq1,$j,1)=&Mutation($MutationArray,$Prefix."\t$iv-1\t$j\t",substr($seq1,$j,1),$mutationreport);
			}
		}
		for my $j(0..($l2-1))
		{
			my $rate=$$MutationSite{int($usernum*$j/$l2)};
			if(rand(1)<$rate)
			{
				substr($seq2,$j,1)=&Mutation($MutationArray,$Prefix."\t$iv-2\t$j\t",substr($seq2,$j,1),$mutationreport);
			}
		}
		if($Library==1)	#single strain library preparation
		{
			if(rand(1)>$lost)
			{
				if(int(rand(2))==0)
				{
					$seq1=&ReverseComplement($seq1);
				}
				push @$fragment,$seq1;
				push @$name," $iv-1 start=$s1 end=$e1 length=$l1";
			}
			if(rand(1)>$lost)
			{
				if(int(rand(2))==0)
				{
					$seq2=&ReverseComplement($seq2);
				}
				push @$fragment,$seq2;
				push @$name," $iv-2 start=$s2 end=$e2 length=$l2";
			}
		}
		else
		{
			if(rand(1)<$Efficiency)    #seq2 amplify
			{
				$seq2.=&ReverseComplement(substr($seq1,0,($s2-$s1)));
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq2=&ReverseComplement($seq2);
					}
					push @$name," $iv\-2\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
					push @$fragment,$seq2;
				}
				$seq1=substr($seq2,0,($e2-$s1+1));
				if(int(rand(2))==0)
				{
					$seq1=&ReverseComplement($seq1);
				}
				push @$name," $iv\-1\ start\=$s1 end\=$e2 length\=".($e2-$s1+1);
				push @$fragment,$seq1;
			}
			else		#seq2 not amplify
			{
				if(rand(1)>$lost)
				{
					if(int(rand(2))==0)
					{
						$seq2=&ReverseComplement($seq2);
					}
					push @$name," $iv\-2\ start\=$s2 end\=$e2 length\=$l2";
					push @$fragment,$seq2;
				}
				$seq1=substr($seq1,($s2-$s1),$l2);
				if(int(rand(2))==0)
				{
					$seq1=&ReverseComplement($seq1);
				}
				push @$name," $iv\-1 start\=$s2 end\=$e2 length\=$l2";
				push @$fragment,$seq1;
			}
		}
	}
}
1;
sub Mutation
{
	my ($array,$prefix,$char,$report)=@_;
	my $re;
	my @kind=keys %{$$array{$char}};
	if(scalar(@kind)==0)
	{
		$re=$char;	#not change
	}
	elsif(scalar(@kind)==1)
	{
		$re=$kind[0];
		push @$report,"$prefix$char\t$re";
	}
	elsif(scalar(@kind)==2)
	{
		if(rand(1)<$$array{$char}{$kind[0]})
		{
			$re=$kind[0];
			push @$report,"$prefix$char\t$re";
		}
		else
		{
			$re=$kind[1];
			push @$report,"$prefix$char\t$re";
		}
	}
	else
	{
		my $random=rand(1);
		if($random<$$array{$char}{$kind[0]})
		{
			$re=$kind[0];
			push @$report,"$prefix$char\t$re";
		}
		elsif($random<($$array{$char}{$kind[0]}+$$array{$char}{$kind[1]}))
		{
			$re=$kind[1];
			push @$report,"$prefix$char\t$re";
		}
		else
		{
			$re=$kind[2];
			push @$report,"$prefix$char\t$re";
		}
	}
	return $re;
}
1;
