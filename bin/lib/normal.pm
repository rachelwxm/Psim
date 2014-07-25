=head1 normal

Copyright 2014, Rachel Wu. All Rights Resered.  
This program is free software. You may copy or redistribute 
it under the same terms as Perl itself.  

=head2 Version

=head4 Author: Wu Xiaomeng(Rachel) rachelwu123@gmail.com

=head4 Version: 1.5

=head4 Update: 2014-04-20

=head2 USAGE

Used to return an array of "size" number outcomes generated from Normal distribution with mean "mean", sd "sd" etc al. User can limit the minimum and/or maximum limit of outcomes that format as "minimum+", "maximum-" or "minimum+maximum-". For example, limit could be like "25+60-" what means the outcomes should larger than 25 and small than 60. Moreover, the distribution will not strictly meet the normal distribution. User can fix the number of outcomes or the whole length equal to "size"*"mean".

=head2 PARAMETERS

	size->number of outputs
	mean->average of Normal distribution
	sd->standard deviation
	limit->only return nums that above limit(+) or below limit(-)[optional]
	type->fixed type if set 'limit' [optional]
	    1 for fixed number, variable whole length
	    otherwise for fixed whole length

=head2 EXAMPLE

	#generate 10000 numbers of average of 1200 and sd of 100
	&normal("size",10000,"mean",1200,"sd",100)
	#generate 10000 numbers of average of 1200 and sd of 100, the numbers should smaller than 1500, fix 10000 outputs
	&normal("size",10000,"mean",1200,"sd",100,"limit","1500-","type",1)
	#generate 10000 numbers of average of 1200 and sd of 100, the numbers should smaller than 1500, fix the sum of outcomes equal to 10000*1200
	&normal("size",10000,"mean",1200,"sd",100,"limit","1500-","type",2)
	&normal("size",10000,"mean",1200,"sd",100,"limit","1500-")
	
	my @array2=&normal("size",$num2,"mean",$mean2,"sd",$sd2,"limit","$hash{fragmean}\+","type",$type);

=cut

package normal;
require Exporter;
use warnings;
use strict;

our @ISA=qw(Exporter);
our @EXPORT=qw(normal);
our $VERSION=1.5;

sub normal
{
	my %hash=@_;
	my @r;
	my $wl=$hash{size}*$hash{mean};
	if(defined $hash{limit})
	{
		my ($lm,$lm1,$lm2);
		if($hash{limit}=~/\+/ && $hash{limit}=~/\-/)
		{
			if($hash{limit}=~/(\d*\.?\d+)\+/)
			{
				$lm1=$1;
			}
			if($hash{limit}=~/(\d*\.?\d+)\-/)
			{
				$lm2=$1;
			}
		}
		else
		{
			$lm=$hash{limit};
			$lm=~s/\+//;
			$lm=~s/\-//;
		}			
		if(defined $hash{type} && $hash{type} eq 1)
		{
			for(my $i=0;$i<$hash{size};$i++)
			{
				my $yita1=rand(1);
				my $yita2=rand(1);
				while($yita1==0){$yita1=rand(1);}
				while($yita2==0){$yita2=rand(1);}
				my $x=sqrt(-2*log($yita1)/log(2.718281828459045))*sin(2*3.1415926*$yita2);
				$x=$hash{mean}+$hash{sd}*$x;
				$x=int($x+0.5);

				if($hash{limit}=~/\+/ && $hash{limit}=~/\-/)
				{
					if($x>=$lm1 && $x<$lm2)
					{
						push @r,$x;
					}
					else
					{
						$i--;
						next;
					}
				}
				elsif($hash{limit}=~/\+/ && !($hash{limit}=~/\-/))
				{
					if($x>=$lm)
					{
						push @r,$x;
					}
					else
					{
						$i--;
						next;
					}
				}
				elsif($hash{limit}=~/\-/)
				{
					if($x<=$lm)
					{
						push @r,$x;
					}
					else
					{
						$i--;
						next;
					}
				}
			}
		}
		else
		{
			if($hash{limit}=~/-/)
			{
				$hash{limit}=~s/\-//;
				my $l=0;
				while($l<$wl)
				{
					my $yita1=rand(1);
					my $yita2=rand(1);
					while($yita1==0)
					{   
						$yita1=rand(1);
					}
					while($yita2==0)
					{   
						$yita2=rand(1);
					}
					my $x=sqrt(-2*log($yita1)/log(2.718281828459045))*sin(2*3.1415926*$yita2);
					$x=$hash{mean}+$hash{sd}*$x;
					$x=int($x+0.5);
					if($x<=$hash{limit} && $x>=0)
					{
						push @r,$x;
						$l+=$x;
					}
				}
			}
			else
			{
				$hash{limit}=~s/\+//;
				my $l=0;
				while($l<$wl)
				{   
					my $yita1=rand(1);
					my $yita2=rand(1);
					while($yita1==0)
					{   
						$yita1=rand(1);
					}
					while($yita2==0)
					{   
						$yita2=rand(1);
					}
					my $x=sqrt(-2*log($yita1)/log(2.718281828459045))*sin(2*3.1415926*$yita2);
					$x=$hash{mean}+$hash{sd}*$x;
					$x=int($x+0.5);
					if($x>=$hash{limit} && $x>=0)
					{   
						push @r,$x;
						$l+=$x;
					}
				}
			}
		}
	}
	else
	{
		for(my $i=0;$i<$hash{size};$i++)
		{
			my $yita1=rand(1);
			my $yita2=rand(1);
			while($yita1==0)
			{
				$yita1=rand(1);
			}
			while($yita2==0)
			{
				$yita2=rand(1);
			}
			my $x=sqrt(-2*log($yita1)/log(2.718281828459045))*sin(2*3.1415926*$yita2);
			$x=$hash{mean}+$hash{sd}*$x;
			$x=int($x+0.5);
			push @r,$x;
		}
	}
	return @r;
}
1;
