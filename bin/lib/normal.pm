=pod

=head1 normal
parameter:
	size    number of output
	mean    average num
	sd      sd of nums
	limit   only return nums that above limit(+) or below limit(-)
	type	1 for fixed number, variable whole length
Exmaple
&normal(10000,1200,100,"1000-",1)

=cut

package normal;
require Exporter;
use warnings;
use strict;

our @ISA=qw(Exporter);
our @EXPORT=qw(normal);
our $VERSION=1.0;

sub normal
{
	my ($size,$mean,$sd,$limit,$type)=@_;
 #   print "TEST:from $info mean $mean sd $sd\n";
	my @r;
	my $wl=$size*$mean;
	if(defined $limit)
	{
		my $lm=$limit;
        $lm=~s/\+//;
        $lm=~s/\-//;
        if(defined $type && $type eq 1)
		{
			for(my $i=0;$i<$size;$i++)
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
				$x=$mean+$sd*$x;
				$x=int($x+0.5);
				
   #             print "TEST\tlm $lm\tlimit $limit\n";
				if(($limit=~/\+/ && $x<$lm) || ($limit=~/\-/ && $x>$lm) || $x<0)
				{
					$i--;
					next;
				}
				else
				{
					push @r,$x;
				}
			}
		}
		else
		{
			if($limit=~/-/)
			{
				$limit=~s/\-//;
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
					$x=$mean+$sd*$x;
					$x=int($x+0.5);
					if($x<=$limit && $x>=0)
					{
						push @r,$x;
						$l+=$x;
					}
				}
			}
			else
			{
				$limit=~s/\+//;
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
					$x=$mean+$sd*$x;
					$x=int($x+0.5);
					if($x>=$limit && $x>=0)
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
		for(my $i=0;$i<$size;$i++)
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
			$x=$mean+$sd*$x;
			$x=int($x+0.5);
			push @r,$x;
		}
	}
	return @r;
}
1;
