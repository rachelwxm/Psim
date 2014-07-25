=head1 PowerLaw

Copyright 2014, Rachel Wu. All Rights Resered.  
This program is free software. You may copy or redistribute 
it under the same terms as Perl itself.  

=head2 Version

=head4 Author: Wu Xiaomeng(Rachel) rachelwu123@gmail.com

=head4 Version: 1.5

=head4 Update: 2014-04-21

=head2 USAGE

E<verbar>  I<In statistics, a power law is a functional relationship between two quantities, where one quantity varies as a power of another. For instance, the number of cities having a certain population size is found to vary as a power of the size of the population. Empirical power-law distributions hold only approximately or over a limited range.>   

E<verbar> I< One attribute of power laws is their scale invariance. Given the probability density function f(x)=a*x^b, where xE<gt>0 and bE<lt>0.>  --wikipedia

Used to generate random number(s) fit power law distribution.

=head2 PARAMETERS

	number of outcomes
	parameter a, default 0.35
	power b, default -2
	upper limit

=head2 EXAMPLE

	#generate 10000 outcomes fit PowerLaw distribution with the paratemers of a=0.35 and b=-2, the outcomes should smaller than 72
	&PowerLaw(10000,0.35,-2,72)

=cut

package PowerLaw;
use Exporter;
use warnings;
use strict;

our @ISA=qw(Exporter);
our @EXPORT=qw(PowerLaw);
our $VERSION=1.5;

sub PowerLaw
{
	my ($num,$a,$b,$lim)=@_;
	die "wrong parameter b\n" if($b  && ($b>0));
	die "wrong parameter a\n" if($a  && $a>1);
	$a||=0.35;
	$b||=-2;
	$lim||=int(sqrt($a*10000));
	my @train;
	for(my $n=$lim;$n>0;$n--)
	{
		my $freq=int(10000*$a*$n**$b+0.5);
		$freq=1 if($freq < 1);
		for(1..$freq)
		{
			push @train,$n;
		}
	}
	my $tra=scalar(@train);
	my @return;
	if($num=~/r/)
	{
		$num=~s/[^0-9.]//g;
		my $cal=0;
		while($cal<$num)
		{
			my $i=int(rand($tra));
			my $this=$train[$i];
			push @return,$this;
			$cal+=$this;
		}
		my $over=pop @return;
		$cal-=$over;
		push @return,($num-$cal);
	}
	else
	{
		$num=~s/[^0-9.]//g;
		for(1..$num)
		{
			my $i=int(rand($tra));
			my $this=$train[$i];
			push @return,$this;
		}
	}
	return @return;
}
