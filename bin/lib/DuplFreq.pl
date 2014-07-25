use warnings;
use strict;
use PowerLaw;
use Getopt::Long;

my ($num,$pa,$pb,$lim);
GetOptions(
	"num:s"	=>	\$num,
	"a:s"	=>	\$pa,
	"b:s"	=>	\$pb,
	"lim:s"	=>	\$lim,
);

die "
P(X=x)=a*x^b

perl $0 [-num] [-a] [-b] [-lim]
    num	number of digits(if user want to limit the whole reads num, please set this parameter as rNUM)
	a	parameter a (default 0.35)
	b	power b (default -2)
	lim	hight limit(default 59)

output:
	frequency of each duplication times (sorted)

examples:
	perl DuplFreq.pl -num 10000 > default-1
	perl DuplFreq.pl -num 10000 -lim 73 > lim-73
\n" if(!$num);



#my @return=&PowerLaw("number",$num,"a",$pa,"b",$pb,"limit",$lim);
my @return=&PowerLaw($num,$pa,$pb,$lim);

#foreach my $n(sort{$a<=>$b} @return)
#{
#	print "$n\n";
#}

#my @sort=sort{$a<=>$b} @return;

#my $s=1;
my $count=0;
#foreach(@sort)
#{
#	if($_ == $s)
#	{
#		$count++;
#	}
#	else
#	{
#		print "$s\t$count\n";
#		$s=$_;
#		$count=1;
#	}
#}
#print "$s\t$count\n";

foreach (@return)
{
	$count+=$_;
}
print "count is $count\n";
