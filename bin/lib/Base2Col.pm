=pod
=head1 Base2Col
convert bases to color int reads

usage: base2col(\$seq)
=cut

package Base2Col;
require Exporter;
use warnings;
use strict;
our @ISA=qw(Exporter);
our @EXPORT=qw(Base2Col);
our @VERSION=1.0;
my @code=([0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]);
my @bases=qw/A C G T/;
my @othercode=(".",4,5,6);
my %colcode=();
my %decode=();
foreach my $i(0..3)
{
	foreach my $j(0..3)
	{
		$decode{$code[$i]->[$j]}->{$bases[$i]}=$bases[$j];
	}
}
foreach my $i(0..3)
{
	foreach my $j(0..3)
	{
		$colcode{"$bases[$i]$bases[$j]"}=$code[$i]->[$j];
	}

}
sub Base2Col
{
	my $reads=shift;
	my $Header=shift;
	my @bases=split '',uc($$reads);
	my $col_code = (defined $Header) ? $Header : "G";
	my $current_base = $col_code;
	my $code = '';
	for(my $i=0;$i<@bases;$i++)
	{
		$code = (exists $colcode{"$current_base$bases[$i]"})?$colcode{"$current_base$bases[$i]"}:$othercode[int(rand(@othercode))];
		$col_code .= $code;
		$current_base = $bases[$i];
	}
	$$reads=$col_code;
	1;
}
1;
