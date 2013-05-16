=pod

=head1 ReverseComplement

=cut

package ReverseComplement;
use Exporter;
use warnings;
use strict;

our @ISA=qw(Exporter);
our @EXPORT=qw(ReverseComplement);
our $VERSION=1.0;

sub ReverseComplement
{
	my $seq=shift;
	my $rev=reverse($seq);
	$rev=~tr/ATGCatgc/TACGtacg/;
	return $rev;
}
1;
