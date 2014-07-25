=head1 ReverseComplement

Copyright 2014, Rachel Wu. All Rights Resered.  
This program is free software. You may copy or redistribute 
it under the same terms as Perl itself.  

=head2 Version

=head4 Author: Wu Xiaomeng(Rachel) rachelwu123@gmail.com

=head4 Version: 1.5

=head4 Update: 2014-04-20

=head2 USAGE

Used to translate nucleotide sequence to  reverse complement seqeunce.

=head2 EXAMPLE
	
	&ReverseComplement($sequence)

=cut

package ReverseComplement;
use Exporter;
use warnings;
use strict;

our @ISA=qw(Exporter);
our @EXPORT=qw(ReverseComplement);
our $VERSION=1.5;

sub ReverseComplement
{
	my $seq=shift;
	my $rev=reverse($seq);
	$rev=~tr/ATGCatgc/TACGtacg/;
	return $rev;
}
1;
