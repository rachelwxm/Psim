=head1 PCRduplication

Copyright 2014, Rachel Wu. All Rights Resered.  
This program is free software. You may copy or redistribute 
it under the same terms as Perl itself.  

=head2 Version

=head4 Author: Wu Xiaomeng(Rachel) rachelwu123@gmail.com

=head4 Version: 1.5

=head4 Update: 2014-04-24

=head2 USAGE

Used to generate PCR duplication sequences. 
Amplification efficiency fit normal distribution. 

=head2 PARAMETERS
	
	DNA fragment array(array reference)
	Library DNA that have amplified(array reference)
	information of each fragment(array reference)
	information of each library DNA(array reference)
	average amplify multiples
	standard diviation of amplify multiples

=head2 EXAMPLE

	&PCRduplication(\@Fragment,\@Library,\@PaloName,\@LibraryName,$AmpMean,$AmpSD)

=cut

package PCRduplication;
use Exporter;
use warnings;
use strict;
use normal;

our @ISA=qw(Exporter);
our @EXPORT=qw(PCRduplication);
our $VERSION=1.5;

sub PCRduplication
{
	my ($Fragment,$Library,$name,$libraryname,$AmpMean,$AmpSD)=@_;
	my $NumFragment=scalar(@$Fragment);
	my @AmpNum=normal("size",$NumFragment,"mean",$AmpMean,"sd",$AmpSD);
	for my $i(0..($NumFragment-1))
	{
		for(1..$AmpNum[$i])
		{
			push @$Library,$$Fragment[$i];
			push @$libraryname,$$name[$i];
		}
	}
}
1;
		
