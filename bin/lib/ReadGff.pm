=head1 READGFF

Copyright 2014, Rachel Wu. All Rights Resered.  
This program is free software. You may copy or redistribute 
it under the same terms as Perl itself.  

=head2 Version

=head4 Author: Wu Xiaomeng(Rachel) rachelwu123@gmail.com

=head4 Version: 1.5

=head4 Update: 2014-04-20

=head2 USAGE

read gff file and get exon/cdna region and gene info for Tag-Seq or RNA-Seq

=head2 PARAMETERS

	gff->Gff file name
	info->gff info hash(hash reference)

=head2 EXAMPLE

	&ReadGff("gff",$GFF,"info",\%DigestRegion);

=cut

package ReadGff;
require Exporter;
use warnings;
use strict;

our @ISA=qw(Exporter);
our @EXPORT=qw(ReadGff);
our $VERSION=1.5;

sub ReadGff
{
	my %hash=@_;
	open (IN,"<$hash{gff}") || die "Can't open $hash{gff}!\n";
	my $Exon="";
	my ($Chr,$Symbol,$Start,$End,$Strand,$Id);
	while(<IN>)
	{
		next if(/^\#/);
		chomp;
		my @t=split /\t/;
		if($t[2]=~/transcript/)
		{
			if($Exon ne "")
			{
				${$hash{info}}{$Chr}{$Id}=$Exon;
			}
			($Chr,$Symbol,$Start,$End,$Strand)=@t[0,2,3,4,6];
			$Exon="$Strand\:$Start\,$End";
			$Id=$1 if($t[8]=~/transcript\ \"(.*)\"/);
		}
		elsif($t[2]=~/exon/)
		{
			$Exon.="\;$t[3]\,$t[4]";
		}
	}
	${$hash{info}}{$Chr}{$Id}=$Exon;
}
1;
