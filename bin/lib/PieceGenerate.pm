=pod

=head1 PieceGenerate
Used to shear the reference genome randomly
parameter:
	whole genome length
	type of library
		0 for others
		$Insert\:length($Linker) for Roche 454 PE
	data coverage
	lamada of distribution
	sd of distribution
	length limit(optional)

example:
	&PieceGenerate($referencelength,$type,$coverage,$fragmean,$fragsd,$limit)

=cut

package PieceGenerate;
use Exporter;
use warnings;
use strict;
use normal;

our @ISA=qw(Exporter);
our @EXPORT=qw(PieceGenerate);
our $VERSION=1.0;

sub PieceGenerate
{
	#	print "TEST: CALLED PIECEGENERATE\n";
	my ($WL,$type,$circle,$Coverage,$Lamada,$sd,$limit)=@_;
	my @StartLeng;
	my $lengthlimit=$WL*$Coverage;
	my $num=int(($WL*$Coverage/$Lamada)+0.5);
	my @length=&normal($num,$Lamada,$sd,$limit);
	if($circle==0)
	{
		$WL=$WL-int($Lamada*1.5);
	}
	if($type==0) #other library types
	{
		foreach my $l(@length)
		{
			my $StartPoint=int(rand($WL));
			push @StartLeng,"$StartPoint\_$l";
		}
	}
	else #Roche 454 PE library
	{
		my ($InsertMean,$InsertSD,$LinkerLeng)=split /\:/,$type;
		my @insertleng=&normal($num,$InsertMean,$InsertSD);
		for my $i(0..$#length)
		{
			#length[$i]=fragment length
			#insertleng[$i]=insert size
			my $InsertSite=int(rand($WL));
			my $InsertEnd=$insertleng[$i]+$InsertSite-1;
			my $low=$insertleng[$i]+$InsertSite-$length[$i];
			my $StartPoint=int(rand($length[$i]+$LinkerLeng))+$low;
			if($StartPoint<($low+$LinkerLeng))
			{
				#reutrn 3, fragment start site, end site and whole fragment size. gap fill with linker sequence
				push @StartLeng,"3\_$StartPoint\_$InsertEnd\_$length[$i]";
			}
			elsif($StartPoint<=$InsertEnd)
			{
				#return 0,fragment start site, insert end site, insert start site and fragment size
				push @StartLeng,"0\_$StartPoint\_$InsertEnd\_$InsertSite\_$length[$i]";
			}
			else
			{
				#reutrn 5, insert start site, end site and whole fragment size. gap fill with linker sequence
				my $end=$InsertSite+$length[$i]+$StartPoint-$LinkerLeng-$InsertEnd;
				push @StartLeng,"5\_$InsertSite\_$end\_$length[$i]";
			}
		}
	}
	return @StartLeng;
}
1;
