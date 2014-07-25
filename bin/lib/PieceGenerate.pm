=head1 PieceGenerate

Copyright 2014, Rachel Wu. All Rights Resered.  
This program is free software. You may copy or redistribute 
it under the same terms as Perl itself.  

=head2 Version

=head4 Author: Wu Xiaomeng(Rachel) rachelwu123@gmail.com

=head4 Version: 1.5

=head4 Update: 2014-04-20

=head2 USAGE

Used to shear the reference genome into pieces.
In Version 1.4, damage DNA length consist of two normal distribution.
And add parameter a, k, mu et al.

=head2 PARAMETERS

	WL->seqeucne length
	type->$ltype
	    0 for others
	    $Insert\:length($Linker) for Roche 454 PE
	circle->reference sequence circle or not
	Coverage->average sequencing coverage
	readmean->average read legnth
	fragmean->average shared fragment length
	fragsd->standard deviation of shared fragment length
	fraglim->length limit of fragments (optional)
	duplp->Duplication Percent
	dupll->Duplication Times Limit(optional)
	dupla->Power Law parameter a
	duplb->Power Law parameter b
	dpa->Damage parameter A
	dpk->Damage parameter K

=head2 EXAMPLE

	#damage DNA library perparation
	@StartLeng=&PieceGenerate("WL",length($$Sequence),"circle",$Circle,"Coverage",$Coverage,"fragmean",$FragmentMean,"fragsd",$FragmentSD,"fraglim",$FragLim,"dpa",$DamageParaA,"dpk",$DamageParaK,"dpmu",$DamageParaMu);

	@StartLeng=&PieceGenerate("WL",length($$Sequence),"type",$ltype,"circle",$Circle,"Coverage",$Coverage,"readmean",$readmean,"fragmean",$FragmentMean,"fragsd",$FragmentSD,"fraglim",$FragLim,"duplp",$duplpercent,"dupll",$dupllimit,"dupla",$a,"duplb",$b);

	@StartLeng=&PieceGenerate("WL",length($$Sequence),"type",$ltype,"circle",$Circle,"Coverage",$Coverage,"readmean",$readmean,"fragmean",$FragmentMean,"fragsd",$FragmentSD,"fraglim",$FragLim,"duplp",$duplpercent,"dupll",$dupllimit);

=cut

package PieceGenerate;
use Exporter;
use warnings;
use strict;
use normal;
use PowerLaw;

our @ISA=qw(Exporter);
our @EXPORT=qw(PieceGenerate);
our $VERSION=1.5;

sub PieceGenerate
{
	my %hash=@_;
#%hash keys: 
#	WL->length($sequence)
#	type->$ltype
#	circle->$circle
#	Coverage->$Coverage
#	readmean->$readmean
#	fragmean->$fragmentmean
#	fragsd->fragmentsd
#	fraglim->fragmentlimit
#	duplp->DuplicationPercent
#	dupll->DuplicationTimesLimit
#	dupla->PowerLawparametera
#	duplb->PowerLawparameterb
#	dpa->DamageParaA
#	dpk->DamageParaK
#	dpmu->DamageParaMu

	my @StartLeng;
	my @length;
	my $lengthlimit=$hash{WL}*$hash{Coverage};
	my $num;
	if($hash{dpa})
	{
		$hash{type}=0;
		$num=int(($hash{WL}*$hash{Coverage}/$hash{fragmean})+0.5);
		my $half=int($num/2+0.4);
		my $last=$num-$half;
		my $num1=$half;
		my $num2=int($last/$hash{dpk}+0.4);
		my $mean1=$hash{fragmean};
		my $sd1=$hash{dpa}*$hash{fragsd};
		my $mean2=$hash{fragmean}-$hash{fragsd}*sqrt(2*log($hash{dpa}*$hash{dpk}))+$hash{dpmu};
		my $sd2=$hash{fragsd};

		my $limhalf="$hash{fraglim}$hash{fragmean}\-";
		my $type=1;
		my @array1=&normal("size",$num1,"mean",$mean1,"sd",$sd1,"limit",$limhalf,"type",$type);
		my @array2=&normal("size",$num2,"mean",$mean2,"sd",$sd2,"limit","$hash{fragmean}\+","type",$type);

		push @length,@array1;
		push @length,@array2;
		for(1..($last-$num2))
		{
			my $ran=int(rand($num2));
			push @length,$array2[$ran];
		}
	}
	else
	{
		if($hash{duplp})
		{
			$num=int((1-$hash{duplp})*$hash{WL}*$hash{Coverage}/$hash{readmean}+0.5);
		}
		else
		{
			$num=int(($hash{WL}*$hash{Coverage}/$hash{readmean})+0.5);
		}
		@length=&normal("size",$num,"mean",$hash{fragmean},"sd",$hash{fragsd},"limit",$hash{fraglim});
	}

	if($hash{type} eq 0) #other library types
	{
		$hash{WL}=$hash{WL}-int($hash{fragmean}*1.5) if($hash{circle} eq 0);
		foreach my $l(@length)
		{
			my $StartPoint=int(rand($hash{WL}));
			push @StartLeng,"$StartPoint\_$l";
		}
	}
	else #Roche 454 PE library
	{
		my ($InsertMean,$InsertSD,$LinkerLeng)=split /\:/,$hash{type};
		my @insertleng=&normal("size",$num,"mean",$InsertMean,"sd",$InsertSD);
		$hash{WL}=$hash{WL}-int($InsertMean*1.5) if($hash{circle} eq 0);
		for my $i(0..$#length)
		{
			my $InsertSite=int(rand($hash{WL}));
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
				my $end=$InsertSite+$length[$i]+$StartPoint-$LinkerLeng-$InsertEnd-1;
				push @StartLeng,"5\_$InsertSite\_$end\_$length[$i]";
			}
		}
	}
	if($hash{duplp})
	{
		my $duplfrag=int($hash{duplp}/(1-$hash{duplp})*$num+0.5);
		$hash{dupla}||=0.35;
		$hash{duplb}||=-2;
		my @dupltimes=&PowerLaw("r".$duplfrag,$hash{dupla},$hash{duplb},$hash{dupll});
		foreach my $dt(@dupltimes)
		{
			my $j=int(rand($num));
			for(1..$dt)
			{
				push @StartLeng,$StartLeng[$j];
			}
		}
	}
	return @StartLeng;
}
1;
