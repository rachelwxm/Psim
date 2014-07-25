=head1 SequencingError

Copyright 2014, Rachel Wu. All Rights Resered.  
This program is free software. You may copy or redistribute 
it under the same terms as Perl itself.  

=head2 Version

=head4 Author: Wu Xiaomeng(Rachel) rachelwu123@gmail.com

=head4 Version: 1.5

=head4 Update: 2014-04-24

=head2 USAGE

Used for producing Sequencing Error,such as single base error, insert and delete.

user can set error type, error rate and the probabilities of different error types or use the default value(0.0005,0.34,0.33,0.33)

=head2 PARAMETERS

input:

	ReadSeq(string)
	ErrorRate		the probability of all sequencing error type
	    SingleBaseError	percentage of single base error
	    Insert			percentage of insert
	    Deletion		percentage of deletion
	reads name
	error array(array reference)

output:

	Error Read Seq

=head2 EXAMPLE

	&SequencingError($read,$Error,$prefix,\@ERROR);
	#$error=$ErrorRate:$SingleBaseError:$Insert:$Deletion

=cut

my %Error;
$Error{"A"}{"0"}="T";
$Error{"A"}{"1"}="G";
$Error{"A"}{"2"}="C";
$Error{"a"}{"0"}="T";
$Error{"a"}{"1"}="G";
$Error{"a"}{"2"}="C";
$Error{"T"}{"0"}="A";
$Error{"T"}{"1"}="G";
$Error{"T"}{"2"}="C";
$Error{"t"}{"0"}="A";
$Error{"t"}{"1"}="G";
$Error{"t"}{"2"}="C";
$Error{"G"}{"0"}="A";
$Error{"G"}{"1"}="T";
$Error{"G"}{"2"}="C";
$Error{"g"}{"0"}="A";
$Error{"g"}{"1"}="T";
$Error{"g"}{"2"}="C";
$Error{"C"}{"0"}="A";
$Error{"C"}{"1"}="T";
$Error{"C"}{"2"}="G";
$Error{"c"}{"0"}="A";
$Error{"c"}{"1"}="T";
$Error{"c"}{"2"}="G";

my %InsertSE=(
		"0" =>  "A",
		"1" =>  "T",
		"2" =>  "G",
		"3" =>  "C",
		);

package SequencingError;
require Exporter;
use warnings;
use strict;

our @ISA=qw(Exporter);
our @EXPORT=qw(SequencingError);
our @VERSION=1.5;

sub SequencingError
{
	my ($ReadSeq,$error,$prefix,$ErrorReport)=@_;
	my $ReadLength=length($ReadSeq);
	my @ErrorSites=&NumOfErrorSite($error,$ReadLength,rand(1));
	my $NumOfError=scalar(@ErrorSites);
	if($NumOfError==4 || $NumOfError==0)
	{
		return $ReadSeq;
	}
	else
	{
		foreach my $site(@ErrorSites)
		{
			my $ori=substr($ReadSeq,$site,1);
			my $new=$Error{$ori}{int(rand(3))};
			substr($ReadSeq,$site,1)=$new;
			my $info="$site\t$ori\t$new";
			push @$ErrorReport, "$prefix\t$info";
		}
		return $ReadSeq;
	}
}

sub NumOfErrorSite
{
	my ($ErrorRate,$ReadLength,$Probability)=@_;
	my $ZeroErrorSite=(1-$ErrorRate)**$ReadLength;
	my $OneErrorSite=$ReadLength*((1-$ErrorRate)**($ReadLength-1))*$ErrorRate;
	my $TwoErrorSite=$ReadLength*($ReadLength-1)/2*(1-$ErrorRate)**($ReadLength-2)*($ErrorRate**2);
	my $ThreeErrorSite=$ReadLength*($ReadLength-1)*($ReadLength-2)/3/2*((1-$ErrorRate)**($ReadLength-3))*($ErrorRate**3);
	if($Probability<$ZeroErrorSite)
	{
		return 0,0,0,0;
	}
	elsif(($Probability>=$ZeroErrorSite) && ($Probability<($ZeroErrorSite+$OneErrorSite)))
	{
		my $site=int(rand($ReadLength));   #error site
			return $site;
	}
	elsif($Probability>=($ZeroErrorSite+$OneErrorSite) && $Probability<($ZeroErrorSite+$OneErrorSite+$TwoErrorSite))
	{    #this read has two error bases
		my $site1=int(rand($ReadLength));
		my $site2=int(rand($ReadLength));
		while($site2==$site1){$site2=int(rand($ReadLength));}
		return $site1,$site2;
	}
	else
	{
		my $site1=int(rand($ReadLength));
		my $site2=int(rand($ReadLength));
		my $site3=int(rand($ReadLength));
		while($site2==$site1 || $site2==$site3 || $site1==$site3)
		{
			$site2=int(rand($ReadLength));
			$site3=int(rand($ReadLength))
		}
		return $site1,$site2,$site3;
	}
}
