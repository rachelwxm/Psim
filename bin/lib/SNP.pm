=pod

=head1 SNP
Used for producing SNP ref seq with given sites,bases and rate
User can offer the SNP sites, bases and rate, or random generate
the SNP base in default rate
edit reference directly on ref seq

Paraments:
    Ref Seq
    SNP Type
        1 for user set
        0 for random setting
    SNP information file
        SNP Sites   SNP Base    rate
        such as: 1462   A,T 0.3,0.2
                 18209  C   0.4
    SNP Rate	percentage of SNP bases
    SNP report array
Example:
    &SNP(\$ref,$SNPreport,$prefix,1,file,\@SNP)
    &SNP(\$ref,$SNPreport,$prefix,0,$SNPrate,\@SNP)
=cut

package SNP;
require Exporter;
use warnings;
use strict;

our @ISA=qw(Exporter);
our @EXPORT=qw(SNP);
our $VERSION=1.0;

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
#my $n=0;

sub SNP
{
#	$n+=1;
	my ($RefSeq,$SNPType,$prefix,$SNP,$report)=@_;
#	print "TEST times: $n\n";
#	map {print "$_\n"} @_;
	if($SNPType eq "1")
	{
		open(SNPFILE,"<$SNP") || die "$!\nCannot open SNP information file $SNP\n";
		while(<SNPFILE>)
		{
			chomp;
			my ($SNPsite,$SNPbase,$SNPrate)=split /\s+/;
			my $r=rand(1);
			my $Ini=substr($$RefSeq,($SNPsite+1),1);
			if ($SNPbase=~/,/)
			{
				my @SNPbase=split /,/,$SNPbase;
				my @SNPrate=split /,/,$SNPrate;
				my $SNPnum=scalar(@SNPbase);
				die "Wrong format of SNP file\n" if($SNPnum!=scalar(@SNPrate));
				if($r<=$SNPrate[0])
				{
					substr($$RefSeq,($SNPsite+1),1)=$SNPbase[0];
					push @$report,"$prefix\t$SNPsite\t$Ini\t$SNPbase[0]";
				}
				elsif($r>$SNPrate[0] || $r<=($SNPrate[0]+$SNPrate[1]))
				{
					substr($$RefSeq,($SNPsite+1),1)=$SNPbase[1];
					push @$report,"$prefix\t$SNPsite\t$Ini\t$SNPbase[1]";
				}
				if(exists $SNPrate[2])
				{
					if($r>($SNPrate[0]+$SNPrate[1]) && $r<=($SNPrate[0]+$SNPrate[1]+$SNPrate[2]))
					{
						substr($$RefSeq,($SNPsite+1),1)=$SNPbase[2];
						push @$report,"$prefix\t$SNPsite\t$Ini\t$SNPbase[2]";
					}
				}
			}
			else
			{
				if($r<=$SNPrate)
				{
					substr($$RefSeq,($SNPsite+1),1)=$SNPbase;
					push @$report, "$prefix\t$SNPsite\t$Ini\t$SNPbase";
				}
			}
		}
		close SNPFILE;
	}
	else
	{
		my $RefLength=length($$RefSeq);
		my $SNPCount=int($SNP*$RefLength+0.5);
		my %SNPSite;
		for(1..$SNPCount)
		{
			my $Site=int(rand($RefLength));
			my $Ini=substr($$RefSeq,$Site,1);
			my $a=int(rand(3));
			while(!exists $Error{$Ini}{$a})
			{
				$Site=int(rand($RefLength));
				$Ini=substr($RefSeq,$Site,1);
			}
			substr($$RefSeq,$Site,1)=$Error{$Ini}{$a};
			push @$report,"$prefix\t$Site\t$Ini\t$Error{$Ini}{$a}";
		}
	}
}

1;
