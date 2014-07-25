=head1 SNP

Copyright 2014, Rachel Wu. All Rights Resered.  
This program is free software. You may copy or redistribute 
it under the same terms as Perl itself.  

=head2 Version

=head4 Author: Wu Xiaomeng(Rachel) rachelwu123@gmail.com

=head4 Version: 1.5

=head4 Update: 2014-04-17

=head2 USAGE

Used for generate SNP ref seq with given sites, bases
       and rate. User can offer the SNP sites, bases and 
       rate, or random generate the SNP bases in default rate.
       edit reference directly on ref seq.

=head2 PARAMETERS

	seq-> reference sequence(string reference)
	seqname->sequence name
	snptype->SNP type
	    1 for user set
	    0 for random setting	
	snp->SNP site & base file or SNP Rate	
	    SNP file format as: SNP Sites  SNP Base  rate
	      such as: 1462   A,T 0.3,0.2
	               18209  C   0.4
	    SNP Rate: percentage of SNP bases
	snpinfo->SNP report array(array reference)

=head2 EXAMPLE

	#$snptype=1 and $SNP contain user set snp info
	&SNP("seq",\$SEQUENCE,"snptype",$snptype,"seqname",$SeqName,"snp",$SNP,"snpinfo",\@SNP);
	#$snptype=0 and $SNP=SNP rate(0~1)
	&SNP("seq",\$SEQUENCE,"snptype",$snptype,"seqname",$SeqName,"snp",$SNP,"snpinfo",\@SNP);

=cut

package SNP;
require Exporter;
use warnings;
use strict;

our @ISA=qw(Exporter);
our @EXPORT=qw(SNP);
our $VERSION=1.5;

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

sub SNP
{
	my %hash=@_;
##====PARAMETER EXPLANATION====
#seq-> reference sequence(string reference)
#snptype->snptype
#seqname->sequence name
#snp->SNP info
#snpinfo->SNP report array(array reference)
	if($hash{snptype} eq "1")
	{
		open(SNPFILE,"<$hash{snp}") || die "$!\nCannot open SNP information file $hash{snp}\n";
		while(<SNPFILE>)
		{
			#print "EN $_";
			#print "round\t";
			chomp;
			my ($Reference,$SNPsite,$SNPbase,$SNPrate)=split /\s+/;
			if($Reference eq $hash{seqname})
			{
				my $r=rand(1);##print $r,"\t";
				my $Ini=substr(${$hash{seq}},($SNPsite-1),1);
				if ($SNPbase=~/,/)
				{
					#print "multi\t";
					my @SNPbase=split /,/,$SNPbase;
					my @SNPrate=split /,/,$SNPrate;
					my $SNPnum=scalar(@SNPbase);
					die "Wrong format of SNP file\n" if($SNPnum!=scalar(@SNPrate));
					if($r<=$SNPrate[0])
					{
						#print "yes one\t";
						substr(${$hash{seq}},($SNPsite-1),1)=$SNPbase[0];
						push @{$hash{snpinfo}},"$hash{seqname}\t$SNPsite\t$Ini\t$SNPbase[0]";
					}
					elsif($r>$SNPrate[0] || $r<=($SNPrate[0]+$SNPrate[1]))
					{
						#print "yes two\t";
						substr(${$hash{seq}},($SNPsite-1),1)=$SNPbase[1];
						push @{$hash{snpinfo}},"$hash{seqname}\t$SNPsite\t$Ini\t$SNPbase[1]";
					}
					elsif(exists $SNPrate[2])
					{
						if($r>($SNPrate[0]+$SNPrate[1]) && $r<=($SNPrate[0]+$SNPrate[1]+$SNPrate[2]))
						{
						#print "yes three\t";
							substr(${$hash{seq}},($SNPsite-1),1)=$SNPbase[2];
							push @{$hash{snpinfo}},"$hash{seqname}\t$SNPsite\t$Ini\t$SNPbase[2]";
						}
					}
				}
				else
				{
					#print "single\t";
					if($r<=$SNPrate)
					{
						#print "yes\t";
						substr(${$hash{seq}},($SNPsite-1),1)=$SNPbase;
						push @{$hash{snpinfo}}, "$hash{seqname}\t$SNPsite\t$Ini\t$SNPbase";
					}
				}
				#print "\n";
			}
		}
		close SNPFILE;
	}
	else
	{
		my $RefLength=length(${$hash{seq}});
		my $SNPCount=int($hash{snp}*$RefLength+0.5);
		my %SNPSite;
		for(1..$SNPCount)
		{
			my $Site=int(rand($RefLength));
#		print "test reflength $RefLength\n";
			my $Ini=substr(${$hash{seq}},$Site,1);
			my $a=int(rand(3));
			while(!exists $Error{$Ini}{$a})
			{
				$Site=int(rand($RefLength));
#				print "test site is $Site\n";
				$Ini=substr(${$hash{seq}},$Site,1);
			}
			substr(${$hash{seq}},$Site,1)=$Error{$Ini}{$a};
			push @{$hash{snpinfo}},"$hash{seqname}\t$Site\t$Ini\t$Error{$Ini}{$a}";
		}
	}
}

1;
