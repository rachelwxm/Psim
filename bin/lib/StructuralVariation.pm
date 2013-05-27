=pod

=head1 StructuralVariation
This package is used for generating SV based on the reference sequence
Include:
	a.deletion          region(s)
	b.insertion         site(s),seq(s)
	c.inversion         region(s)
	d.tandem repeat     region(s),times
	e.translocation     ori_region,aim_begin_site
The sites mean to ori ref seq

SV can be randomly generate with average length and SV rate or set by
specific parameters like below:
	#deletion	begin_site2,length;begin_site2,leng..... different deletions separated by ";"
	deletion	1,25;400,88
	#insertion	insert_site1,seq;insert_site2,path2...
	insertion	50,sequence1;300,sequence2
	#inversion	begin_site1,length;begin_site2,length...
	invertion	100,20
	#tandem_repeat	begin_site1,length,times;begin_site2,length,times...
	tandem_repeat	150,100,4;320,30,2
	#translocation	begin_site1,length,insert site1;...
	translocation	200,50,310

in the main program:
    &StrcutVaretation(\$Sequence,$SVReport,$SeqName,$SVType,$SV,\@SV)
	$type--1 for input SV config file
		   0 for randomly generate with parameters
Example:
	$StrcutVaretation(\$refseq,$SVreport,"Chr1",0,"0.3:3000",\@SV)

=cut

package StructuralVariation;
require Exporter;
use normal;
use warnings;
use strict;

our @ISA=qw(Exporter);
our @EXPORT=qw(StructuralVariation);
our $VERSION=1.0;

sub StructuralVariation
{
#&StrcutVaretation(\$Sequence,$SeqName,$SVType,$SV,\@SV)
	my ($seq,$prefix,$SVtype,$SV,$report)=@_;
	my ($num_del,$num_inse,$num_inve,$num_rep,$num_tran,@del,@inse,@inve,@rep,@tran,%site);
#SV information from specific file
	if($SVtype==1)
	{
		open SVFILE,"<$SV" || die "cannot open SV file $SV!\n";
		while(<SVFILE>)
		{
			chomp;
			my @s=split /s\+/;
			if($s[1])
			{
				if($s[0]~~/del/i)
				{
					@del=split /\;/,$s[1];
					foreach(@del)
					{
						my @tmp1=split /,/;
						$site{$tmp1[0]}="del_".$tmp1[1];
					}
				}
				elsif($s[0]~~/inser/i)
				{
					@inse=split /\;/,$s[1];
					foreach(@inse)
					{
						my @tmp2=split /,/;
						$site{$tmp2[0]}="ins_".$tmp2[1];
					}
				}
				elsif($s[0]~~/inver/i)
				{
					@inve=split /\;/,$s[1];
					foreach(@inve)
					{
						my @tmp3=split /,/;
						$site{$tmp3[0]}="inv_".$tmp3[1];
					}
				}
				elsif($s[0]~~/repeat/i)
				{
					@rep=split /\;/,$s[1]; 
					foreach(@rep)
					{
						my @tmp4=split /,/;
						$site{$tmp4[0]}="rep_".$tmp4[1]."_".$tmp4[2];
					}
				}
				my @tran_seq;
				if($num_tran~~/tran/i)
				{
					@tran=split /\;/,$s[1];
					foreach my $i(0..$#tran)
					{
						my @tmp5=split /,/, $tran[$i];
						$tran_seq[$i]=substr($$seq,$tmp5[0],$tmp5[1]);    
						$site{$tmp5[0]}="del_".$tmp5[1];
						$site{$tmp5[2]}="ins_".$tran_seq[$i];
					}
				}
				@tran_seq=qw();
			}
		}
		close SVFILE;
	}
	else
	{
		my ($SVcoverage,$SVaveleng)=split /\:/,$SV;
		my $SVleng=0;
		my $WholeLength=length($$seq);
		my $LengthLimit=$WholeLength*$SVcoverage;
		for my $i(1..int($LengthLimit/$SVaveleng+10))
		{
#			my @length=&normal(1,$SVaveleng,1000,1000);
			my @length=&normal(1,$SVaveleng,1000);
			my $length=shift @length;
			my $start=int(rand($WholeLength));
			my $type=int(rand(4));
			if($type==0)  #deletion
			{
				push @$report, "$prefix\tdeletion\t$start\t$length";
				$site{$start}="del_".$length;

			}
			elsif($type==1)     #inversion
			{
				push @$report, "$prefix\tinversion\t$start\t$length";
				$site{$start}="inv_".$length;
			}
			elsif($type==2)         #tandem repeat
			{
				my $repe=int(rand(8))+1;
				push @$report, "$prefix\trepeat\t$start\t$length\t$repe";
				$site{$start}="rep_".$length."_".$repe;
			}
			else            #translocation
			{
				my $newsite=int(rand($WholeLength));
				push @$report, "$prefix\ttranslocation\t$start\t$length\t$newsite";
				$site{$start}="del_".$length;
				$site{$newsite}="ins_".substr($$seq,$start,$length);
			}
			$SVleng+=$length;
			last if($SVleng>$LengthLimit)
		}
	}

#handle SV
	my $site_pre;
	my %contig;
	my @site_sort=sort {$a<=>$b} keys %site;
	foreach (@site_sort)
	{
		if($site_pre)
		{
			$contig{$site_pre}=substr ($$seq,$site_pre,($_-$site_pre));
		}
		$site_pre=$_;
	}
	$contig{$site_sort[-1]}=substr($$seq,$site_sort[-1],(length($$seq)-$site_sort[-1]));

	while (my ($key,$value)=each %site)
	{
		if($value=~/del/)
		{
			substr($contig{$key},0,(split/_/,$value)[1])="";
		}
		elsif($value=~/ins/)
		{
			$contig{$key}=(split/_/,$value)[1].$contig{$key};
		}
		elsif($value=~/inv/)
		{
			my $tmp=substr($contig{$key},0,(split/_/,$value)[1]);
			substr($contig{$key},0,(split/_/,$value)[1])=reverse $tmp;
		}
		elsif($value=~/rep/)
		{ 
			my $tmp=substr($contig{$key},0,(split/_/,$value)[1]);
			my $rep_str=$tmp x ((split /_/,$value)[2]-1);
			$contig{$key}=$rep_str.$contig{$key};
		}
	}
	my $return_seq;
	if($site_sort[0]==0)
	{
		foreach(@site_sort)
		{
			$return_seq.=$contig{$_};
		}
	}
	else
	{
		$return_seq=substr($$seq,0,$site_sort[0]-1);
		foreach(@site_sort)
		{
			$return_seq.=$contig{$_};
		}
	}
	return $return_seq;
}

