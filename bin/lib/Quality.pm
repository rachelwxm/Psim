=pod 

=head1 Quality
Used for producing quality score for reads.

FASTQ quality format:
 SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
 ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
 ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
 .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
 !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
 |                         |    |        |                              |                     |
33                        59   64       73                            104                   126
 0........................26...31.......40                                
                          -5....0........9.............................40 
                                0........9.............................40 
                                   3.....9.............................40 
 0........................26...31........41                               

  S - Sanger        Phred+33,  raw reads typically (0, 40)
  X - Solexa        Solexa+64, raw reads typically (-5, 40)
  I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
  J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
	  with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
	  (Note: See discussion above).
  L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

*above fastq information from wikipedia http://en.wikipedia.org/wiki/FASTQ_format
parameter:
	sequencing technology(Illumina, Roche or SOLiD)
	\$quality
	length

example:
	&Quality($Command,\$quality,$readleng,\@qual_th);

=cut

package Quality;
require Exporter;
use warnings;
use strict;

our @ISA=qw(Exporter);
our @EXPORT=qw(Quality);
our @VERSION=1.0;

my @conv_table;
for (-64..64)
{
	$conv_table[$_+64] = chr(int(33 + 10*log(1+10**($_/10.0))/log(10)+.499));
}

sub Quality
{
	my ($Type,$quality,$length,$qual_th)=@_;
	#convert Dec code to ASCII code for random quality value (solexa or std/standard quality)
	if($Type =~ /solid/i)
	{
		my $q1=int(rand(40-$$qual_th[1])+$$qual_th[1]+1);
		my $q2=int(rand(40-$$qual_th[2])+$$qual_th[2]+1);
		if ($length>$$qual_th[0])
		{
			for (my $i=0;$i<$length-$$qual_th[0];$i++)
			{
				my $q_tmp = int(rand($$qual_th[3])+1+$q1);
				$q_tmp = 40 if ($q_tmp>40);
				$$quality .= sprintf("%02d ",$q_tmp);
			}
			for (my $i=$length-$$qual_th[0];$i<$length;$i++)
			{
				my $q_tmp = int(rand($$qual_th[3])+1+$q2);
				$q_tmp = 40 if ($q_tmp>40);
				$$quality .= sprintf("%02d ",$q_tmp);
			}
		}
	}
	else
	{
		if ($length>$$qual_th[0])
		{
			my $q1=int(rand(40-$$qual_th[1])+$$qual_th[1]+1);
			my $q2=int(rand(40-$$qual_th[2])+$$qual_th[2]+1);
			for (my $i=0;$i<$length-$$qual_th[0];$i++)
			{
				if ($Type =~ /roche/i)
				{
					my $q_tmp = int(rand($$qual_th[3])+$q1);
					$q_tmp = 40 if ($q_tmp>40);
					$$quality .= $conv_table[$q_tmp+64]; #std quality
				}
				elsif ($Type =~ /illumina/i)
				{
					my $q_tmp = int(rand($$qual_th[3])+$q1);
					$q_tmp = 40 if ($q_tmp>40);
					$$quality .= chr($q_tmp+64); #solexa quality
				}
			}
			for (my $i=$length-$$qual_th[0];$i<$length;$i++)
			{
				if ($Type =~ /roche/i)
				{
					my $q_tmp = int(rand($$qual_th[3])+$q2);
					$q_tmp = 40 if ($q_tmp>40);
					$$quality .= $conv_table[$q_tmp+64]; #std quality
				}
				elsif ($Type =~ /illumina/i)
				{
					my $q_tmp = int(rand($$qual_th[3])+$q2);
					$q_tmp = 40 if ($q_tmp>40);
					$$quality .= chr($q_tmp+64); #solexa quality
				}
			}
		}
	}
}
1;
