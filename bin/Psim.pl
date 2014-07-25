#!/usr/bin/perl -w
use strict;
use Math::Random qw(random_poisson);
use Getopt::Long;
use Term::ANSIColor;
use FindBin qw($Bin $Script);
use lib "$FindBin::Bin/lib";
use 5.010;
use Data::Dumper;
use File::Basename qw(basename dirname);
use Overhanging;
use normal;
use PieceGenerate;
use ReverseComplement;
use SNP;
use StructuralVariation;
use Base2Col;
use Quality;
use SequencingError;
use ReadGff;
use PCRduplication;

#==========================USER PARAMETERS==============================
my ($Help,$Command,$Reference,$Circle,$SNP,$SV,$Coverage,$PE);
my ($FragmentMean,$FragmentSD,$FragLim);
my ($ReadLeng,$ReadSD);
my ($LostSingle,$DoublePercent,$Mismatch,$Damage,$MutationSite,$MutationArray,$AmpMean,$AmpSD,%MutationRate,%MutationArray,$Efficiency);
my ($Adapter,$Output,$Dir,$Library);
my ($Insert,$Linker,$Error,$adapter1,$adapter2);
my ($Qmean,$Qsd,$Qtype,$Header,$Example);
my ($GFF,$CDNA,$CovcDNA);
my $count=0;
my ($ploidy,$print);
my ($duplpercent,$dupllimit,$duplpara);
my $noamp;
$Command=$ARGV[0];
GetOptions(
	"noamp"		=>\$noamp,
	"h"			=>\$Help,
	"help"		=>\$Help,
	"ref:s"		=>\$Reference,
	"circle"	=>\$Circle,
	"snp:s"		=>\$SNP,
	"sv:s"		=>\$SV,
	"cov:s"		=>\$Coverage,
	"pe"		=>\$PE,
	"insert:s"	=>\$Insert,
	"linker:s"	=>\$Linker,
	"fragmean:s"=>\$FragmentMean,
	"fragsd:s"	=>\$FragmentSD,
	"fraglim:s"	=>\$FragLim,
	"read:s"	=>\$ReadLeng,
	"readsd:s"	=>\$ReadSD,
	"adapter:s"	=>\$Adapter,
	"damage"	=>\$Damage,
	"overhang:s"=>\$Mismatch,
	"library:s"	=>\$Library,
	"mutsite:s"	=>\$MutationSite,
	"mutarray:s"=>\$MutationArray,
	"ds:s"		=>\$DoublePercent,
	"lost:s"	=>\$LostSingle,
	"dir:s"		=>\$Dir,
	"output:s"	=>\$Output,
	"error:s"	=>\$Error,
	"qmean:s"	=>\$Qmean,
	"qsd:s"		=>\$Qsd,
	"header:s"	=>\$Header,
	"example"	=>\$Example,
	"gff:s"		=>\$GFF,
	"cdna"		=>\$CDNA,
	"covcdna:s"	=>\$CovcDNA,
	"qtype:s"	=>\$Qtype,
	"ampmean:s"	=>\$AmpMean,
	"ampsd:s"	=>\$AmpSD,
	"effic:s"	=>\$Efficiency,
	"plo:s"		=>\$ploidy,
	"pri"		=>\$print,
	"duplp:s"	=>\$duplpercent,
	"dupll:s"	=>\$dupllimit,
	"duplpara:s"=>\$duplpara,
	);

my $USAGE="
        =================================================================
      || WARNINGS:                                                       ||
      ||     YOU SHOLD HAVE ENOUGH FREE SPACE AND MEMORY!                ||
      ||     if input reference sequence size is X and Y-fold coverages. ||
      || You shold have at least X*Y free space and X memory.            ||
        =================================================================
USAGE:   Used to generate NGS data. The sequencing platform include Illumina
         Roche 454 and SOLiD. 
         Input reference sequence could be one seq or multiple seqs. And the 
         sequencing sample can have overhang end like palogenome.

Author:  RachelWu\(rachelwu123\@gmail.com\)\, BENM\(BinxiaoFeng\@gmail.com\)
Version: v1.6
Update:  2014-04-17

Run:     perl $0 [command] [options]

command: variation  generate variation (SNP, InDel and SV) sequence(s) aginest reference
         illumina   generate Illumina sequencing data
         roche      generate Roche 454 FLX sequencing data
         solid      generate SOLiD sequencing data
         example    show running examples
****************************************************************************
v1.5\t2014-04-17\tRachel Wu
v1.6\t2014-07-25\tRachel Wu\tfix bug in quality.pm
\n";

my $variation="
*******************************************************************************************************
USAGE:  generate polyploid reference(s) with variation (SNP, InDel and SV) aginest given ref seq
Run:    perl $0 variation [options]
OPTIONS:
    --ref <FILE>          reference sequence file
    --plo <NUM>           ploidy of reference sample(integer)
    --circle              the reference sequence is in circle
    --snp <NUM|FILE>      random snp rate(default=0.001) or specific snp site and rate
                          file(format refer to ../example/snp_example.txt)
    --sv <NUM|FILE>       structure variation rate and average length(default=0.1:3000)
                          or specific sv type and site file(format refer to ../example/sv_example.txt)

  **output parameters**
    --dir <CHAR>          output directory(default=\".\/\")
    --output <CHAR>       prefix of output file(default=\"VariationSequenceByPsim.fa\")
*******************************************************************************************************
";

my $illumina="
*******************************************************************************************************
USAGE:  Simulate Illumina sequencing data
Run:    perl $0 illumina [options]
OPTIONS:
    --ref <FILE>          reference sequence file
    --plo <NUM>           ploidy of reference sample(integer)
    --pri                 print each variant reference sequence(default: no)
    --circle              the reference sequence is in circle
    --snp <NUM|FILE>      random snp rate(default=0.001) or specific snp site and rate
                          file(format refer to ../example/snp_example.txt)
    --sv <NUM|FILE>       structure variation rate and average length(default=0.1:3000)
                          or specific sv type and site file(format refer to ../example/sv_example.txt)
    --cov <NUM>           sequencing coverage of reference sequence(s)(default=3)
    --pe                  PE sequencing or not
    --fragmean <INT>      average length of library fragment(default=200)
    --fragsd <INT>        standard deviation of library fragment(default=10)
    --fraglim <INT>       limit length of fragment library(\"20+\" means must above 20nt, and
                          \"240-\" means must shorter than 240nt,if(-damage) this para default=20+)
    --read <INT>          average length of reads(default=100)
    --adapter <CHAR>      adapter sequence(default seq refer to ../example/adapterIllumina_example.txt, 
                          split two adapters by \"\:\")
    --error <NUM>         sequencing error rate of single base error (default=\"0.0005\")
    --qtype <CHAR>        quality type(\!=Sanger format, \@=illumina 1.3~1.8- format, default=\!)
    --qmean <NUM>         peak value of quality score(default=37)
    --qsd <NUM>           standard deviation of quality score(default=1)
    --duplp <NUM>         duplication reads percent
    --dupll <NUM>         duplication times limit
    --duplpara <NUM,NUM>  duplication power low parameters, format as a-b\(such as \"0.3,-2\"\)
                          warnings: a should be less than 1 and b should be a negative number

  **damage sequence parameters**
    --damage              the sequencing sample have overhanging and injured end
    --overhang <INT>      overhang length range of double strain(default=\"3\_20\")
    --library <INT>       library preparation type(1=single strain,2=double strain,default=1)
    --ds <NUM>            rate of lost one strain of double strain library(default=0.1)
    --lost <NUM>          lost rate single strain while library=1(default=0.5)
    --ampmean <NUM>       average amplification times(default=850)
    --ampsd <NUM>         standard deviation of amplification times
    --mutarray <FILE>     mutation rate array file(format refer to ../example/mutarray_example, 
                          default  100% C->T)
    --mutsite <FILE>      mutation possibility along fragment site(format refer to ../example/mutsite_example,
                          default 0.01)
    --effic <NUM>         efficiency of fill-in reaction(default=0.5)

  **RNA-SEQ sequencing parameters**
    --cdna                RNA-SEQ sequencing
    --gff <FILE>          input gene gff file
    --covcdna <INT>       mean coverage of cDNA sequences(default=10 if --cdna)

  **output parameters**
    --dir <CHAR>          output directory(default=\".\/\")
    --output <CHAR>       prefix of output file(default=\"Psim_Illumina\")
*******************************************************************************************************
";

my $roche="
*******************************************************************************************************
USAGE:  Roche 454 FLX sequencing data
Run:    perl $0 roche [options]
OPTIONS:
    --ref <FILE>          reference sequence file
    --plo <NUM>           ploidy of reference sample(integer)
    --pri                 print each variant reference sequence(default: no)
    --circle              the reference sequence is in circle
    --snp <NUM|FILE>      random snp rate(default=0.001) or specific snp site and rate
                          file(format refer to ../example/snp_example.txt)
    --sv <NUM|FILE>       structure variation rate and average length(default=0.1:3000)
                          or specific sv type and site file(format refer to ../example/sv_example.txt)
    --cov <NUM>           sequencing coverage of reference sequence(s)(default=3)
    --pe                  PE sequencing or not
    --insert <INT>        average insert size and sd(default=8000:30 if --pe)
    --linker <CHAR>       PE insert seq(default=\"ATAACTTCGTATAATGTATGCTATACGAAGTTAT\")
    --fragmean <INT>      average length of library fragment(equal to reads length,default=450)
    --fragsd <INT>        standard deviation of library fragment(default=50)
    --fraglim <INT>       limit length of fragment library(\"20+\" means must above 20nt, and
                          \"800-\" means must shorter than 800nt,if(-damage) this para default=20+)
    --error <NUM>         sequencing error rate of single base error (default=\"0.0005\")
    --qmean <NUM>         peak value of quality score(default=37)
    --qsd <NUM>           standard deviation of quality score(default=1)
    --duplp <NUM>         duplication reads percent
    --dupll <NUM>         duplication times limit
    --duplpara <NUM,NUM>  duplication power low parameters, format as a-b\(such as \"0.3,-2\"\)
                          warnings: a should be less than 1 and b should be a negative number

  **damage sequence parameters**
    --damage              the sequencing sample have overhanging and injured end
    --overhang <INT>      overhang length range of double strain(default=\"3\_20\")
    --library <INT>       library preparation type(1=single strain,2=double strain,default=1)
    --ds <NUM>            rate of lost one strain of double strain library(default=0.1)
    --lost <NUM>          lost rate single strain while library=1(default=0.5)
    --ampmean <NUM>       average amplification times(default=850)
    --ampsd <NUM>         standard deviation of amplification times
    --mutarray <FILE>     mutation rate array file(format refer to ../example/mutarray_example, 
                          default  100% C->T)
    --mutsite <FILE>      mutation possibility along fragment site(format refer to ../example/mutsite_example,
                          default 0.01)
    --effic <NUM>         efficiency of fill-in reaction(default=0.5)

  **RNA-SEQ sequencing parameters**
    --cdna                RNA-SEQ sequencing
    --gff <FILE>          input gene gff file
    --covcdna <INT>       mean coverage of cDNA sequences(default=10 if --cdna)

  **output parameters**
    --dir <CHAR>          output directory(default=\".\/\")
    --output <CHAR>       prefix of output file(default=\"Psim_Roche\")
*******************************************************************************************************
";

my $solid="
*******************************************************************************************************
USAGE:  SOLiD sequencing data
Run:    perl $0 solid [options]
OPTIONS:
    --ref <FILE>          reference sequence file
    --plo <NUM>           ploidy of reference sample(integer)
    --pri                 print each variant reference sequence(default: no)
    --circle              the reference sequence is in circle
    --snp <NUM|FILE>      random snp rate(default=0.001) or specific snp site and rate
                          file(format refer to ../example/snp_example.txt)
    --sv <NUM|FILE>       structure variation rate and average length(default=0.1:3000)
                          or specific sv type and site file(format refer to ../example/sv_example.txt)
    --cov <NUM>           sequencing coverage of reference sequence(s)(default=3)
    --pe                  PE sequencing or not
    --fragmean <INT>      average length of library fragment(default=100)
    --fragsd <INT>        standard deviation of library fragment(default=10)
    --fraglim <INT>       limit length of fragment library(\"20+\" means must above 20nt, and
                          \"240-\" means must shorter than 240nt,if(-damage) this para default=20+)
    --read <INT>          average length of reads(default=50)
    --error <NUM>         sequencing error rate of single base error (default=\"0.0005\")
    --qmean <NUM>         peak value of quality score(default=37)
    --qsd <NUM>           standard deviation of quality score(default=1)
    --header <CHAR>       sequencing header base(default=G)
    --duplp <NUM>         duplication reads percent
    --dupll <NUM>         duplication times limit
    --duplpara <NUM,NUM>  duplication power low parameters, format as a-b\(such as \"0.3,-2\"\)
                          warnings: a should be less than 1 and b should be a negative number

  **damage sequence parameters**
    --damage              the sequencing sample have overhanging and injured end
    --overhang <INT>      overhang length range of double strain(default=\"3\_20\")
    --library <INT>       library preparation type(1=single strain,2=double strain,default=1)
    --ds <NUM>            rate of lost one strain of double strain library(default=0.1)
    --lost <NUM>          lost rate single strain while library=1(default=0.5)
    --ampmean <NUM>       average amplification times(default=850)
    --ampsd <NUM>         standard deviation of amplification times
    --mutarray <FILE>     mutation rate array file(format refer to ../example/mutarray_example, 
                          default  100% C->T)
    --mutsite <FILE>      mutation possibility along fragment site(format refer to ../example/mutsite_example,
                          default 0.01)
    --effic <NUM>         efficiency of fill-in reaction(default=0.5)

  **RNA-SEQ sequencing parameters**
    --cdna                RNA-SEQ sequencing
    --gff <FILE>          input gene gff file
    --covcdna <INT>       mean coverage of cDNA sequences(default=10 if --cdna)

  **output parameters**
    --dir <CHAR>          output directory(default=\".\/\")
    --output <CHAR>       prefix of output file(default=\"Psim_SOLiD\")
*******************************************************************************************************
";

my $ExampleInfo="======================================EXAMPLE================================================
#randomly generate reference sequence(s) 
perl GenerateSeq.pl 

#generate SNP infomation file(snp-config)
perl lib/SNPInfoGenerate.pl  -f ../example/example.fa -p 0.05 -r 0.7,0.3 > ../example/snp-config

#generate sequencing data accroding to poluploid(assume 4) sequences based on one reference sequence(../example/example.fa), print each variant sequences to a new file
perl Psim.pl illumina --ref ../example/example.fa  --plo 4 --cov 10 --snp ../example/snp-config --sv 0 --pri --error 0.01

#generate roche SE sequencing data with no SV or SNP.
#save the output files to ../output/
perl Psim.pl roche --ref ../example/example.fa --snp 0 --sv 0 --dir ../output/

#generate illumina PE sequencing data.
#reference sequence is circle.
#sv information in the file of ../example/sv_example.txt
#snp information in the file of ../example/snp_example.txt
#coverage of reference is 10X
#Sequencing with sequencing errors (~0.5%) for each reads via stochastic probability.
#save the output files to ../output/
perl Psim.pl illumina --ref ../example/example.fa --snp ../example/snp_example.txt --sv ../example/sv_example.txt 
 --pe --circle 1 --cov 10 --error 0.005 --dir ../output/

#sv information
#one reference sequence named with \"ChromosomeOne\". It has a 2600bp deletion start from site 300, a 3300bp deletion start from site 900000, a 3000bp inversion start from site 37000, a 1200bp repeat start from stie 26000 and repeat time of 4, a 5000bp translocation with original start site of 43000 and new start site of 60000 and insert \"TTTTTTGGGGGGGGGCCCCA\" to the site of 16350.
#in the sv config file:
 ChromosomeOne	deletion	300,2600;900000,3300
 ChromosomeOne	inversion	37000,2600
 ChromosomeOne	tandem_repeat	26000,1200,4
 ChromosomeOne	translocation	43000,5000,60000	
 ChromosomeOne	insertion	16350,TTTTTTGGGGGGGGGCCCCA

#snp information
#one reference sequence named with \"ChromosomeOne\".
#The base at site 1462 has 30% possibility of mutation into \"A\" and 20% possibility of mutation into \"T\"
#The base at site 18209 has 40% possibility of mutation into \"C\"
#The base at site 2840 has 20% possibility of mutation into \"A\", 20% possibility of mutation into \"T\" and 30% possibility of mutation into \"G\"
#in the snp config file:
 ChromosomeOne	1462	A,T	0.3,0.2
 ChromosomeOne	18209	C	0.4
 ChromosomeOne	2840	A,T,G	0.2,0.2,0.3

#damage1
#reference sequence is ../example/example.fa
#coverage of 10%
#SE sequencing
#average length of DNA fragments is 30~40 with the shortest length of 20. And the max length could be reach to about 200bp.
 |     ..
 |    .  .
 |   .     .
 |  .         .
 |               .
 |                  .
 |                      .
 --------------------------->
#overhang length of DNA double strain is 3bp to 20bp of random distribution.
#Single strain library preparation is adopted
#Add damage by nucleotide substitution according to training stochastic matrix.
#PCR duplication, ~1000x.
#lost rate single strain is 50%.
#efficiency of fill-in reaction is 50%
#WARNING:MAKE SURE THAT YOU HAVE ENOUGH MEMARY SPACE MORE THAN ReferenceSize*Coverage*Duplication*3
 perl Psim.pl illumina --damage --ref ../example/example.fa --cov 0.1 --fragmean 35 --fragsd 48 --fraglim 20+ --overhang 3_20 --library 1 --lost 0.5 --ampmean 1000 --mutarray ../example/mutarray_example.txt --mutsite ../example/mutsite_example.txt --effic 0.5 --dir ../output/palo/

#RNA-SEQ
#Illumina PE sequencing
#randomly generate snp with the rate of 0.001
 perl Psim.pl illumina --ref ../example/example.fa --cdna --gff ../example/gff_example --covcdna 10 --snp 0.001 --dir ../output/
======================================EXAMPLE END============================================
";
#============================CHECK PARAMETER===============================
die color ("bold magenta"),$USAGE,color("reset"),"\n" if(!$Command || $Help);
die $ExampleInfo if($Example || $Command=~/example/i);
die "Wrong para \n" if(!($Command=~/example/i || $Command=~/illumina/i || $Command=~/roche/i || $Command=~/variation/i || $Command=~/solid/i));
$ploidy||=1;
die "Wrong ploidy parameter!\n" if($ploidy<1);
$ploidy=int($ploidy);
$SNP=0.001 if(!defined $SNP);
$SV="0.1:3000" if(!defined $SV);
my $SNPReport="SNPReportByPsim.txt";
my $SVReport="SVReportByPsim.txt";
$Coverage||=3;
if(defined $Circle){$Circle=1;}
else{$Circle=0}
$Error="0.0005" if(!defined $Error);
die "Wrong SNP parameter \n$USAGE" if((!($SNP=~/\d+\.?\d*/) && !(-e $SNP)) || ($SNP=~/\d+\.\d*/ && $SNP>1));
die "Wrong SV parameter \n$USAGE" if( !(-e $SV) && (!$SV=~/\d*\.?\d+\:\d+/));
$Qtype||="!";
$Qmean||=37;
$Qsd||=1;
my ($l1,$l2);
if($Command=~/variation/i)
{
	die color ("bold magenta"),$variation,color("reset"),"\n" if(!$Reference || $Help);
	die "Cannot open reference sequence file\n" if(!(-e $Reference));
	$Output||="VariationSequenceByPsim";
}
#######130917-18-59-here
die "Wrong duplication percent parameter. It should be less than 1\n" if($duplpercent && $duplpercent>=1);
die "duplication power law parameters should be seoarate by \"\,\" and a\<1 and b\<0\n" if($duplpara && (!($duplpara=~/,/)  || (split /,/,$duplpara)[0]>=1 || (split /,/,$duplpara)[1]>=0));

if($Damage)
{
	die color("red"),"\nOooopssssssssss\nAs to my point, Roche 454 PE sequencing technoly is not suitable for palogenome sequencing!\nPalogenome are sheared into short pieces after long time digestion.\nHow do you cut the gel while library preparation???\nSorry, I don't wanna simulate this process~\n\n",color("reset") if($Command=~/roche/i && $PE);
	$FragmentMean||=55;
	$FragmentSD||=28;
	$FragLim||="27+";
	$Mismatch||="3_20";
	$Efficiency||=0.5;
	if(defined $MutationSite)
	{
		open MUTATIONSITE,"<$MutationSite" || die "cannot open mutation rate along site file $MutationSite\n";
		my $mutationsite=<MUTATIONSITE>;
		chomp $mutationsite;
		my @rate=split /\,/,$mutationsite;
		for my $i(0..$#rate)
		{
			$MutationRate{$i}=$rate[$i];
		}
		close MUTATIONSITE;
	}
	else
	{
		$MutationRate{0}=0.01;
	}
	if(defined $MutationArray)
	{
		open MUTATIONARRAY,"<$MutationArray" || die "cannot open mutation array file $MutationArray\n";
		while(<MUTATIONARRAY>)
		{
			if(!/^#/)
			{
				chomp;
				my @array=split /\;/;
				foreach(@array)
				{
					my @mutinfo=split /\,/;
					$MutationArray{$mutinfo[0]}{$mutinfo[1]}=$mutinfo[2];
				}
				last;
			}
		}
		close MUTATIONARRAY;
	}
	else
	{
		$MutationArray{"C"}{"T"}=1;
	}
	$Library||=1;
	$AmpMean||=850;
	$AmpSD||=0.0212*$AmpMean+10.974;
	if($Library==1)
	{
		$LostSingle||=0.5;
	}
	else
	{
		$DoublePercent||=0.9;
	}
}
if($Command=~/illumina/i)
{
	die color ("cyan"),$illumina,color("reset"),"\n" if(!$Reference || $Help);
	die "Cannot open reference sequence file\n" if(!(-e $Reference));
	$FragmentMean||=200;
	$FragmentSD||=10;
	$ReadLeng||=100;
	$Adapter||="AAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT:GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGCAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC";
	die "Wrong adapter sequence \n$USAGE" if(!$Adapter=~/\w+\:?\w*/);
	($adapter1,$adapter2)=split /\:/,$Adapter;
	$l1=length($adapter1);
	$l2=length($adapter2);
	$Output||="Psim_Illumina";
}
elsif($Command=~/roche/i)
{
	die color ("blue"),$roche,color("reset"),"\n" if(!$Reference || $Help);
	die "Cannot open reference sequence file\n" if(!(-e $Reference));
	$Insert||="8000:30" if($PE);
	$FragmentMean||=450;
	$Linker||="ATAACTTCGTATAATGTATGCTATACGAAGTTAT";
	$FragmentSD||=50;
	$Output||="Psim_Roche454";
}
elsif($Command=~/solid/i)
{
	die color ("green"),$solid,color("reset"),"\n" if(!$Reference || $Help);
	die "Cannot open reference sequence file\n" if(!(-e $Reference));
	$FragmentMean||=100;
	$FragmentSD||=10;
	$ReadLeng||=50;
	$Output||="Psim_SOLiD";
	$Header||="G";
}
$Dir||="./";
die  "$Dir do not exist!\nplease make directory $Dir first\n" if(!(-d $Dir));
my $flowcell=1;
my $titlenumber=1;
my $n=0;

my %DigestRegion;  #Two-dimensional hash, key1--chromosome, key2--transcript ID, value--strand:start site:end site
my @rna_poisson;   #each RNA sequence sequencing depth
my @cDNAFragment=qw();
if($CDNA)
{
	$CovcDNA||=10;
	&ReadGff("gff",$GFF,"info",\%DigestRegion);
	my $NumTrans=scalar(keys %DigestRegion);
	die "NO TRANSCRIPT INFORMATION IN $GFF!\n" if($NumTrans==0);
	@rna_poisson=random_poisson($NumTrans,$CovcDNA);
	$Circle=0;
}
if(($SNP ne 0) || ($SV ne 0))
{
	open NEWSEQ,">$Dir"."NewReferenceSequence.fa";
}

my $info="Simulation Sequencing Parameters
===========================================================
reference data $Reference\n";
if($Circle eq 1){$info.="circle reference\n"}
else{$info.="liner reference\n"}
$info.="simulate $ploidy\-ploid with variation of\n";
if($SNP eq 0){$info.="    no snp sites\n"}
elsif($SNP=~/\d+\.?\d*/){$info.="    snp rate $SNP\n"}
else{$info.="    snp rate file $SNP\n"}
if($SV eq 0){$info.="    no sv\n"}
elsif($SV=~/\:/ ){$info.="    sv rate and average length $SV\n"}
else{$info.="    sv file $SV\n"}
$info.="coverage of sequencing data is $Coverage\n";
if($Command=~/illumina/i && $PE){$info.="Illumina PE sequencing\nadapter sequence is $Adapter\n";}
elsif($Command=~/illumina/i && !$PE){$info.="Illumina SE sequencing\nadapter sequence is $Adapter\n";}
elsif($Command=~/roche/i && $PE){$info.="Roche PE sequencing\nlinker sequence is $Linker\n";}
elsif($Command=~/roche/i && !$PE){$info.="Roche PE sequencing\n";}
elsif($Command=~/solid/i && !$PE){$info.="SOLiD SE sequencing\nheader base is $Header\n";}
else{$info.="SOLiD PE sequencing\n";}
$info.="average fragment length is $FragmentMean\nsd of fragment length is $FragmentSD\n";
$info.="fragment length limit is $FragLim\n" if(defined $FragLim);
$info.="read length is $ReadLeng\n" if($Command=~/illumina/i || $Command=~/solid/i);
if($Error ne 0){$info.="sequencing error rate is $Error\n"}
if($Damage){$info.="simulation palogenome\noverhang range is $Mismatch\n$Library stain library preparation method\n"}
if($DoublePercent){$info.="double strain both kept rate is $DoublePercent\n"}
elsif($LostSingle){$info.="single strain lost rate is $LostSingle\n\n";}
if($CDNA){$info.="RNA-seq\ninput gff file is $GFF\nmean coverage is $CovcDNA\n\n";}
$info.="output directory is $Dir\noutput sequencing fiel prefix is $Output
============================================================\n";
print color("yellow"),$info,color("reset") if($Reference && $Command && !$Help);

#============================CHECK PARAMETER END===============================

#============================HANDLE INFORMATION================================

if(!($Command=~/variation/i))   #illumina: fq file format, roche & solid: fa $ qual
{
	open OUT1,">$Dir$Output.fasta" if(!($Command=~/illumina/i));
	open OUT2,">$Dir$Output.qual" if(!($Command=~/illumina/i));
	open OUT1,">$Dir$Output.fastq" if($Command=~/illumina/i && !$PE);
	open OUT1,">$Dir$Output\-1.fastq" if($Command=~/illumina/i && $PE);
	open OUT2,">$Dir$Output\-2.fastq" if($Command=~/illumina/i && $PE);
	open NAME,">$Dir"."NameRecordByPsim.txt";      #relationship of reads name and reads information
}

my (@SV,@SNP,@MUTATION,@NAME,@ERROR);
if($SNP ne 0)
{
	open SNPOUT,">$Dir"."SNPReportByPsim.txt";
}
if($SV ne 0)
{
	open SVOUT,">$Dir"."SVReportByPsim.txt";
}
if( !($Command=~/variation/i) && (defined $Error) && ($Error ne 0))
{
	open ERROR,">$Dir"."SequencingErrorByPsim.txt";
}
if($MutationSite)
{
	open MUTATION,">$Dir"."MutationByPsim.txt";
}

#==============================MAIN PROGRAM===============================
my @Qualitys;
	my $Sequence="";
	my $SeqName="";
local $/=">";  #v1.6
my $m;  #count from 1 to ploidy number
open REFERENCE,"<$Reference";
<REFERENCE>;
while(<REFERENCE>)
{
	next if(/^>$/);
	s/\>//;
	my @tmp=split /\n/,$_;
	map {chomp $_} @tmp;
	$SeqName=shift @tmp;
	$Sequence=join "",@tmp;
	if($Sequence ne "" && length($Sequence)>0)
	{
		for($m=1;$m<=$ploidy;$m++)
		{
			$/="\n";
			&Main($Sequence);
			$/=">";
		}
	}
}
close REFERENCE;
$Sequence="";    #Release memory
$/="\n";

#=============================MAIN SUBROUTINE===============================
sub Main
{
#STEP ONE:	SNP generating
	my $SEQUENCE=shift;
	if($SNP ne 0)
	{
		my $snptype=0;
		#if $SNP is a decimal==> set snp rate, otherwise snp information from file
		$snptype=1 if(!($SNP=~/\d\.\d*/));
		&SNP("seq",\$SEQUENCE,"snptype",$snptype,"seqname",$SeqName,"snp",$SNP,"snpinfo",\@SNP);
		print scalar(@SNP),"\n";
		map{print "$_\n"} @SNP;
	}
#RNA-seq sequencing
	if($CDNA)
	{
		my $n=0;  #used to get sequencing depth from array @rna_poisson
		my @name=split /\s+/,$SeqName;
		my $CHR=$name[0]; #in GFF file, the sequence name should be the content before the first space of fa title
		if(exists $DigestRegion{$CHR})
		{
			foreach my $ID(keys %{$DigestRegion{$CHR}})
			{
				my $transcript="";
				my @a=split /\:/,$DigestRegion{$CHR}{$ID};
				my $Strand=$a[0];
				my @se=split /\;/,$a[1];
				if($Strand eq "+")
				{
					foreach(@se)
					{
						my ($start,$end)=split /\,/;
						$transcript.=substr($SEQUENCE,$start,($end-$start+1));
					}
				}
				elsif($Strand eq "-")
				{
					foreach(@se)
					{
						my ($start,$end)=split /\,/;
						$transcript=substr($SEQUENCE,$start,($end-$start+1)).$transcript;
						$transcript=&ReverseComplement($transcript);
					}
				}
				&Library(\$transcript,$rna_poisson[$n]);
				$n++;
			}
		}
		else
		{
			die "Wrong Reference Name format\ndo not exists annotation information in gff file\n";
		}
	}
	else
	{
#STEP TWO:	SV
		if($SV ne 0)
		{
			my $SVType=1;
			$SVType=0 if($SV=~/\d*\.?\d+\:\d+/);
			$SEQUENCE=&StructuralVariation("seq",\$SEQUENCE,"name",$SeqName,"type",$SVType,"sv",$SV,"report",\@SV);
		}

#STEP THREE:	LIBRARY PREPARATION AND SEQUENCING
#====LIBRARY PREPARATION
		&Library(\$SEQUENCE,$Coverage) if(!($Command=~/variation/));
	}
	if($Command=~/variation/ || ((($SNP ne 0) || ($SV ne 0))  && $print))
	{
		print NEWSEQ ">$SeqName-$m\n";
		print NEWSEQ $SEQUENCE."\n";
	}
	if($SNP ne 0)
	{
		print SNPOUT "========$SeqName-$m snp information========\n";
		map {print SNPOUT "$_\n"} @SNP;
		@SNP=qw//;
		print SNPOUT "\n"
	}
	if($SV ne 0)
	{
		print SVOUT "========$SeqName-$m sv information========\n";
		map {print SVOUT "$_\n"} @SV;
		@SV=qw//;
		print SVOUT "\n";
	}
	if( !($Command=~/variation/i) && (defined $Error) && ($Error ne 0))
	{
		print ERROR "========$SeqName-$m sequencing error information========\n";
		map {print ERROR "$_\n"} @ERROR;
		@ERROR=qw//;
		print ERROR "\n";
	}
	if($MutationSite)
	{
		print MUTATION "========$SeqName-$m sequencing error information========\n";
		print MUTATION "chrom\tread\tsite_of_frag\tori\tmut\n";
		map {print MUTATION "$_\n"} @MUTATION;
		@MUTATION=qw//;
		print MUTATION "\n";
	}
}
sub Library
{
	my ($Sequence,$Coverage)=@_;
	my $ltype=0;
	my @Fragment;
	my ($insert,$insertsd);
	if($Command=~/roche/i && $PE)
	{
		$ltype=$Insert.":".length($Linker);
		($insert,$insertsd)=split /\:/,$Insert;
	}
	my $readmean;
	if($Command=~/roche/i)
	{
		$readmean=$FragmentMean;
	}
	else
	{
		$readmean=$ReadLeng;
	}
	my @StartLeng;
	if($Damage)
	{
		my $DamageParaA=1.5;
		my $DamageParaK=1.2;
		my $DamageParaMu=15;
		@StartLeng=&PieceGenerate("WL",length($$Sequence),"circle",$Circle,"Coverage",$Coverage,"fragmean",$FragmentMean,"fragsd",$FragmentSD,"fraglim",$FragLim,"dpa",$DamageParaA,"dpk",$DamageParaK,"dpmu",$DamageParaMu);
	}
	elsif($duplpercent && $duplpara && $duplpara=~/,/)
	{
		$FragmentMean||=200;
		$FragmentSD||=10;
		my ($a,$b)=split/,/,$duplpara;
		my @StartLeng=&PieceGenerate("WL",length($$Sequence),"type",$ltype,"circle",$Circle,"Coverage",$Coverage,"readmean",$readmean,"fragmean",$FragmentMean,"fragsd",$FragmentSD,"fraglim",$FragLim,"duplp",$duplpercent,"dupll",$dupllimit,"dupla",$a,"duplb",$b);
	}
	else
	{
		$FragmentMean||=200;
		$FragmentSD||=10;
		@StartLeng=&PieceGenerate("WL",length($$Sequence),"type",$ltype,"circle",$Circle,"Coverage",$Coverage,"readmean",$readmean,"fragmean",$FragmentMean,"fragsd",$FragmentSD,"fraglim",$FragLim,"duplp",$duplpercent,"dupll",$dupllimit);
	}
	$$Sequence.=substr($$Sequence,0,int($FragmentMean*2));
	$$Sequence.=substr($$Sequence,0,int($insert*2)) if($Command=~/roche/i && $PE);

#====DAMAGE SEQUENCING
	if($Damage)
	{
		my $lost=$LostSingle;
		$lost=$DoublePercent if($Library==2);
		my @PaloName;   #keep information of each fragment
		&Overhanging($Library,$Sequence,\@StartLeng,$Mismatch,$SeqName,$lost,\%MutationRate,\%MutationArray,\@Fragment,\@PaloName,\@MUTATION,$Efficiency);
		my (@Library,@LibraryName);
#@Fragment--DNA library before amplify, @Library--DNA library after amplify
#if the library no need to amplify, @Library=@Fragment
		if($noamp)   
		{
			@Library=@Fragment;
			@LibraryName=@PaloName;
		}
		else   #PCR amplification
		{
			&PCRduplication(\@Fragment,\@Library,\@PaloName,\@LibraryName,$AmpMean,$AmpSD)
		}
		my $ReadsNum=scalar(@Library);
		$ReadsNum=2*$ReadsNum if($PE);
		&Quality("cmd",$Command,"array",\@Qualitys,"readleng",$ReadLeng,"mean",$Qmean,"sd",$Qsd,"num",$ReadsNum,"type",$Qtype);
#====SEQUENCING
		for my $i(0..$#Library)
		{
			if($Command=~/illumina/i && !$PE){
				&SequencingIlluminaSE($Library[$i],$LibraryName[$i],\$flowcell,\$titlenumber,\$n)}
			elsif($Command=~/illumina/i && $PE){
				&SequencingIlluminaPE($Library[$i],$LibraryName[$i],\$flowcell,\$titlenumber,\$n)}
			elsif($Command=~/roche/i){
				&SequencingRoche($Library[$i],$LibraryName[$i],\$titlenumber)}
			elsif($Command=~/solid/i && !$PE){
				&SequencingSOLiDSE($Library[$i],$LibraryName[$i],\$flowcell,\$titlenumber,\$n)}
			elsif($Command=~/solid/i && $PE){
				&SequencingSOLiDPE($Library[$i],$LibraryName[$i],\$flowcell,\$titlenumber,\$n)}
		}
	}
	else
	{
#====QUALITY
		if($Command=~/illumina/i || $Command=~/solid/i)
		{
			my $ReadsNum=scalar(@StartLeng);
			$ReadsNum=2*$ReadsNum if($PE);
#&Quality($Command,\@Qualitys,$ReadLeng,$Qmean,$Qsd,$ReadsNum,$Qtype);
			&Quality("cmd",$Command,"array",\@Qualitys,"readleng",$ReadLeng,"mean",$Qmean,"sd",$Qsd,"num",$ReadsNum,"type",$Qtype);
		}
		else
		{
			my @ReadsLeng;
			map {my @temp=split /\_/;push @ReadsLeng,$temp[-1]} @StartLeng;
			&Quality("cmd",$Command,"array",\@Qualitys,"readleng",\@ReadsLeng,"mean",$Qmean,"sd",$Qsd);
		}
##====GENERATE LIBRARY
		for my $i(0..$#StartLeng)
		{
			my $sequence;
			if($Command=~/roche/i && $PE)
			{
				#push @StartLeng,"3\_$StartPoint\_$InsertEnd\_$length[$i]";
				#push @StartLeng,"0\_$StartPoint\_$InsertEnd\_$InsertSite\_$length[$i]";
				#push @StartLeng,"5\_$InsertSite\_$end\_$length[$i]";
				my @temp_roche=split /\_/,$StartLeng[$i];
				if($temp_roche[0]==3)
				{
					my $orileng=$temp_roche[2]-$temp_roche[1]+1;
					my $linkeradd=$temp_roche[3]-$orileng;
					$sequence=substr($$Sequence,$temp_roche[1],$orileng).substr($Linker,0,$linkeradd);
					my $readinfo="$SeqName\-$m\-3end ".($i+1)." start=$temp_roche[1] end=$temp_roche[2] fragmentlength=$temp_roche[3]";
##====SEQUENCING
					&SequencingRoche($sequence,$readinfo,\$titlenumber) if($Command=~/roche/i);
				}
				elsif($temp_roche[0]==0)
				{
					my $ori1=$temp_roche[2]-$temp_roche[1]+1;
					my $ori2=$temp_roche[4]-length($Linker)-$ori1;
					$sequence=substr($$Sequence,$temp_roche[1],$ori1);
					$sequence.=$Linker;
					$sequence.=substr($$Sequence,$temp_roche[3],$ori2);
					my $readinfo="$SeqName\-$m\-bothPE ".($i+1)." start=$temp_roche[3] end=$temp_roche[2] fragmentlength=$temp_roche[4]";
##====SEQUENCING
					&SequencingRoche($sequence,$readinfo,\$titlenumber) if($Command=~/roche/i);
				}
				else
				{
					my $orileng=$temp_roche[2]-$temp_roche[1]+1;
					my $linkeradd=$orileng-$temp_roche[3];
					$sequence=substr($Linker,$linkeradd).substr($$Sequence,$temp_roche[1],$orileng);
					my $readinfo="$SeqName-$m-5end ".($i+1)." start=$temp_roche[1] end=$temp_roche[2] fragmentlength=$temp_roche[3]";
##====SEQUENCING
					&SequencingRoche($sequence,$readinfo,\$titlenumber) if($Command=~/roche/i);
				}
			}
##====SEQUENCING
			else
			{
				my ($start,$length)=split /\_/,$StartLeng[$i];
				$sequence=substr($$Sequence,$start,$length);
				my $readinfo="$SeqName-$m ".($i+1)." start=$start end=".($start+$length-1)." fragmentlength=$length";
				if($Command=~/illumina/i && !$PE){
					&SequencingIlluminaSE($sequence,$readinfo,\$flowcell,\$titlenumber,\$n)}
				elsif($Command=~/illumina/i && $PE){
					&SequencingIlluminaPE($sequence,$readinfo,\$flowcell,\$titlenumber,\$n)}
				elsif($Command=~/roche/i){
					&SequencingRoche($sequence,$readinfo,\$titlenumber)}
				elsif($Command=~/solid/i && !$PE){
					&SequencingSOLiDSE($sequence,$readinfo,\$flowcell,\$titlenumber,\$n)}
				elsif($Command=~/solid/i && $PE){
					&SequencingSOLiDPE($sequence,$readinfo,\$flowcell,\$titlenumber,\$n)}
			}
		}
	}
}

my %xy;
sub SequencingIlluminaSE
{
	my ($seq,$readinfo,$f,$t,$n)=@_;
	$seq.=$adapter2;
	my $x=int(rand(20000))+1;
	my $y=int(rand(21500))+1;
#@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
	my $read=substr($seq,0,$ReadLeng);
	if(!($read=~/N/i))
	{
		my $prefix="\@RachelWu\:1\:BENM\:$$f\:$$t\:$x\:$y 0\:N\:0\:ATCACG";
		print NAME "$prefix\t$readinfo\n";

		$read=&SequencingError($read,$Error,$prefix,\@ERROR);

		print OUT1 "$prefix\n";
		print OUT1	"$read\n\+\n";
		print OUT1 "$Qualitys[$count]\n";

		$$n+=1;
		if($$n>400)
		{
			$$t+=1;
			$$n=1;
		}
		$$f+=1 if($$t>120);
	}
	$count+=1;
}

sub SequencingIlluminaPE
{
	my ($reads,$readinfo,$f,$t,$n)=@_;
	my $readseq1=$reads.$adapter2;
	my $readseq2=&ReverseComplement($adapter1.$reads);

	my $x=int(rand(20000))+1;
	my $y=int(rand(21500))+1;
	#@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
	my $prefix="\@RachelWu\:1\:BENM\:$$f\:$$t\:$x\:$y";
	my $pe1=substr($readseq1,0,$ReadLeng);
	my $pe2=substr($readseq2,0,$ReadLeng);

	if(!($pe1=~/N/i))
	{
		$pe1=&SequencingError($pe1,$Error,$prefix." 1\:N\:0\:ATCACG",\@ERROR);
		print OUT1 "$prefix"." 1\:N\:0\:ATCACG\n$pe1\n";
		print OUT1 "+\n"."$Qualitys[$count]\n";

	}
	if(!($pe2=~/N/i))
	{
		$pe2=&SequencingError($pe2,$Error,$prefix." 2\:N\:0\:ATCACG",\@ERROR);
		print OUT2 "$prefix"." 2\:N\:0\:ATCACG\n$pe2\n";
		print OUT2 "+\n"."$Qualitys[$count+1]\n";
	}
	if(!($pe1=~/N/i && $pe2=~/N/i))
	{
		print NAME "$prefix\t$readinfo\n";
		$$n+=1;
		if($$n>400)
		{
			$$t+=1;
			$$n=1;
		}
		$$f+=1 if($$t>120);
	}
	$count+=2;
}

sub SequencingRoche
{
	my ($reads,$readinfo,$n)=@_;
	if(!($reads=~/N/i))
	{
		my $prefix="\>Psim_Roche454\.$$n";
		print NAME "$prefix\t$readinfo\n";

		&SequencingError($reads,$Error,$prefix,\@ERROR);
		print OUT1 "$prefix\n$reads\n";
		print OUT2 "$prefix\n$Qualitys[$$n-1]\n";
	}
	$$n+=1;
}

sub SequencingSOLiDSE 
{
	my ($seq,$readinfo,$f,$t,$n)=@_;
	my $read=substr($seq,0,$ReadLeng);
	if(!($read=~/N/i))
	{
		my $prefix="\>RachelWu\_$$f\_$$t\_$$n\_F3";
		print NAME "$prefix\t$readinfo\n";
		$read=&SequencingError($read,$Error,$prefix,\@ERROR);
		&Base2Col(\$read,$Header);
		print OUT1 "$prefix\n$read\n";
		print OUT2 "$prefix\n$Qualitys[$count]\n";
		$$n+=1;
		if($$n>400)
		{
			$$t+=1;
			$$n=1;
		}
	}
	$count+=1;
}
sub SequencingSOLiDPE
{
	my ($reads,$readinfo,$f,$t,$n)=@_;
	my $readseq1=$reads;
	my $readseq2=&ReverseComplement($reads);

	my $pe1=substr($readseq1,0,$ReadLeng);
	my $pe2=substr($readseq2,0,$ReadLeng);
	my $prefix="\>RachelWu\_$$f\_$$t\_$$n";
	if(!($pe1=~/N/i))
	{
		$pe1=&SequencingError($pe1,$Error,$prefix."\_F3",\@ERROR);
		&Base2Col(\$pe1,$Header);
		print OUT1 "$prefix"."\_F3\n$pe1\n";
		print OUT2 "$prefix"."\_F3\n$Qualitys[$count]\n";
	}
	if(!($pe2=~/N/i))
	{
		$pe2=&SequencingError($pe2,$Error,$prefix."\_R3",\@ERROR);
		&Base2Col(\$pe2,$Header);

		print OUT1 "$prefix"."\_R3\n$pe2\n";
		print OUT2 "$prefix"."\_R3\n$Qualitys[$count+1]\n";
	}
	if(!($pe1=~/N/i && $pe2=~/N/i))
	{
		print NAME "$prefix\t$readinfo\n";
		$$n+=1;
		if($$n>400)
		{
			$$t+=1;
			$$n=1;
		}
	}
	$count+=2;
}
close NEWSEQ;
close NAME;
close OUT1;
close OUT2 || print "Simulation Done\n";
close SNPOUT || print "";
close SVOUT || print "";
close ERROR|| print "";
close MUTATION || print "";
