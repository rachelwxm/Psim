#!/usr/bin/perl -w
use strict;
use Math::Random qw(random_poisson);
use Getopt::Long;
use Term::ANSIColor;
use FindBin qw($Bin $Script);
use lib "$FindBin::Bin/lib";
use 5.010;
use constant PI=>3.14159265;
use Data::Dumper;
use File::Basename qw(basename dirname);
use MismatchDoubleStrain;
use normal;
use PieceGenerate;
use ReverseComplement;
use SNP;
use StructuralVariation;
use Base2Col;
use Quality;
use SequencingError;
use ReadGff;

#======================DEFAULT PARAMETERS========================
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

my @code=([0,1,2,3],[1,0,3,2],[2,3,0,1],[3,2,1,0]);
my @bases=qw/A C G T/;
my @othercode=(".",4,5,6);
my %colcode=();
my %decode=();
foreach my $i(0..3)
{
	foreach my $j(0..3)
	{
		$decode{$code[$i]->[$j]}->{$bases[$i]}=$bases[$j];
	}
}
foreach my $i(0..3)
{
	foreach my $j(0..3)
	{
		$colcode{"$bases[$i]$bases[$j]"}=$code[$i]->[$j];
	}
}

#==========================USER PARAMETERS==============================
my ($Help,$Command,$Reference,$Circle,$SNP,$SV,$Coverage,$PE,$FragmentMean,$FragmentSD,$FragLim,$ReadLeng,$ReadSD,$LostSingle,$DoublePercent,$Mismatch,$Palo,$Adapter,$Output,$Dir,$Library,$Mutation);
my ($Insert,$Linker,$Error,$adapter1,$adapter2,$Qmean,$Qsd,$Header,$Example,$GFF,$CDNA,$CovcDNA,$Qtype);
my $count=0;
my $SNPReport="SNPReportByPsim.txt";
my $SVReport="SVReportByPsim.txt";
$Command=$ARGV[0];
GetOptions(
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
	"palo"		=>\$Palo,
	"mismatch:s"=>\$Mismatch,
	"library:s"	=>\$Library,
	"mutation:s"=>\$Mutation,
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
	);

my $USAGE="
       =============================================================
      || WARNINGS:                                                 ||
      ||     YOU SHOLD HAVE ENOUGH FREE SPACE AND MEMORY!          ||
      ||     if input reference sequence size is X, you shold have ||
      || at least 3X free space and X memory.                      ||
       =============================================================
USAGE:   Used to generate NGS data. The sequencing platform include Illumina
         Roche 454 and SOLiD. 
         Input reference sequence could be one seq or multiple seqs. And the 
         sequencing sample can have mismatch end like palogenome.

Author:  RachelWu\(rachelwu123\@gmail.com\)\, BENM\(BinxiaoFeng\@gmail.com\)
Version: v1.0

Run:     perl $0 [command] [options]

command: illumina   generate Illumina sequencing data
         roche      generate Roche 454 sequencing data
         solid      generate SOLiD sequencing data
         example    show running examples
****************************************************************************
\n";

my $illumina="
*******************************************************************************************************
USAGE:  perl $0 illumina [options]
OPTIONS:
    --ref <FILE>          reference sequence file
    --circle              the reference sequence is in circle
    --snp <NUM|FILE>      random snp rate(default=0.001) or specific snp site and rate
                          file(format refer to ../example/snp_example.txt)
    --sv <NUM|FILE>       structure variction rate and average length(default=0.1:3000)
                          or specific sv type and site file(format refer to ../example/sv_example.txt)
    --cov <NUM>           sequencing coverage of reference sequence(s)(default=3)
    --pe                  PE sequencing or not
    --fragmean <INT>      average length of library fragment(default=200)
    --fragsd <INT>        standard deviation of library fragment(default=10)
    --fraglim <INT>       limit length of fragment library(\"20+\" means must above 20nt, and
                          \"240-\" means must shorter than 240nt,if(-palo) this para default=20+)
    --read <INT>          average length of reads(default=100)
    --adapter <CHAR>      adapter sequence(default seq refer to ../example/adapterIllumina_example.txt, 
                          split two adapters by \"\:\")
    --error <NUM>         sequencing error rate and possiblities of single base error, insert and deletion
                          (default=\"0.0005:0.34:0.33:0.33\")
    --qtype <CHAR>        quality type(\!=Sanger format, \@=illumina 1.3~1.8- format, default=\!)
    --qmean <NUM>         peak value of quality score(default=37)
    --qsd <NUM>           standard deciation of quality score(default=1)

  **palogenome sequencing parameters**
    --palo                the sequencing sample have mismatch and injured end
    --mismatch <INT>      mismatch length range of double strain(default=\"3\_20\")
    --library <INT>       library preparation type(1=single strain,2=double strain,default=1)
    --mutation <NUM>      mutation rate of methylation and deamination at mismatch end(default=0.01)
    --ds <NUM>            percent of two strain both be sequenced while library=2(default=0.9)
    --lost <NUM>          lost rate single strain while library=1(default=0.5)

  **RNA-SEQ sequencing parameters**
    --cdna                RNA-SEQ sequencing
    --gff <FILE>          input gene gff file
    --covcdna <INT>       mean coverage of cDNA sequences(default=10 if --cdna)

  **output parameters**
    --dir <CHAR>          output directory(default=\".\/\")
    --output <CHAR>       prefix of output file(default=\"SimulateIlluminaSequencingByPsim\")
*******************************************************************************************************
";

my $roche="
*******************************************************************************************************
USAGE:  perl $0 roche [options]
OPTIONS:
    --ref <FILE>          reference sequence file
    --circle              the reference sequence is in circle
    --snp <NUM|FILE>      random snp rate(default=0.001) or specific snp site and rate
                          file(format refer to ../example/snp_example.txt)
    --sv <NUM|FILE>       structure variction rate and average length(default=0.1:3000)
                          or specific sv type and site file(format refer to ../example/sv_example.txt)
    --cov <NUM>           sequencing coverage of reference sequence(s)(default=3)
    --pe                  PE sequencing or not
    --insert <INT>        average insert size and sd(default=8000:30 if --pe)
    --linker <CHAR>       PE insert seq(default=\"ATAACTTCGTATAATGTATGCTATACGAAGTTAT\")
    --fragmean <INT>      average length of library fragment(equal to reads length,default=450)
    --fragsd <INT>        standard deviation of library fragment(default=50)
    --fraglim <INT>       limit length of fragment library(\"20+\" means must above 20nt, and
                          \"800-\" means must shorter than 800nt,if(-palo) this para default=20+)
    --error <NUM>         sequencing error rate and possiblities of single base error, insert and deletion
                          (default=\"0.0005:0.34:0.33:0.33\")
    --qmean <NUM>         peak value of quality score(default=37)
    --qsd <NUM>           standard deciation of quality score(default=1)

  **palogenome sequencing parameters**
    --palo                the sequencing sample have mismatch end
    --mismatch <INT>      mismatch length range of double strain(default=\"3\_20\")
    --library <INT>       library preparation type(1=single strain,2=double strain,default=1)
    --mutation <NUM>      mutation rate of methylation and deamination at mismatch end(default=0.01)
    --ds <NUM>            percent of two strain both be sequenced while library=2(default=0.9)
    --lost <NUM>          lost rate single strain while library=1(default=0.5)

  **RNA-SEQ sequencing parameters**
    --cdna                RNA-SEQ sequencing
    --gff <FILE>          input gene gff file
    --covcdna <INT>       mean coverage of cDNA sequences(default=10 if --cdna)

  **output parameters**
    --dir <CHAR>          output directory(default=\".\/\")
    --output <CHAR>       prefix of output file(default=\"Simulate454SequencingByPsim\")
*******************************************************************************************************
";

my $solid="
*******************************************************************************************************
USAGE:  perl $0 solid [options]
OPTIONS:
    --ref <FILE>          reference sequence file
    --circle              the reference sequence is in circle
    --snp <NUM|FILE>      random snp rate(default=0.001) or specific snp site and rate
                          file(format refer to ../example/snp_example.txt)
    --sv <NUM|FILE>       structure variction rate and average length(default=0.1:3000)
                          or specific sv type and site file(format refer to ../example/sv_example.txt)
    --cov <NUM>           sequencing coverage of reference sequence(s)(default=3)
    --pe                  PE sequencing or not
    --fragmean <INT>      average length of library fragment(default=100)
    --fragsd <INT>        standard deviation of library fragment(default=10)
    --fraglim <INT>       limit length of fragment library(\"20+\" means must above 20nt, and
                          \"240-\" means must shorter than 240nt,if(-palo) this para default=20+)
    --read <INT>          average length of reads(default=50)
    --error <NUM>         sequencing error rate and possiblities of single base error, insert and deletion
                          (default=\"0.0005:0.34:0.33:0.33\")
    --qmean <NUM>         peak value of quality score(default=37)
    --qsd <NUM>           standard deciation of quality score(default=1)
    --header <CHAR>       sequencing header base(default=G)

  **palogenome sequencing parameters**
    --palo                the sequencing sample have mismatch end
    --mismatch <INT>      mismatch length range of double strain(default=\"3\_20\")
    --library <INT>       library preparation type(1=single strain,2=double strain,default=1)
    --mutation <NUM>      mutation rate of methylation and deamination at mismatch end(default=0.01)
    --ds <NUM>            percent of two strain both be sequenced while library=2(default=0.9)
    --lost <NUM>          lost rate single strain while library=1(default=0.5)

  **RNA-SEQ sequencing parameters**
    --cdna                RNA-SEQ sequencing
    --gff <FILE>          input gene gff file
    --covcdna <INT>       mean coverage of cDNA sequences(default=10 if --cdna)

  **output parameters**
    --dir <CHAR>          output directory(default=\".\/\")
    --output <CHAR>       prefix of output file(default=\"SimulateSOLiDSequencingByPsim\")
*******************************************************************************************************
";

my $ExampleInfo="======================================EXAMPLE================================================
#generate a random reference file length of 10000 bp
perl GenerateSeq.pl  -l 10000  > MyNewGenerateRefFile.fa\n
#generate several(assume 3) references with average length of 10000 bp, and seq pre is \"example\"
perl GenerateSeq.pl -n 3  -l 10000 -auto -pre example > MyNewGenerateRefFile.fa\n
#generate several(assume 3) references with length of 10000 bp, 20000 bp and 16000 bp
perl GenerateSeq.pl -n 3  -l 10000,20000,16000 > MyNewGenerateRefFile.fa\n
#generate illumina PE sequencing data using default parameter
perl Psim.pl illumina --ref ../example/example.fa --pe \n
#generate illumina PE palogenome sequencing data with set parameters
perl Psim.pl illumina --snp 0.03 --sv 0.3 --cov 0.01 --pe --fragmean 35 --fragsd 48 --fraglim 20+ 
--read 100 --palo --mismatch 3_20 --library 1 --mutation 0.01 --ds 0.9 --lost 0.5 --dir ../output
--output palo_sim_test
======================================EXAMPLE END============================================\n
";
#============================CHECK PARAMETER===============================
die color ("bold magenta"),$USAGE,color("reset"),"\n" if(!$Command || $Help);
die $ExampleInfo if($Example || $Command=~/example/i);
$SNP=0.001 if(!defined $SNP);
$SV="0.1:3000" if(!defined $SV);
$Coverage||=3;
if(defined $Circle){$Circle=1;}
else{$Circle=0}
$Error="0.0005:0.34:0.33:0.33" if(!defined $Error);
die "Wrong error parameter \n" if(($Error=~tr/:/:/)!=3);
die "Wrong SNP parameter \n$USAGE" if((!($SNP=~/\d+\.?\d*/) && !(-e $SNP)) || ($SNP=~/\d+\.?\d*/ && $SNP>1));
die "Wrong SV parameter \n$USAGE" if( !(-e $SV) && (!$SV=~/\d*\.?\d+\:\d+/));
$Qtype||="!";
$Qmean||=37;
$Qsd||=1;
my ($l1,$l2);
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

if($Palo)
{
	die color("red"),"\nOooopssssssssss\nAs to my point, Roche 454 PE sequencing technoly is not suitable for palogenome sequencing!\nPalogenome are sheared into short pieces after long time digestion.\nHow do you cut the gel while library preparation???\nSorry, I don't wanna simulate this process~\n\n",color("reset") if($Command=~/roche/i && $PE);
	$FragLim||="20+";
	$LostSingle||=0.5;
	$DoublePercent||=0.9;
	$Mismatch||="3_20";
	$Mutation||=0,01;
	$Library||=1;
}
$Dir||="./";
die  "$Dir do not exist!\nplease make directory $Dir first\n" if(!(-d $Dir));
my $flowcell=1;
my $titlenumber=1;
my $n=0;

my (%DigestRegion,@rna_poisson);
my @cDNAFragment=qw();
#my ($gene_num,$mRNA_num)=(0,0);
if($CDNA)
{
	$CovcDNA||=10;
	&ReadGff($GFF,\%DigestRegion);
	my $NumcDNA=$CovcDNA*scalar(keys %DigestRegion);
	die "NO TRANSCRIPT INFORMATION IN $GFF!\n" if($NumcDNA==0);
	@rna_poisson=random_poisson($NumcDNA,$CovcDNA);
	$Circle=0;
}

my $info="Simulation Sequencing Parameters
===========================================================
reference data $Reference\n";
if($Circle){$info.="circle reference\n"}
else{$info.="liner reference\n"}
if($SNP eq 0){$info.="no snp sites\n"}
elsif($SNP=~/\d+\.?\d*/){$info.="snp rate $SNP\n"}
else{$info.="snp rate file $SNP\n"}
if($SV eq 0){$info.="no sv\n"}
elsif($SV=~/\:/ ){$info.="sv rate and average length $SV\n"}
else{$info.="sv file $SV\n"}
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
if($Error=~/\:/){$info.="error rate and possible of each error type is $Error\n"}
else{$info.="error file is $Error\n";}
if($Palo){$info.="simulation palogenome\nmismatch range is $Mismatch\n$Library stain library preparation method\nmutation rate is $Mutation\n"}
if($DoublePercent){$info.="double strain both kept rate is $DoublePercent\n"}
elsif($LostSingle){$info.="single strain lost rate is $LostSingle\n\n";}
if($CDNA){$info.="RNA-seq\ninput gff file is $GFF\nmean coverage is $CovcDNA\n\n";}
$info.="output directory is $Dir\noutput sequencing fiel prefix is $Output
============================================================\n";
print color("yellow"),$info,color("reset") if($Reference && $Command && !$Help);

#============================CHECK PARAMETER END===============================

#============================HANDLE INFORMATION================================

open OUT1,">$Dir$Output.fasta" if(!($Command=~/illumina/i));
open OUT2,">$Dir$Output.qual" if(!($Command=~/illumina/i));
open OUT1,">$Dir$Output.fastq" if($Command=~/illumina/i && !$PE);
open OUT1,">$Dir$Output\-1.fastq" if($Command=~/illumina/i && $PE);
open OUT2,">$Dir$Output\-2.fastq" if($Command=~/illumina/i && $PE);

open NAME,">$Dir"."NameRecordByPsim.txt";
my (@SV,@SNP,@MUTATION,@NAME,@ERROR);

#==============================MAIN PROGRAM===============================
my @Qualitys;
my $Sequence="";
my $SeqName="";
open REFERENCE,"<$Reference";
while(<REFERENCE>)
{
	if(/^>/)
	{
		&Main if($Sequence ne "" && length($Sequence)>0);
		chomp;
		$SeqName=$_;
		$SeqName=~s/^>//;
		$Sequence="";
	}
	else
	{
		chomp;
		$Sequence.=$_;
	}
}
close REFERENCE;
&Main() if (length($Sequence)>0);
$Sequence="";

#=============================MAIN SUBROUTINE===============================
sub Main
{
#RNA-seq sequencing
	if($CDNA)
	{
		my $n=0;
		my @name=split /\s+/,$SeqName;
		$name[0]=~s/^>//;
		my $CHR=$name[0];
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
						$transcript.=substr($Sequence,$start,($end-$start+1));
					}
				}
				elsif($Strand eq "-")
				{
					foreach(@se)
					{
						my ($start,$end)=split /\,/;
						$transcript=substr($Sequence,$start,($end-$start+1)).$transcript;
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
#STEP ONE:	SNP
		if($SNP ne 0)
		{
			my $snptype=0;
			$snptype=1 if(!($SNP=~/\d\.?\d*/));
			&SNP(\$Sequence,$snptype,$SeqName,$SNP,\@SNP);
		}

#STEP TWO:	SV
		if($SV ne 0)
		{
			my $SVType=1;
			$SVType=0 if($SV=~/\d*\.?\d+\:\d+/);
			&StructuralVariation(\$Sequence,$SeqName,$SVType,$SV,\@SV)
		}

#STEP THREE:	GENERATE LIBRARY FRAGMENT AND SEQUENCING
#====GENERATE LIBRARY
		&Library(\$Sequence,$Coverage);
	}
}
sub Library
{
	my ($Sequence,$Coverage)=@_;
    my $ltype=0;
	my @Fragment;
	$ltype=$Insert.":".length($Linker) if($Command=~/roche/i && $PE);
	$$Sequence.=substr($$Sequence,0,int($FragmentMean*2));
	my @StartLeng=&PieceGenerate(length($$Sequence),$ltype,$Circle,$Coverage,$FragmentMean,$FragmentSD,$FragLim);
#====QUALITY
	if($Command=~/illumina/ || $Command=~/solid/)
	{
		my $ReadsNum=scalar(@StartLeng);
		$ReadsNum=2*$ReadsNum if($PE);
		&Quality($Command,\@Qualitys,$ReadLeng,$Qmean,$Qsd,$ReadsNum,$Qtype);
	}
	else
	{
		my @ReadsLeng;
		map {my @temp=split /\_/;push @ReadsLeng,$temp[-1]} @StartLeng;
		&Quality($Command,\@Qualitys,\@ReadsLeng,$Qmean,$Qsd);
	}
	if($Palo)
	{
		my $lost=$LostSingle;
		$lost=$DoublePercent if($Library==2);
		my @PaloName;
		&MismatchDoubleStrain($Library,$Sequence,\@StartLeng,$Mismatch,$SeqName,$lost,$Mutation,\@Fragment,\@PaloName,\@MUTATION);
#====SEQUENCING
		for my $i(0..$#Fragment)
		{
			&SequencingIlluminaSE($Fragment[$i],$PaloName[$i],\$flowcell,\$titlenumber,\$n) if($Command=~/illumina/i && !$PE);
			&SequencingIlluminaPE($Fragment[$i],$PaloName[$i],\$flowcell,\$titlenumber,\$n) if($Command=~/illumina/i && $PE);
			&SequencingRoche($Fragment[$i],$PaloName[$i],\$titlenumber) if($Command=~/roche/i);
			&SequencingSOLiDSE($Fragment[$i],$PaloName[$i],\$flowcell,\$titlenumber,\$n) if($Command=~/solid/i && !$PE);
			&SequencingSOLiDPE($Fragment[$i],$PaloName[$i],\$flowcell,\$titlenumber,\$n) if($Command=~/solid/i && $PE);
		}
	}
#====GENERATE LIBRARY
	else
	{
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
					my $orileng=$temp_roche[3]-$temp_roche[2]+1;
					my $linkeradd=$temp_roche[3]-$orileng;
					$sequence=substr($$Sequence,$temp_roche[1],$orileng).substr($Linker,0,$linkeradd);
					my $readinfo="$SeqName-3end $i start=$temp_roche[1] end=$temp_roche[2] length=$temp_roche[3]";
#====SEQUENCING
					&SequencingRoche($sequence,$readinfo,\$titlenumber) if($Command=~/roche/i);
				}
				elsif($temp_roche[0]==0)
				{
					my $ori1=$temp_roche[4]-$temp_roche[2]+1;
					my $ori2=$temp_roche[4]-length($Linker);
					$sequence=substr($$Sequence,$temp_roche[1],$ori1).$Linker.substr($$Sequence,$temp_roche[3],$ori2);
					my $readinfo="$SeqName-bothPE $i start=$temp_roche[3] end=$temp_roche[2] length=$temp_roche[4]";
#====SEQUENCING
					&SequencingRoche($sequence,$readinfo,\$titlenumber) if($Command=~/roche/i);
				}
				else
				{
					my $orileng=$temp_roche[2]-$temp_roche[1]+1;
					my $linkeradd="-".($temp_roche[3]-$orileng);
					$sequence=substr($Linker,$linkeradd).substr($$Sequence,$temp_roche[1],$orileng);
					my $readinfo="$SeqName-5end $i start=$temp_roche[1] end=$temp_roche[2] length=$temp_roche[3]";
#====SEQUENCING
					&SequencingRoche($sequence,$readinfo,\$titlenumber) if($Command=~/roche/i);
				}
			}
#====SEQUENCING
			else
			{
				my ($start,$length)=split /\_/,$StartLeng[$i];
#			print "TEST\tstart=$start\tlength=$length\n";
				$sequence=substr($$Sequence,$start,$length);
				my $readinfo="$SeqName $i start=$start end=".($start+$length-1)." length=$length";
				&SequencingIlluminaSE($sequence,$readinfo,\$flowcell,\$titlenumber,\$n) if($Command=~/illumina/i && !$PE);
				&SequencingIlluminaPE($sequence,$readinfo,\$flowcell,\$titlenumber,\$n) if($Command=~/illumina/i && $PE);
				&SequencingRoche($sequence,$readinfo,\$titlenumber) if($Command=~/roche/i);
				&SequencingSOLiDSE($sequence,$readinfo,\$flowcell,\$titlenumber,\$n) if($Command=~/solid/i && !$PE);
				&SequencingSOLiDPE($sequence,$readinfo,\$flowcell,\$titlenumber,\$n) if($Command=~/solid/i && $PE);
			}
		}
	}
}
sub SequencingIlluminaSE
{
	my ($seq,$readinfo,$f,$t,$n)=@_;
	$seq.=$adapter2;
	my $x=int(rand(20000))+1;
	my $y=int(rand(21500))+1;
	my $prefix="\@RachelWu\:1\:BENM\:$$f\:$$t\:$x\:$y 0\:N\:0\:ATCACG";
	print NAME "$prefix\t$readinfo\n";
	my $read=substr($seq,0,$ReadLeng);
	$read=&SequencingError($read,$Error,$prefix,\@ERROR);

	print OUT1 "$prefix\n$read\n\+\n$Qualitys[$count]\n";

	$$n+=1;
	if($$n>400)
	{
		$$t+=1;
		$$n=1;
	}
	$$f+=1 if($$t>120);
	$count+=1;
}

sub SequencingIlluminaPE
{
	my ($reads,$readinfo,$f,$t,$n)=@_;
	my $readseq1=$reads.$adapter2;
	my $readseq2=&ReverseComplement($adapter1.$reads);

	my $x=int(rand(20000))+1;
	my $y=int(rand(21500))+1;
	my $prefix="\@RachelWu\:$$f\:$$t\:$x\:$y\#0";
	print NAME "$prefix\t$readinfo\n";

	my $pe1=substr($readseq1,0,$ReadLeng);
	my $pe2=substr($readseq2,0,$ReadLeng);

	$pe1=&SequencingError($pe1,$Error,$prefix." 1\:N\:0\:ATCACG",\@ERROR);
	$pe2=&SequencingError($pe2,$Error,$prefix." 2\:N\:0\:ATCACG",\@ERROR);

	print OUT1 "$prefix"." 1\:N\:0\:ATCACG\n$pe1\n";
	print OUT1 "+\n"."$Qualitys[$count]\n";
	print OUT2 "$prefix"." 2\:N\:0\:ATCACG\n$pe2\n";
	print OUT2 "+\n"."$Qualitys[$count+1]\n";

	$$n+=1;
	if($$n>400)
	{
		$$t+=1;
		$$n=1;
	}
	$$f+=1 if($$t>120);
	$count+=2;
}

sub SequencingRoche
{
	my ($reads,$readinfo,$n)=@_;
	my $prefix="\@Psim_Roche454\.$$n";
	print NAME "$prefix\t$readinfo\n";

	&SequencingError($reads,$Error,$prefix,\@ERROR);
	print OUT1 "$prefix\n$reads\n";
	print OUT2 "$prefix\n$Qualitys[$n-1]\n";
	$$n+=1;
}

sub SequencingSOLiDSE 
{
	my ($seq,$readinfo,$f,$t,$n)=@_;
	my $prefix="\>RachelWu\_$$f\_$$t\_$$n\_F3";
	print NAME "$prefix\t$readinfo\n";
	my $read=substr($seq,0,$ReadLeng);
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
	$count+=1;
}
sub SequencingSOLiDPE
{
	my ($reads,$readinfo,$f,$t,$n)=@_;
	my $prefix="\>RachelWu\_$$f\_$$t\_$$n";
	print NAME "$prefix\t$readinfo\n";
	my $readseq1=$reads;
	my $readseq2=&ReverseComplement($reads);

	my $pe1=substr($readseq1,0,$ReadLeng);
	my $pe2=substr($readseq2,0,$ReadLeng);

	$pe1=&SequencingError($pe1,$Error,$prefix."\_F3",\@ERROR);
	$pe2=&SequencingError($pe2,$Error,$prefix."\_R3",\@ERROR);
	&Base2Col(\$pe1,$Header);
	&Base2Col(\$pe2,$Header);

	print OUT1 "$prefix"."\_F3\n$pe1\n";
	print OUT2 "$prefix"."\_F3\n$Qualitys[$count]\n";
	print OUT1 "$prefix"."\_R3\n$pe2\n";
	print OUT2 "$prefix"."\_R3\n$Qualitys[$count+1]\n";
	$$n+=1;
	if($$n>400)
	{
		$$t+=1;
		$$n=1;
	}
	$count+=2;
}

if($SNP ne 0)
{
	open SNPOUT,">$Dir"."SNPReportByPsim.txt" if($SNP ne 0);
	map {print SNPOUT "$_\n"} @SNP;
	close SNPOUT;
}
if($SV ne 0)
{
	open SVOUT,">$Dir"."SVReportByPsim.txt" if($SV ne 0);
	map {print SVOUT "$_\n"} @SV;
	close SVOUT;
}
if($Error ne 0)
{
	open ERROR,">$Dir"."SequencingErrorByPsim.txt" if($Error ne 0);
	map {print ERROR "$_\n"} @ERROR;
	close ERROR;
}
if($Mutation)
{
	open MUTATION,">$Dir"."MutationByPsim.txt" if($Mutation);
	map {print "$_\n"} @MUTATION;
	close MUTATION;
}
close NAME;
close OUT1;
close OUT2 || print "Simulation Done\n";
