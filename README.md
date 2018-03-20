## PSim--a pilot study sequencing simulator.
**Author**: Rachel Wu(rachelwu123@gmail.com ), BENM Feng(BinxiaoFeng@gmail.com )  
**Version**: 1.0.0  
**Update**: 2013.06.16  
**Depends**: perl, Math-Random    
### Introduction  
PSim is a pilot study sequencing simulator. It can simulate NGS sequencing data such as illuminaTM, 454 FLXTM , SOLiDTM , Ion TorrentTM or PacBioTM etc. It satisfied variety requires, using this tool to provide irregular or regular reads length, paired-end reads with hieratical insert span, single strand specified library, BAC liked loop sequences, pooling strategy, variation by randomly generated or custom introduced including single nucleotide polymorphism (SNP), small insertion and deletion (InDel) and structural variation (deletion, insertion, inversion, tandem repeat and translocation), DNA damage such as methylation or deamination, particular sequencing error models, library construction, pollution mixture, metagenome with multi-index or barcode, random passion/bimodal distribution, it also varied reference/ancestral sequence based on Bayesian probability models. Beyond these, Monte Carlo method is presented for pilot study by sequencing simulator PSim, with wide utility in generating DNA/mRNA sequencing with different sequencing data format.   

The real sequencing library preparation process is simulated in PSim. The process include template DNA variation,\ *DNA damage, \*DNA fragment pre-amplification, \*single-chain library or double-strain library, randomly fragmentation, add adapter or PE linker, generation of quality value and so on.   
\*especially for damage DNA sequencing steps   
### UPDATE
Version 1.2:	function of PCR duplication added   
Version 1.3:	the reference sequence could have "N". And reads include "N" will be ignored    
Version 1.4:	length of damage consist of two Normal-like distribution, and have some hidden parameters   
Version 1.5:	add pod in each pm file   
### Obtaining and Installing
PSim is available from https://sourceforge.net/projects/pilotsim/. PSim depends on perl and Math-Random module.   
####Installation of perl:  
1)  Perl always has been pre-installed on Linux system. There is no need to reinstall;   
2)  Windows system user can download activeperl from http://www.activestate.com/activeperl.     
#### Installation of Math-Ranodm:   
1)  Automatically install with CPAN (need to be root)    

	#Linux user
	perl –MCPAN –e shell
	cpan>install Math::Random
	cpan>q
	#Windows user
	ppm
	ppm>install Math::Random

2)  Manual installation   
Download Math-Random source code from http://search.cpan.org/~grommel/Math-Random-0.70/Random.pm and then type in a terminal:     

	tar zxvf Math-Random-0.70.tar.gz
	cd Math-Random-0.70
	perl Makefile.PL
	make
	make test
	make install
	#manually set the installation path if you don’t have root privilege
	tar zxvf Math-Random-0.70.tar.gz
	cd Math-Random-0.70
	perl Makefile.PL PREFIX=YOURDIR
	make
	make test
	make install

### Functions
#### 1) Randomly generate reference sequence   
Generate specified length or random length of Poisson distribution of sequence(s) which could be used ad reference sequence(s).      

	perl GenerateSeq.pl [-n] [-l] [-pre] [-auto] > OUTFILE
	  -n NUM        	     number of sequence (default=1)
	  -l NUM or NUM,NUM...	length of each sequence, separate different length by "," 
	  -pre CHAR             prefix of your sequence(s)
	  -auto                 randomly generate n sequences for Poisson distribution with average length of “–l” 

**Example**:  

	#generate five sequences of 3000bp, 4000bp, 6000bp, 5000bp and 8000bp. Prefix of sequences is “test-1”. And the sequences file saved as “test-1.fa”
	perl GenerateSeq.pl –n 5 –l 3000,4000,6000,5000,8000 –pre test-1 > test-1.fa
	#generate five sequences of average length of 5000bp. Pre ix of sequences is “test-2”. And the sequences file saved as “test-2.fa”
	perl GenerateSeq.pl –n 5 –l 5000 –pre test-2 –auto > test-2.fa

#### 2) Structure Variation   
Structural variation is the variation in structure of an organism's chromosome. It consists of many kinds of variation, such as deletion, insertion, tandem repeat, reversion and translocation. Typically a structure variation affects a sequence length about 1Kb to 3Mb. While, indel of short segment is also included in PSim. To specify variation location(s) and type(s), variation information file like /example/sv\_example.txt should be offered while running PSim.    

The variation information file format as: the first column with reference name, the second column with variation type (include deletion、insertion、tandem\_repeat、reversion and translocation), the third column with location, length or other information. Such as  

	deletion	site1, length1; site2, length2……
	insertion	insert_site1,sequence1; insert_site2,sequence2……
	inversion	site1, length1; site2, length2……
	tandem_repeat	site1, length1, repeat_times1 ;site2, length2, repeat_times2……
	translocation	site1, length1, insert_site1; site2, length2, insert_site2……
**Parameter**:   

	--sv <NUM|FILE>	structure variation rate and average length (default= 0.1:3000) or specific sv type and site file(format refer to ../example/sv_example.txt)
**Example**:   

	#variation as SV information file (like DIR/PSim/example/sv_example.txt)
	--sv DIR/PSim/input/sv_example.txt
	#randomly generate SV of average length of 3000bp and coverage of 10%
	--sv 0.1:3000
	#no SV
	--sv 0

#### 3) SNP   
To specify SNP site(s) and base(s), mutation information file like /example/snp\_example.txt should be offered while running PSim. The mutation information file format as: the first column with reference name, the second column with mutation site, the third column with possible mutation base(s) and the fourth column with mutation rate(S). Separate different sites and rates in column three and four by comma.     
**Parameter**:    

	--snp <NUM|FILE>	random snp rate (default=0.001) or specific snp site and rate file (format refer to ../example/snp_examp le.txt)
**Example**:  

	#mutation as SNP information file DIR/PSim/example/snp_example.txt
	--snp DIR/PSim/example/snp_example.txt
	#randomly generate SNP of 10% coverage
	--snp 0.001
	#no SNP
	--snp 0
#### 4) Library preparation    
Library preparation process includes random fragmentation of library sample and adapter or PE linker added. For fixed read length sequencing methods, most or all the lengths of fragments should be set above read length. And for non-fixed read length sequencing methods, the length of fragment is the read length.    

In this section, user need to specify the average length of fragment (-fragmean), stander deviation (-fragsd) and length limit (-fraglim, optional). Adapter sequences are added by Illumina like sequencing methods and linker sequence is added by Roche 454 PE like sequencing methods.     
**Parameter**: 

	--cov <NUM>	sequencing coverage of reference sequence(s) (default=3)
	--fragmean <INT>	average length of library fragment (default=200)
	--fragsd <INT>	standard deviation of library fragment (default=10)
	--fraglim <INT>	limit length of fragment library ("20+" means must above 20nt, and "240-" means must shorter than 240nt,if(-damage) this parameter default=20+)
	--adapter <CHAR>	adapter sequence (default sequence refer to ../example/adapterIllumina_example.txt, split two adapters by ":")
	While adopting Roche 454 PE sequencing method:
	--insert <INT>    average insert size and sd(default=8000:30 if --pe)
	--linker <CHAR>   PE insert sequence (default=\"ATAACTTCGTATAATGTATGCT ATACGAAGTTAT\")
**Example**:  

	#simulating sequencing data of 150bp fragments average length, 6 stander deviation, 5 folds coverage and use the default adapter sequences
	--cov 5 --fragmean 150 --fragsd 6 
	#simulating Roche 454 PE sequencing data of 8000bp insert size, 35 stander deviation of insert length, 4 folds coverage, 650bp fragments average length, 50 stander deviation of fragment length and use the default linker sequence
	--pe --cov 4 --insert 8000:35 --fragmean 650 --fragsd 50
#### 5) RNA-seq   
Simulate RNA-seq based on given GFF (General Feature Format) file. Refer to Sanger website (http://www.sanger.ac.uk/resources/software/gff/spec.html) for GFF file format detail.    
**Parameter**:   

	--cdna				RNA-SEQ sequencing
	--gff <FILE>		input gene gff file
	--covcdna <INT>	mean coverage of cDNA sequences(default=10 if --cdna)
**Example**:  

	#simulating RNA-seq result with 20-fold coverage, gff file at DIR/PSim/example/gff_exampl e.txt
	--cdna --gff DIR/PSim/example/gff_example.txt –covcdna 20
##### 6) Damage DNA   
##### 7) Sequencing    
For fixed read length sequencing methods (such as Illumina and SOLiD), user should set read length with the parameter of “-read”.     
### OUTPUT FILES
	NewReferenceSequence.fa--new reference sequence with snp or(and) sv
	output.fa, output.qual or output.fastq--sequencing reads and quality scores
	NameRecordByPSim.txt--sequencing name whit reads info
	SVReportByPSim.txt--SV report
	SNPReportByPSim.txt--SNP report
	SequencingErrorByPSim.txt--sequencing error report
	MutationByPSim.txt--palo genome damage report

