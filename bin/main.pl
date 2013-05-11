#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Term::ANSIColor;
use FindBin qw($Bin);
use lib "$Bin/lib/";
use 5.010;






die "
USAGE:   ****************
Author:  RachelWu\(rachelwu123\@gmail.com\)\,BENM\(binxiaofeng\@gmail.com\)
Version: *******

Run:     perl $0 [command] [options]

command: Illumina
         Roche454
         SoLid
         Palo
\n" if(!$ARGV[0]);

my $command=$ARGV[0];
given($command)
{

