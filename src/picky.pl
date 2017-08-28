#!/usr/bin/env perl
$|++;

#####
#
# Picky - Structural Variants Pipeline for (ONT) long read
#
# Copyright (c) 2016-2017  Chee-Hong WONG  The Jackson Laboratory
#
#####

#####
# picky.pl
#
# hashFq Examples:
#   perl picky.pl hashFq -pfile 2016-08-24-R9-WTD-R009.pass.fastq -ffile 2016-08-24-R9-WTD-R009.fail.fastq -oprefix WTD09
#   perl picky.pl hashFq -pfile 2016-08-24-R9-WTD-R009.pass.fastq -oprefix WTD09P
#   perl picky.pl hashFq -ffile 2016-08-24-R9-WTD-R009.fail.fastq -oprefix WTD09F
#
# selectRep Examples:
#   perl picky.pl selectRep --in WTD09.2D.subset.v755.hg19.maf
####   last-755/src/lastal -r1 -q1 -a0 -b2 -v -P 28 -Q 1 hg19.lastdb ${RUNTYPE}.fastq 2>${RUNTYPE}.v755.hg19.log | perl picky.pl selectRep --in -
#
# callSV Examples:
#   perl picky.pl callSV --in WTD09.2D.subset.v755.hg19.multiple.align --fastq WTD09.2D.subset.fastq --genome hg19.fa --exclude=chrY --exclude=chrM --readseq --correctminus --removehomopolymerdeletion 2>&1 | tee WTD09.2D.subset.homopolymer.log
#   perl picky.pl callSV --in WTD09.2D.subset.v755.hg19.multiple.align --fastq WTD09.2D.subset.fastq --genome hg19.fa --exclude=chrY --exclude=chrM --readseq --correctminus --removehomopolymerdeletion --sam 2>&1 | tee WTD09.2D.subset.homopolymer.log
#
#####

#####
# TODO:
# 1. command to generate the lastal command for processing
#

use strict;
use Getopt::Long;
use hashFastq;
#use selectAlignment; # non-threading version
use selectAlignmentMT; # threading version
use callSV;

my $G_USAGE = "
$0 <command> -h

<command> [hashFq, selectRep, callSV]
hashFq : hash read uuids to friendly ids
selectRep : select representative alignments for read
callSV    : call structural variants
";

my $command = undef;
if (!defined $ARGV[0] || substr($ARGV[0],0,1) eq '-' || substr($ARGV[0],0,2) eq '--') {
	die("Please specify the command.\n",$G_USAGE);
}
$command = shift @ARGV; $command = lc $command;

if ('hashfq' eq $command) {
	runHashFastq();
} elsif ('selectrep' eq $command) {
	runSelectRepresentativeAlignments();
} elsif ('callsv' eq $command) {
	runCallStructuralVariants();
} else {
	print $G_USAGE;
}

exit 0;

