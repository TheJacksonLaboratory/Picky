#!/usr/bin/env perl
#####
#
# Jackson Laboratory Non-Commercial License
# See the LICENSE file (LICENSE.txt) for license rights and limitations
#
# Picky - Structural Variants Pipeline for long read
#
# Created Aug 16, 2016
# Copyright (c) 2016-2017  Chee-Hong WONG
#                          Genome Technologies
#                          The Jackson Laboratory
#
#####

#####
# picky.pl
#
# hashFq Examples:
#   picky.pl hashFq -pfile 2016-08-24-R9-WTD-R009.pass.fastq -ffile 2016-08-24-R9-WTD-R009.fail.fastq -oprefix WTD09
#   picky.pl hashFq -pfile 2016-08-24-R9-WTD-R009.pass.fastq -oprefix WTD09P
#   picky.pl hashFq -ffile 2016-08-24-R9-WTD-R009.fail.fastq -oprefix WTD09F
#
# selectRep Examples:
#   last-755/src/lastal -r1 -q1 -a0 -b2 -P4 -Q1 hg19.lastdb ${RUNTYPE}.fastq 1>${RUNTYPE}.v755.hg19.maf 2>${RUNTYPE}.v755.hg19.log ; \
#   cat ${RUNTYPE}.v755.hg19.maf | ./picky.pl selectRep --thread 4 --preload 6 1>${RUNTYPE}.v755.hg19.selectRep.align 2>${RUNTYPE}.v755.hg19.selectRep.log
#   .. or ..
#   last-755/src/lastal -r1 -q1 -a0 -b2 -P4 -Q1 hg19.lastdb ${RUNTYPE}.fastq 2>${RUNTYPE}.v755.hg19.log \
#   | ./picky.pl selectRep --thread 4 --preload 6 1>${RUNTYPE}.v755.hg19.selectRep.align 2>${RUNTYPE}.v755.hg19.selectRep.log
#
# callSV Examples:
#   cat ${RUNTYPE}.v755.hg19.selectRep.align | ./picky.pl callSV --oprefix ${RUNTYPE}.v755.hg19.selectRep --fastq ${RUNTYPE}.fastq \
#     --genome hg19.fa --removehomopolymerdeletion --exclude=chrY --exclude=chrM --sam 2>${RUNTYPE}.v755.hg19.callSV.log
#   .. or ..
#   cat ${RUNTYPE}.v755.hg19.selectRep.align | ./picky.pl callSV --oprefix ${RUNTYPE}.v755.hg19.selectRep --fastq ${RUNTYPE}.fastq \
#     --genome hg19.fa --removehomopolymerdeletion --exclude=chrY --exclude=chrM 2>${RUNTYPE}.v755.hg19.callSV.log
#   .. or ..
#   cat ${RUNTYPE}.v755.hg19.selectRep.align | ./picky.pl callSV --oprefix ${RUNTYPE}.v755.hg19.selectRep --fastq ${RUNTYPE}.fastq \
#     --exclude=chrY --exclude=chrM 2>${RUNTYPE}.v755.hg19.callSV.log
#
# lastParam Examples:
#   ./Picky.pl lastParam
#
#####

#####
# TODO:
# 1. command to generate the lastal command for processing
#

use strict;
use Getopt::Long;
use hashFastq;
use selectAlignmentMT; # threading version
use callSV;
use xlsToVCF;
use PBSCluster;

my $G_USAGE = "
$0 <command> -h

<command> [hashFq, selectRep, callSV]
hashFq    : hash read uuids to friendly ids
lastParam : Last parameters for alignment
selectRep : select representative alignments for read
callSV    : call structural variants
xls2vcf   : convert Picky sv xls file to vcf
preparepbs: chunk last fastq file and write pbs script for cluster submission
script    : write a bash-script for single fastq processing
";

my $command = undef;
if (!defined $ARGV[0] || substr($ARGV[0],0,1) eq '-' || substr($ARGV[0],0,2) eq '--') {
	die("Please specify the command.\n",$G_USAGE);
}
$command = shift @ARGV; $command = lc $command;

# auto-flush for both stderr and stdout
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if ('hashfq' eq $command) {
	runHashFastq();
} elsif ('selectrep' eq $command) {
	runSelectRepresentativeAlignments();
} elsif ('callsv' eq $command) {
	runCallStructuralVariants();
} elsif ('lastparam' eq $command) {
	#printf "-r1 -q1 -a0 -b2 -v -Q1"; #"-P<threads>"
	printf "-C2 -K2 -r1 -q3 -a2 -b1 -v -Q1"; #"-P<threads>"
} elsif ('xls2vcf' eq $command) {
	runXLStoVCF();
} elsif ('preparepbs' eq $command) {
	runPrepareForPBSCluster();
} elsif ('script' eq $command) {
	runWriteScript();
} else {
	print $G_USAGE;
}

exit 0;
