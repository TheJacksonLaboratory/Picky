#!/usr/bin/env perl
#####
#
# Jackson Laboratory Non-Commercial License
# See the LICENSE file (LICENSE.txt) for license rights and limitations
# https://github.com/TheJacksonLaboratory/Picky/blob/master/LICENSE.md
#
# Picky - Structural Variants Pipeline for long read
#
# Repository: https://github.com/TheJacksonLaboratory/Picky
# Documentation: https://github.com/TheJacksonLaboratory/Picky/wiki
#
# Created Aug 16, 2016
# Copyright (c) 2016-2017  Chee-Hong WONG
#                          Genome Technologies
#                          The Jackson Laboratory
#
#####

use strict;
use Getopt::Long;
use FindBin;
use lib $FindBin::Bin;
use hashFastq;
use selectAlignmentMT; # threading version
use callSV;
use xlsToVCF;
use PBSCluster;
use samToAlign;

my $G_USAGE = "
$0 <command> -h

<command> [hashFq, selectRep, callSV]
hashFq    : hash read uuids to friendly ids
lastParam : Last parameters for alignment
selectRep : select representative alignments for read
callSV    : call structural variants
xls2vcf   : convert Picky sv xls file to vcf
sam2align : convert sam to align format
preparepbs: chunk last fastq file and write pbs script for cluster submission
script    : write a bash-script for single fastq processing


      License = JACKSON LABORATORY NON-COMMERCIAL LICENSE AGREEMENT
                https://github.com/TheJacksonLaboratory/Picky/blob/master/LICENSE.md
Documentation = https://github.com/TheJacksonLaboratory/Picky/wiki
   Repository = https://github.com/TheJacksonLaboratory/Picky/
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
} elsif ('sam2align' eq $command) {
	runSamToAlign();
} elsif ('preparepbs' eq $command) {
	runPrepareForPBSCluster();
} elsif ('script' eq $command) {
	runWriteScript();
} else {
	print $G_USAGE;
}

exit 0;
