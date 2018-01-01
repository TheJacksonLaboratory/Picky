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

package PBSCluster;

BEGIN {
	$VERSION = '0.1';
}

use strict;
use warnings;

use Data::Dumper;
use Storable qw(dclone);
use Getopt::Long;
use base 'Exporter';
our @EXPORT = qw(runPrepareForPBSCluster runWriteScript);

my $G_PBS_CHUNKSIZE = 1000;
my $G_PBS_TEMPLATE_FILE = 'template.pbs';

sub runPrepareForPBSCluster {
	my $G_USAGE = "
$0 preparePBS --fastq <fastq_file> [--chunksize <numberOfReadsPerChunk>] [--template <template_file>]

  --fastq STR      fastq file
  --chunksize INT  number of fastq record per chunk file [default: $G_PBS_CHUNKSIZE]
  --template STR   template file for PBS script [default: $G_PBS_TEMPLATE_FILE]
  --init STR       write a copy of the template to specific file
";

	my $fastqFile = '';
	my $chunksize = $G_PBS_CHUNKSIZE;
	my $template = $G_PBS_TEMPLATE_FILE;
	my $inittemplate = '';
	my $help = 0;
	
	GetOptions (
	"fastq=s"     => \$fastqFile,
	"chunksize=i" => \$chunksize,
	"template=s"  => \$template,
	"init=s"      => \$inittemplate,
	"help!"       => \$help)
	or die("Error in command line arguments\n$G_USAGE");

	die "$G_USAGE" if ($help);

	if ('' ne $inittemplate) {
		_initTemplate($inittemplate);
	} else {
		die "$G_USAGE" if ($chunksize<1);
		die "Template file '$template' does not exists!\n$G_USAGE" if (! -f $template);
		die "Fastq file '$fastqFile' does not exists!\n$G_USAGE" if (! -f $fastqFile);
		_writeFastqChunkAndScript($fastqFile, $chunksize, $template);
	}
}

sub runWriteScript {
	my $fastqFile = '';
	my $help = 0;
	my $original2DScore = 0;
	my $scoreParams = '-C2 -K2 -r1 -q3 -a2 -b1 -v -v';

	my $threads = 4;
	my $preload = 6;
	
my $G_USAGE = "
$0 script --fastq <fastq_file> [--thread <numberOfThreads>] [--preload <preload_fold_of_reads_alignments>]

  --fastq STR    fastq file
  --thread INT   number of threads to use [default: $threads]
  --preload INT  fold of preloading for read alignments [default: $preload]
";

	GetOptions (
	"fastq=s"     => \$fastqFile,
	"thread=i"    => \$threads,
	"preload=i"   => \$preload,
	"2D!"         => \$original2DScore,
	"help!"       => \$help)
	or die("Error in command line arguments\n$G_USAGE");

	die "$G_USAGE" if ($help);
	die "Fastq file '$fastqFile' does not exists!\n$G_USAGE" if (! -f $fastqFile);
	die "Expect Fastq file to end with .fastq but see $fastqFile\n" if ($fastqFile !~ /\.fastq$/);
	my $run = $fastqFile; $run =~ s/\.fastq$//;

	$scoreParams = '-r1 -q1 -a0 -b2 -v' if (0!=$original2DScore);

	print '# general installation', "\n";
	print 'export LASTAL=last-755/src/lastal', "\n";
	print 'export PICKY=./picky.pl', "\n";
	print "\n";
	print '# general configuration', "\n";
	print 'export LASTALDB=hg19.lastdb', "\n";
	print 'export LASTALDBFASTA=hg19.fa', "\n";
	print 'export NTHREADS=',$threads, "\n";
	print "\n";
	print '# FASTQ file', "\n";
	print 'export RUN=',$run, "\n";
	print "\n";
	print '# reads alignments', "\n";
	print 'time (${LASTAL} '.$scoreParams.' -P${NTHREADS} -Q1 ${LASTALDB} ${RUN}.fastq 2>${RUN}.lastal.log | ${PICKY} selectRep --thread ${NTHREADS} --preload '.$preload.' 1>${RUN}.align 2>${RUN}.selectRep.log)', "\n";
	print "\n";
	print '# call SVs', "\n";
	print 'time (cat ${RUN}.align | ${PICKY} callSV --oprefix ${RUN} --fastq ${RUN}.fastq --genome ${LASTALDBFASTA} --exclude=chrM --sam)', "\n";
	print "\n";
	print '# generate VCF format', "\n";
	print '${PICKY} xls2vcf --xls ${RUN}.profile.DEL.xls --xls ${RUN}.profile.INS.xls --xls ${RUN}.profile.INDEL.xls --xls ${RUN}.profile.INV.xls --xls ${RUN}.profile.TTLC.xls --xls ${RUN}.profile.TDSR.xls --xls ${RUN}.profile.TDC.xls > ${RUN}.allsv.vcf', "\n";
}

sub _initTemplate {
	my ($file) = @_;

	my $threads = 16;
	my $preload = 6;
	my $memory = '18GB';
	my $walltime = '01:00:00';
	open OUTFILE, ">$file" || die "Fail to write $file\n$!\n";
	print OUTFILE '#!/bin/bash', "\n";
	print OUTFILE '#PBS -l nodes=1:ppn='.$threads, "\n";
	print OUTFILE '#PBS -l walltime='.$walltime, "\n";
	print OUTFILE '#PBS -l mem='.$memory, "\n";
	print OUTFILE "\n";
	print OUTFILE 'cd "$PBS_O_WORKDIR"', "\n";
	print OUTFILE 'export LASTAL=last-755/src/lastal', "\n";
	print OUTFILE 'export LASTALDB=hg19.lastdb', "\n";
	print OUTFILE 'export LASTALDBFASTA=hg19.fa', "\n";
	print OUTFILE 'export PICKY=./picky.pl', "\n";
	print OUTFILE 'export RUN=', "\n";
	print OUTFILE 'time (${LASTAL} -v -C2 -K2 -r1 -q3 -a2 -b1 -v -P'.$threads.' -Q1 ${LASTALDB} ${RUN}.fastq 2>${RUN}.lastal.log | ${PICKY} selectRep --thread '.$threads.' --preload '.$preload.' 1>${RUN}.align 2>${RUN}.selectRep.log)', "\n";
	print OUTFILE 'time (cat ${RUN}.align | ${PICKY} callSV --oprefix ${RUN}.sv --fastq ${RUN}.fastq --genome ${LASTALDBFASTA} --exclude=chrM --sam)', "\n";
	close OUTFILE;
}

sub _writeBashScript {
	my ($chunkFastq, $preLinesRef, $postLinesRef) = @_;
	my $prefix = $chunkFastq;
	$prefix=~s/\.fastq$// if ($prefix=~/\.fastq$/);
	my $pbsfile = $prefix.'.pbs';
	open PBSFILE, ">$pbsfile" || die "Fail to open $pbsfile\n$!\n";
	foreach my $line (@{$preLinesRef}) { print PBSFILE $line, "\n"; }
	print PBSFILE 'export RUN='.$prefix, "\n";
	foreach my $line (@{$postLinesRef}) { print PBSFILE $line, "\n"; }
	close PBSFILE;
}

sub _loadTemplate {
	my ($template, $preLinesRef, $postLinesRef) = @_;

	@{$preLinesRef} = ();
	@{$postLinesRef} = ();
	open INFILE, $template || die "Fail to open $template\n$!\n";
	my $destRef = $preLinesRef;
	while (<INFILE>) {
		chomp();
		if (lc($_) eq 'export chunk=') {
			$destRef = $postLinesRef;
		} else {
			push @{$destRef}, $_;
		}
	}
	close INFILE;
}

sub _writeFastqChunkAndScript {
	my ($fastqFile, $chunkSize, $template) = @_;

	my @preLines = ();
	my @postLines = ();
	_loadTemplate($template, \@preLines, \@postLines);

	open INFILE, $fastqFile || die "Fail to open $fastqFile\n$!\n";
	my $prefix = $fastqFile; $prefix=~s/\.fastq$//;
	my $outfile = '';
	my $numRecords = 0;
	my $numChunks = 0;
	while (<INFILE>) {
			my $idLine = $_;
			my $seq = <INFILE>;
			my $delimiter = <INFILE>;
			my $qual = <INFILE>;
			$numRecords++;

			if (1==($numRecords % $chunkSize)) {
					if ($numChunks>0) {
							close OUTFILE;
							_writeBashScript ($outfile, , \@preLines, \@postLines);
							print "Chunk#", ($numChunks+1)," written as ", $outfile, "\n";
					}
					$numChunks++;
					$outfile = sprintf("%s-c%06d.fastq", $prefix, $numChunks);
					open OUTFILE, ">$outfile" || die "Fail to open $outfile\n$!\n";
			}
			print OUTFILE $idLine;
			print OUTFILE $seq;
			print OUTFILE $delimiter;
			print OUTFILE $qual;
	}
	close INFILE;
	if ($numChunks>0) {
			close OUTFILE;
			_writeBashScript ($outfile, , \@preLines, \@postLines);
			print "Chunk#", ($numChunks+1)," written as ", $outfile, "\n";
	}
	print $numRecords, " records into ", $numChunks, " files\n";
	print "To submit the cluster jobs:\n";
	print '    for i in '.$prefix.'-c??????.pbs ; do echo ${i}; qsub ${i}; done', "\n";
}

1;

__END__
