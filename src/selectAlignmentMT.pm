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

package selectAlignmentMT;

BEGIN {
	$VERSION = '0.1';
}

#####
# TODO:
# 1) OO Picker so that we can swap for different algorithms
# 2) decide if  additional data in align help with downstream analysis when piped
#

use strict;
use warnings;

use Data::Dumper;
use Storable qw(dclone freeze thaw);
use Getopt::Long;
use base 'Exporter';
our @EXPORT = qw(runSelectRepresentativeAlignments);

use utilities;
use Thread;
use Thread::Queue; # 3.12;
# WORKAROUND#1: added as substitute does not have a proper multi-threaded Perl using Thread::Queue 3.12 !
use Thread::Semaphore;
# END-WORKAROUND#1
use PickerMT;

my $G_KEEPLINE = 0; # keey maf record lines in memory; NOTE: may consume much memory
my $G_MAX_EG2 = 1.74743e-12; # this is for 10x coverage
my $G_MIN_IDENTITY_VARIANCE = 55; # i.e. %= >=55%

my $G_CACHE_FOLD = 10;
#my $G_NUMBER_OF_WORKERS = 8;
my $G_NUMBER_OF_WORKERS = 1; # for testing purposes

my $G_USAGE = "
$0 selectRep [--thread <numberOfThreads>] [--preload <preloadFold>]

--thread INT   number of threads
--preload INT  Fold of thread count to preload maf records
";

# STDIN for .maf stream
# STDOUT for our .align stream
# STDERR for our logging
sub runSelectRepresentativeAlignments {
	my $numberOfThreads = $G_NUMBER_OF_WORKERS;
	my $preloadFold = $G_CACHE_FOLD;
	my $help = 0;
	
	GetOptions (
	"thread=i" => \$numberOfThreads,
	"preload=i" => \$preloadFold,
	"help!"        => \$help)
	or die("Error in command line arguments\n$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	$numberOfThreads = $G_NUMBER_OF_WORKERS if ($numberOfThreads<1);

	my $maxEG2 = $G_MAX_EG2; #Picker::getMaxEG2();
	my $minIdentityPercent = $G_MIN_IDENTITY_VARIANCE; #Picker::getMinIdentityPercentage();
	
	my $startTime = time;
	
	my %qRequests = ();
	my $qRequests = Thread::Queue->new();
	$qRequests{qRequests} = $qRequests;
	my $requestSlots = Thread::Semaphore->new($numberOfThreads * $preloadFold);
	$qRequests{_hackSlots} = $requestSlots;
	my $qResults = Thread::Queue->new();
	
	# start alignments writer thread
	my $tWriter = Thread->new(\&writePickedAlignments, $qResults, $numberOfThreads);
	
	# start processor worker threads
	for(my $i=0; $i<$numberOfThreads; ++$i) {
		my $tPicker = Thread->new(\&handleAlignmentsPicking, \%qRequests, $qResults);
	}
	
	# start maf record reader in current thread
	my $numRecords = readMAFRecords(\%qRequests, $maxEG2, $minIdentityPercent);
	for(my $i=0; $i<$numberOfThreads; ++$i) {
		$requestSlots->down(); # block until there is a free slot
		$qRequests->enqueue("DONE");
	}
	
	# producer is done, let's signal to all threads
	
	# wait for all threads to finish!
	$_->join for threads->list;
	
	my $totalTime = time - $startTime;
	printf STDERR "%s Total run time = %d secs\n", utilities::getTimeStamp(), $totalTime;
	printf STDERR "%s %d .maf lines processed.\n", utilities::getTimeStamp(), $.;
	printf STDERR "%s %d read(s) processed\n", utilities::getTimeStamp() , $numRecords;
	
	close STDOUT || die "Fail to close STDOUT\n$!\n";
	close STDERR || die "Fail to close STDERR\n$!\n";
}

sub handleAlignmentsPicking {
	my ($qRequests, $qResults) = @_;
	
	while (my $request = $qRequests->{qRequests}->dequeue()) {
		if ("DONE" eq $request) {
			$qResults->enqueue("DONE");
			last;
		}

		my $maxEG2 = $request->{maxEG2};
		my $minIdentityPercent = $request->{minIdentityPercent};
		my $currReadRef = undef;
		my @alignments = ();
		my $linesRef = $request->{lines};
		for (my $i=0; $i<$request->{numAlignments}; ++$i) {
			my %alignment = ();
			my $aline = shift @{$linesRef};
			my $sline = shift @{$linesRef};
			my $qline = shift @{$linesRef};
			my $qualityLine = shift @{$linesRef};

			# let's parse the lines
			utilities::parseMAFRecordGenCigar ($aline, $sline, $qline, $qualityLine, \%alignment, $G_KEEPLINE);
			
			$currReadRef = $alignment{'read'} if (0==$i);
			
			next if ($alignment{'EG2'}>$maxEG2);
			next if ($alignment{'profile'}->{'%='}<$minIdentityPercent);
			push @alignments, \%alignment;
		}
		
		# let's perform alignments picking
		my $currAlignmentsRef = \@alignments;
		my @alignLines = ();
		my @logLines = ();
		PickerMT::processAlignments ($currReadRef, $currAlignmentsRef, $request->{numAlignments}, $request->{ioTime}, \@alignLines, \@logLines);

		# inform writer
		my %results = (requestId=>$request->{requestId}, aligns=>\@alignLines, logs=>\@logLines);
		$qResults->enqueue(\%results);

		$qRequests->{_hackSlots}->up(); # I just free up a slot with job completed
	}
}

sub writePickedAlignments {
	my ($qResults, $numberOfWorkers) = @_;
	
	my $pendingDONE = $numberOfWorkers;
	my %outOfOrderResults = ();
	my $expectedId = 1;
	while (my $result = $qResults->dequeue()) {
		if ("DONE" eq $result) {
			$pendingDONE--;
			last if (0==$pendingDONE);
			next;
		}
		# check if result is in order
		# process in-order result, else cache it
		if ($expectedId == $result->{requestId}) {
			# process as much as possible, including any cached
			$outOfOrderResults{$result->{requestId}} = $result;
			while (exists $outOfOrderResults{$expectedId}) {
				my $currResult = $outOfOrderResults{$expectedId};
				print STDERR join("\n", @{$currResult->{logs}}), "\n" if (scalar(@{$currResult->{logs}})>0);
				print join("\n", @{$currResult->{aligns}}), "\n" if (scalar(@{$currResult->{aligns}})>0);
				delete $outOfOrderResults{$expectedId};
				$expectedId++;
			}
		} else {
			$outOfOrderResults{$result->{requestId}} = $result;
		}
	}
}

sub readMAFRecords {
	my ($qRequests, $maxEG2, $minIdentityPercent) = @_;
	
	my $startTime = time;
	my $currQueryId = '';
	my $currReadRef = undef;
	my $currAlignmentsRef = undef;
	my $endTime = $startTime;
	my $numAlignments = 0;
	my $parseHeader = 1;
	
	my $serialNumber = 0;
	
	while (<STDIN>) {
		if ('#' eq substr($_,0,1)) {
			if (0!=$parseHeader) {
				# LAST version 755
				#
				# a=0 b=2 A=0 B=2 e=42 d=24 x=41 y=9 z=41 D=1e+06 E=174.743
				# R=10 u=0 s=2 S=0 M=0 T=0 m=10 l=1 n=10 k=1 w=1000 t=0.910239 j=3 Q=1
				# hg19.lastdb
				# Reference sequences=25 normal letters=2861343702
				# lambda=0.806198 K=0.071992
				#
				#    A  C  G  T
				# A  1 -1 -1 -1
				# C -1  1 -1 -1
				# G -1 -1  1 -1
				# T -1 -1 -1  1
				#
				# Coordinates are 0-based.  For - strand matches, coordinates
				# in the reverse complement of the 2nd sequence are used.
				#
				# name start alnSize strand seqSize alignment
				#
				# batch 0
				my $line = '';
				do {
					if (/LAST\s+version\s+(\d+)/i) {
						my $version = $1;
						printf "# \@PG_ID\t%s\n", 'lastal';
						printf "# \@PG_PN\t%s\n", 'lastal';
						printf "# \@PG_VN\t%s\n", $version;
					} elsif (/\.lastdb/) {
						chomp();
						my $lastdb = $_; $lastdb =~ s/^#\s+//;
						printf "# \@PG_DB\t%s\n", $lastdb;
					} else {
						# ignore
					}
					$_=<STDIN>;
				} while (defined $_ && '#' eq substr($_,0,1));
				print "# \@PG_END\n";
				$parseHeader = 0;
			} else {
				# Query sequences=1669
				next;
			}
		}
		if (defined $_ && /^a\s+/) {
			# let's process the block
			# a score=1897 EG2=0 E=0
			# s chr6                        160881684 3801 + 171115067 CCCAAGAAAAC
			# s WTD01:50183:2D:P:3:M:L13738        32 3791 +     13738 CCCAAGAAAAT
			# q WTD01:50183:2D:P:3:M:L13738                            4/0825;7<=3
			
			my %alignment = ();
			my $aline = $_;
			my $sline = <STDIN>;
			my $qline = <STDIN>;
			my $qualityLine = <STDIN>;
			
			my $queryId = utilities::getQueryIdValue ($qline);
			
			if ($currQueryId ne $queryId) {
				# end of previous batch of alignment for the same read
				if ('' ne $currQueryId) {
					$endTime = time();
					# process the batch before starting to work on new alignment
					$serialNumber++;
					my %pickRequest = (requestId=>$serialNumber, lines=>$currAlignmentsRef, numAlignments=>$numAlignments, ioTime=>$endTime-$startTime, maxEG2=>$maxEG2, minIdentityPercent=>$minIdentityPercent);
					$qRequests->{_hackSlots}->down(); # block until there is a free slot
					$qRequests->{qRequests}->enqueue(\%pickRequest);
					
					$endTime = time; # we could have been stucked waiting for processor, reset this wait time that affect I/O timing
				}
				
				# initialize new batch
				$startTime = $endTime;
				$numAlignments = 0;
				my @currAlignments = ();
				$currAlignmentsRef = \@currAlignments;
				$currQueryId = $queryId;
			}
			
			$numAlignments++;
			
			push @{$currAlignmentsRef}, $aline, $sline, $qline, $qualityLine;
		} else {
			# do nothing
		}
	}
	if ('' ne $currQueryId) {
		$endTime = time();
		# process the batch before starting to work on new alignment
		$serialNumber++;
		my %pickRequest = (requestId=>$serialNumber, lines=>$currAlignmentsRef, numAlignments=>$numAlignments, ioTime=>$endTime-$startTime, maxEG2=>$maxEG2, minIdentityPercent=>$minIdentityPercent);
		$qRequests->{_hackSlots}->down(); # block until there is a free slot
		$qRequests->{qRequests}->enqueue(\%pickRequest);
	}

	return $serialNumber;
}

1;

__END__
