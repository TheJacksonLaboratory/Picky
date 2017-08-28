#####
#
# Picky - Structural Variants Pipeline for (ONT) long read
#
# Copyright (c) 2016-2017  Chee-Hong WONG  The Jackson Laboratory
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
# WORKAROUND#1: added as institute does not have a proper multi-threaded Perl using Thread::Queue 3.12 !
use Thread::Semaphore;
# END-WORKAROUND#1
use PickerMT;

my $G_KEEPLINE = 0; # keey maf record lines in memory; NOTE: may consume much memory
my $G_MAX_EG2 = 1.74743e-12; # this is for 10x coverage
my $G_MIN_IDENTITY_VARIANCE = 55; # i.e. %= >=55%

my $G_CACHE_FOLD = 10;
#my $G_NUMBER_OF_WORKERS = 8;
my $G_NUMBER_OF_WORKERS = 1; # for testing purposes

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
			my $currLine = shift @{$linesRef};

			# let's parse the lines
			utilities::parseMAFRecord ($aline, $sline, $qline, $qualityLine, \%alignment, $G_KEEPLINE, $currLine);
			
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
	my ($qResults, $numberOfWorkers, $fhAlign, $fhLog) = @_;
	
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
				print $fhLog join("\n", @{$currResult->{logs}}), "\n" if (scalar(@{$currResult->{logs}})>0);
				print $fhAlign join("\n", @{$currResult->{aligns}}), "\n" if (scalar(@{$currResult->{aligns}})>0);
				delete $outOfOrderResults{$expectedId};
				$expectedId++;
			}
		} else {
			$outOfOrderResults{$result->{requestId}} = $result;
		}
	}
}

sub readMAFRecords {
	my ($qRequests, $file, $maxEG2, $minIdentityPercent, $startId, $endId) = @_;
	
	my $startTime = time;
	my $currQueryId = '';
	my $currReadRef = undef;
	my $currAlignmentsRef = undef;
	my $endTime = $startTime;
	my $numAlignments = 0;
	
	my $serialNumber = 0;
	open INFILE, "$file" || die "Fail to open $file\n$!\n";
	while (<INFILE>) {
		next if ('#' eq substr($_,0,1));
		if (/^a\s+/) {
			# let's process the block
			# a score=1897 EG2=0 E=0
			# s chr6                        160881684 3801 + 171115067 CCCAAGAAAAC
			# s WTD01:50183:2D:P:3:M:L13738        32 3791 +     13738 CCCAAGAAAAT
			# q WTD01:50183:2D:P:3:M:L13738                            4/0825;7<=3
			
			my $currLine = $.;
			
			my %alignment = ();
			my $aline = $_;
			my $sline = <INFILE>;
			my $qline = <INFILE>;
			my $qualityLine = <INFILE>;
			
			my $queryId = utilities::getQueryIdValue ($qline);
			next if ($queryId<$startId); # not at the read that we are interested
			last if ($queryId>$endId); # not at the read that we are interested
			
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
			
			push @{$currAlignmentsRef}, $aline, $sline, $qline, $qualityLine, $currLine;
		} else {
			# do nothing
		}
	}
	close INFILE;
	if ('' ne $currQueryId) {
		$endTime = time();
		# process the batch before starting to work on new alignment
		$serialNumber++;
		my %pickRequest = (requestId=>$serialNumber, lines=>$currAlignmentsRef, numAlignments=>$numAlignments, ioTime=>$endTime-$startTime, maxEG2=>$maxEG2, minIdentityPercent=>$minIdentityPercent);
		$qRequests->{_hackSlots}->down(); # block until there is a free slot
		$qRequests->{qRequests}->enqueue(\%pickRequest);
	}
}

my $G_USAGE = "
$0 selectRep --in <alignFile>

  --in STR       .align file
  --out STR      output filename prefix
  --thread INT   number of threads
  --preload INT  Fold of thread count to preload maf records
  --start INT    First read id to be processed
  --end INT      Last read id to be processed
  --block INT    Number of reads to be processed
";

sub runSelectRepresentativeAlignments {
	my $file = undef;
	my $outfile = undef;
	my %rangeMarker = (specified=>0, start=>0, end=>0, block=>0, label=>'');
	my $numberOfThreads = $G_NUMBER_OF_WORKERS;
	my $preloadFold = $G_CACHE_FOLD;

	GetOptions (
	"in=s"   => \$file,
	"out=s"   => \$outfile,
	"thread=i" => \$numberOfThreads,
	"preload=i" => \$preloadFold,
	"start=i" => \$rangeMarker{start},
	"end=i" => \$rangeMarker{end},
	"block=i" => \$rangeMarker{block})
	or die("Error in command line arguments\n$G_USAGE");
	
	die "Fail to open $file\n" if (! -f $file);
	
	$numberOfThreads = $G_NUMBER_OF_WORKERS if ($numberOfThreads<1);
	
	if (!defined $outfile) {
		$outfile = $file;
		if ($outfile =~ /.maf$/) { $outfile =~ s/.maf$/.multiple/; }
		else { $outfile .= '.multiple'; }
	}
	$outfile .= rangeMarker{'label'} if (0!=utilities::makeRange(\%rangeMarker));
	
	my $logfile = $outfile; $outfile .= '.align';
	open my $fhAlign, ">$outfile" || die "Fail to open $outfile\n$!\n";
	$fhAlign->autoflush(1);
	print $fhAlign "# FILE=", $file, "\n";
	
	$logfile .= '.log';
	open my $fhLog, ">$logfile" || die "Fail to open $logfile\n$!\n";
	$fhLog->autoflush(1);

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
	my $tWriter = Thread->new(\&writePickedAlignments, $qResults, $numberOfThreads, $fhAlign, $fhLog);
	
	# start processor worker threads
	for(my $i=0; $i<$numberOfThreads; ++$i) {
		my $tPicker = Thread->new(\&handleAlignmentsPicking, \%qRequests, $qResults);
	}
	
	# start maf record reader in current thread
	readMAFRecords(\%qRequests, $file, $maxEG2, $minIdentityPercent, $rangeMarker{'start'}, $rangeMarker{'end'});
	for(my $i=0; $i<$numberOfThreads; ++$i) {
		$requestSlots->down(); # block until there is a free slot
		$qRequests->enqueue("DONE");
	}
	
	# producer is done, let's signal to all threads
	
	# wait for all threads to finish!
	$_->join for threads->list;
	
	my $totalTime = time - $startTime;
	printf $fhLog "%s Total run time = %d secs\n", utilities::getTimeStamp(), $totalTime;
	
	close $fhAlign;
	close $fhLog;
}

1;

__END__
