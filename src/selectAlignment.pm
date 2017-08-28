#####
#
# Picky - Structural Variants Pipeline for (ONT) long read
#
# Copyright (c) 2016-2017  Chee-Hong WONG  The Jackson Laboratory
#
#####

package selectAlignment;

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
use Storable qw(dclone);
use Getopt::Long;
use base 'Exporter';
our @EXPORT = qw(runSelectRepresentativeAlignments);

use utilities;
use Picker;

my $G_KEEPLINE = 0; # keey maf record lines in memory; NOTE: may consume much memory
my $G_MAX_EG2 = 1.74743e-12; # this is for 10x coverage
my $G_MIN_IDENTITY_VARIANCE = 55; # i.e. %= >=55%

my $G_USAGE = "
$0 selectRep --in <alignFile>

  --in STR       .align file
  --out STR      output filename prefix
  --start INT    First read id to be processed
  --end INT      Last read id to be processed
  --block INT    Number of reads to be processed
";

sub runSelectRepresentativeAlignments {
	my $file = undef;
	my $outfile = undef;
	my %rangeMarker = (specified=>0, start=>0, end=>0, block=>0, label=>'');
	
	GetOptions (
	"in=s"   => \$file,
	"out=s"   => \$outfile,
	"start=i" => \$rangeMarker{start},
	"end=i" => \$rangeMarker{end},
	"block=i" => \$rangeMarker{block})
	or die("Error in command line arguments\n");
	
	open INFILE, "$file" || die "Fail to open $file\n$!\n";
	
	if (!defined $outfile) {
		$outfile = $file;
		if ($outfile =~ /.maf$/) { $outfile =~ s/.maf$/.multiple/; }
		else { $outfile .= '.multiple'; }
	}
	$outfile .= rangeMarker{'label'} if (0!=utilities::makeRange(\%rangeMarker));
	
	my $logfile = $outfile; $outfile .= '.align';
	open OUTA, ">$outfile" || die "Fail to open $outfile\n$!\n";
	print OUTA "# FILE=", $file, "\n";
	
	$logfile .= '.log';
	open OUTL, ">$logfile" || die "Fail to open $logfile\n$!\n";

	my $maxEG2 = $G_MAX_EG2; #Picker::getMaxEG2();
	my $minIdentityPercent = $G_MIN_IDENTITY_VARIANCE; #Picker::getMinIdentityPercentage();
	
	my $currQueryId = '';
	my $currReadRef = undef;
	my @currAlignments = ();
	my $startTime = time;
	my $chunkStartTime = $startTime;
	my $endTime = 0;
	my $numAlignments = 0;
	my $startId = $rangeMarker{'start'};
	my $endId = $rangeMarker{'end'};
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
			
			utilities::parseMAFRecord ($aline, $sline, $qline, $qualityLine, \%alignment, $G_KEEPLINE, $currLine);
			
			if ($currQueryId ne $alignment{'read'}->{'read'}) {
				# end of previous batch of alignment for the same read
				if ('' ne $currQueryId) {
					$endTime = time();
					# process the batch before starting to work on new alignment
					Picker::processAlignments ($currReadRef, \@currAlignments, $numAlignments, $endTime-$startTime, *OUTA, *OUTL);
				}
				
				# initialize new batch
				$startTime = $endTime;
				$numAlignments = 0;
				@currAlignments = ();
				$currQueryId = $alignment{'read'}->{'read'};
				$currReadRef = $alignment{'read'};
			}
			
			$numAlignments++;
			
			next if ($alignment{'EG2'}>$maxEG2);
			next if ($alignment{'profile'}->{'%='}<$minIdentityPercent);
			
			push @currAlignments, \%alignment;
		} else {
			# do nothing
		}
	}
	close INFILE;
	if ('' ne $currQueryId) {
		$endTime = time();
		# process the batch before starting to work on new alignment
		Picker::processAlignments ($currReadRef, \@currAlignments, $numAlignments, $endTime-$startTime, *OUTA, *OUTL);
	}
	
	my $totalTime = time - $chunkStartTime;
	printf OUTL "%s Total run time = %d secs\n", utilities::getTimeStamp(), $totalTime;
	
	close OUTA;
	close OUTL;
}

1;

__END__
