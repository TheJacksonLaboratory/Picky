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

package xlsToVCF;

BEGIN {
	$VERSION = '0.1';
}

use strict;
use warnings;

use Data::Dumper;
use Storable qw(dclone);
use Getopt::Long;
use base 'Exporter';
our @EXPORT = qw(runXLStoVCF);

my $G_XLStoVCF_MERGE_SPAN = 1000;
my $G_XLStoVCF_NON_CONVERGENT_SPAN = 20;
my $G_XLStoVCF_NON_CONVERGENT_SPAN_NEGATE = -1 * $G_XLStoVCF_NON_CONVERGENT_SPAN;
my %G_XLStoVCF_IDPREFIXES = ('DEL'=>'D', 'INS'=>'I', 'INDEL'=>'ID', 'DUP'=>'DP', 'TTLC'=>'TL', 'TDC'=>'TD', 'TDSR'=>'TJ', 'INV'=>'IV', 'INVB'=>'IB', 'TTLC'=>'TL');
my $G_XLStoVCF_REPORT_ISVTYPES = 1;
my $G_XLStoVCF_RE = 2;

my $G_USAGE = "
$0 xls2vcf --xls <picky_xls_file> [--chrom <chromosome>] [--re <minReadsSupport>]

  --xls STR       picky SV xls file
  --chrom STR     restrict output to specified chromosomes [e.g. chr20]
  --re INT        min number of read evidence [default:$G_XLStoVCF_RE]
  --merge         window to merge SV [default: $G_XLStoVCF_MERGE_SPAN bp]
  --converge      window which SVs are considered converged concordantly [default: $G_XLStoVCF_NON_CONVERGENT_SPAN bp]
";

sub runXLStoVCF {
	my @xlsfiles = ();
	my $chromosomes = '';
	my $svtypes = '';
	my $re = $G_XLStoVCF_RE;
	my $help = 0;
	
	GetOptions (
	"xls=s"    => \@xlsfiles,
	"chrom=s"  => \$chromosomes,
	"svtype=s" => \$svtypes,
	"re=i"     => \$re,
	"merge=i"  => \$G_XLStoVCF_MERGE_SPAN,
	"converge=i" => \$G_XLStoVCF_NON_CONVERGENT_SPAN,
	"help!"    => \$help)
	or die("Error in command line arguments\n$G_USAGE");

	die "$G_USAGE" if ($help);

	die "$G_USAGE" if (0==scalar(@xlsfiles));
	die "$G_USAGE" if ($G_XLStoVCF_MERGE_SPAN<=0);
	die "$G_USAGE" if ($G_XLStoVCF_NON_CONVERGENT_SPAN<=0);

	$G_XLStoVCF_NON_CONVERGENT_SPAN_NEGATE = -1 * $G_XLStoVCF_NON_CONVERGENT_SPAN;

	my %chromosomes = ();
	grep { $chromosomes{$_} = 0; } split(/,/, $chromosomes);
	my $chromosomesRef = (0==scalar(keys %chromosomes)) ? undef : \%chromosomes;

	my %svtypes = ();
	grep { $svtypes{uc($_)} = 0; } split(/,/, $svtypes);

	@xlsfiles = split(/,/, join(',', @xlsfiles));

	# TODO: properly partition the 'runs'
	my ($run) = split(/\./, $xlsfiles[0]);
	my @runs = ($run);
	_printVCFHeaders(\@runs);

	_clusterSVsEvidences(\@xlsfiles, $chromosomesRef, $svtypes, \%svtypes, $re, \@runs);
}

sub _printVCFHeaders {
    my ($runsRef) = @_;

    print '##fileformat=VCFv4.1', "\n";
	my($day, $month, $year) = (localtime)[3,4,5];
	my $fileDate = sprintf("%04d%02d%02d", $year+1900, $month, $day);
    print '##fileDate=', $fileDate, "\n";
    print '##ALT=<ID=DEL,Description="Deletion">', "\n";
    print '##ALT=<ID=DUP,Description="Duplication">', "\n";
    print '##ALT=<ID=BND,Description="Breakend">', "\n";
    print '##ALT=<ID=INS,Description="Insertion">', "\n";
    print '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">', "\n";
    print '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">', "\n";
    print '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">', "\n";
    print '##INFO=<ID=ISVTYPE,Number=1,Type=String,Description="Internal type of structural variant">', "\n" if (0!=$G_XLStoVCF_REPORT_ISVTYPES);
    print '##INFO=<ID=BERS,Number=1,Type=String,Description="Breakend replacement string">', "\n";
    print '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">', "\n";
    print '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">', "\n";
    print '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">', "\n";
    print '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">', "\n";
    print '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Distance between the two genomic positions ( END - POS )">', "\n";
    print '##INFO=<ID=RE,Number=1,Type=Integer,Description="Number of the read evidence">', "\n";
    print '##INFO=<ID=GAP,Number=1,Type=Integer,Description="Median number of bases between the two segments of the SV, in case of an insertion this is the size of the insertion">', "\n";
    print '##INFO=<ID=NOTE,Number=1,Type=String,Description="Additional notes">', "\n";
    print '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', "\n";
    print '##FILTER=<ID=singleton,Description="Single read">', "\n";
    print '##FILTER=<ID=CIPOS,Description="The CIPOS is greater or less than ',$G_XLStoVCF_NON_CONVERGENT_SPAN,'">', "\n";
    print '##FILTER=<ID=CIEND,Description="The CIEND is greater or less than ',$G_XLStoVCF_NON_CONVERGENT_SPAN,'">', "\n";

    my @cols = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT');
    push @cols, @{$runsRef};
    print '#', join("\t", @cols), "\n";
}

sub _median {
  my ( $numbersRef ) = @_;
  my $numNumbers = scalar(@{$numbersRef});
  return 'na' if (0==$numNumbers);
  my @numbers = sort { $a<=>$b } @{$numbersRef};
  if (0==($numNumbers % 2)) {
      my $midPoint = $numNumbers / 2;
      return int(($numbers[$midPoint-1]+$numbers[$midPoint])/2);
  } else {
      return $numbers[int($numNumbers/2)];
  }
}

sub _reportClusteredSpanSVs {
	my ($svGroupType, $svsRef, $runsRef, $re, $id) = @_;

    # sort the svs
    my @svs = sort { $a->{SVChrom} cmp $b->{SVChrom} || $a->{SVStart} <=> $b->{SVStart} || $a->{SVEnd} <=> $b->{SVEnd}} @{$svsRef};
    # cluster them into block first, and iterate each block
    my @startClusters = ();
    my $startClusterIndex = 0;
    my $clusterRepStart = $svs[0]->{SVStart};
    my $clusterRepChrom = $svs[0]->{SVChrom};
    for(my $i=0; $i<scalar(@svs); ++$i) {
        if ($clusterRepChrom eq $svs[$i]->{SVChrom} && abs($svs[$i]->{SVStart}-$clusterRepStart)<=$G_XLStoVCF_MERGE_SPAN) {
            $startClusters[$startClusterIndex] = [] if (!defined $startClusters[$startClusterIndex]);
            push @{$startClusters[$startClusterIndex]}, $i;
            $clusterRepStart = $svs[$i]->{SVStart};
        } else {
            $startClusterIndex++;
            $startClusters[$startClusterIndex] = [];
            push @{$startClusters[$startClusterIndex]}, $i;
            $clusterRepStart = $svs[$i]->{SVStart};
            $clusterRepChrom = $svs[$i]->{SVChrom};
        }
    }

    # next, subclusters within each block
    my @startEndClusters = ();
    for(my $j=0; $j<scalar(@startClusters); ++$j) {
        my $startGroupRef = $startClusters[$j];
        my @endClusters = ();
        my $endClusterIndex = 0;
        # need to sort the end
        @{$startGroupRef} = sort { $svs[$a]->{SVEnd} <=> $svs[$b]->{SVEnd} } @{$startGroupRef};
        my $clusterRepEnd = $svs[$startGroupRef->[0]]->{SVEnd};
        for(my $i=0; $i<scalar(@{$startGroupRef}); ++$i) {
            if (abs($svs[$startGroupRef->[$i]]->{SVEnd}-$clusterRepEnd)<=$G_XLStoVCF_MERGE_SPAN) {
                $endClusters[$endClusterIndex] = [] if (!defined $endClusters[$endClusterIndex]);
                push @{$endClusters[$endClusterIndex]}, $startGroupRef->[$i];
                $clusterRepEnd = $svs[$startGroupRef->[$i]]->{SVEnd};
            } else {
                $endClusterIndex++;
                $endClusters[$endClusterIndex] = [];
                push @{$endClusters[$endClusterIndex]}, $startGroupRef->[$i];
                $clusterRepEnd = $svs[$startGroupRef->[$i]]->{SVEnd};
            }
        }

        # we have the proper subclusters form, equivalent to condense to a structural variant
        foreach my $groupRef (@endClusters) {
            my @t = sort { $a<=>$b } @{$groupRef};
            push @startEndClusters, \@t;
        }
    }

    # let's output the vcf results set
    for(my $j=0; $j<scalar(@startEndClusters); ++$j) {
        my $groupRef = $startEndClusters[$j];
		next if (scalar(@{$groupRef})<$re);
        my $memberRef = $svs[$groupRef->[0]];
        my $chrom = $memberRef->{SVChrom};
        my @starts = ();
        my @ends = ();
        my @rnames = ();
		my %svtypes = ();
        for(my $i=0; $i<scalar(@{$groupRef}); ++$i) {
            my $memberRef = $svs[$groupRef->[$i]];
            push @starts, $memberRef->{SVStart};
            if ('INS' eq $memberRef->{SVType}) {
                push @ends, $memberRef->{SVStart} + $memberRef->{SVSpan} - 1;
            } else {
                push @ends, $memberRef->{SVEnd};
            }
            push @rnames, $memberRef->{ReadId};

			$svtypes{$memberRef->{SVType}}++;
        }
        my $repStart = _median(\@starts);
        my $repEnd = _median(\@ends);
        my $ciStartLow = List::Util::min(@starts) - $repStart;
        my $ciStartHigh = List::Util::max(@starts) - $repStart;
        my $ciEndLow = List::Util::min(@ends) - $repEnd;
        my $ciEndHigh = List::Util::max(@ends) - $repEnd;
        my $numRNames = scalar(@rnames);

        my @cols = ();
        push @cols, $chrom; # CHROM
        push @cols, $repStart; # POS
        $id++;
        push @cols, $G_XLStoVCF_IDPREFIXES{$memberRef->{SVType}}.$id; # ID
        push @cols, 'N'; # REF

        push @cols, '<'.$svGroupType.'>';

        push @cols, '.'; # QUAL

        my $notes = '';
        if (1==$numRNames) {
            push @cols, 'singleton';
        } elsif ($ciStartLow<$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciStartHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciEndLow<$G_XLStoVCF_NON_CONVERGENT_SPAN_NEGATE || $ciEndHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN) {
            my @filters = ();
            if ($ciStartLow<$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciStartHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN) {
                push @filters, 'CIPOS';
            }
            if ($ciEndLow<$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciEndHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN) {
                push @filters, 'CIEND';
            }
            #push @cols, join(',', @filters);
            $notes = join('_', @filters);
            push @cols, 'PASS'; # TODO: >=2 ?
        } else {
            push @cols, 'PASS'; # TODO: >=2 ?
        }

        my @infos = ();
        if (0==$ciStartLow && 0==$ciStartHigh && 0==$ciEndLow && 0==$ciEndHigh && $numRNames>1) {
            push @infos, 'PRECISE';
        } else {
            push @infos, 'IMPRECISE';
        }
        push @infos, 'SVMETHOD=picky';
        push @infos, 'END='.$repEnd;
        push @infos, 'SVTYPE='.$svGroupType;
        push @infos, 'RE='.$numRNames;
        push @infos, 'RNAMES='.join(",", @rnames);
        push @infos, 'SVLEN='.($repEnd-$repStart+1);
        if ($numRNames>1) {
            push @infos, 'CIPOS='.$ciStartLow.','.$ciStartHigh;
            push @infos, 'CIEND='.$ciEndLow.','.$ciEndHigh;
        }
        push @infos, 'NOTE='.$notes;
		if (0!=$G_XLStoVCF_REPORT_ISVTYPES) {
			my @svtypes = (); grep { push @svtypes, sprintf("%s(%d)", $_, $svtypes{$_}); } (keys %svtypes);
			push @infos, 'ISVTYPE='.join(',', @svtypes);
		}
        push @cols, join(';', @infos); # INFO

        my $format = 'GT';
        push @cols, $format; # INFO

        # TODO: handle each run in the future
        grep { push @cols, './.'; } @{$runsRef};

        # ready for output
        print join("\t", @cols), "\n";
    }

	return $id;
}

sub _reportClusteredBreakpointSVs {
	my ($svGroupType, $allSVsRef, $runsRef, $re, $id) = @_;

	my @svTypes = ();
	my %svTypes = ();
	foreach my $svRef (@{$allSVsRef}) {
		if (!exists $svTypes{$svRef->{SVType}}) {
			$svTypes{$svRef->{SVType}} = [];
			push @svTypes, $svRef->{SVType};
		}
		push @{$svTypes{$svRef->{SVType}}}, $svRef;
	}

	foreach my $svType (@svTypes) {
		my $svsRef = $svTypes{$svType};

		# sort the svs
		my @svs = sort {$a->{SVChrom1} cmp $b->{SVChrom1} || $a->{SVChrom2} cmp $b->{SVChrom2} || $a->{SVPos1} <=> $b->{SVPos1} ||  $a->{SVPos2} <=> $b->{SVPos2}} @{$svsRef};
		# cluster them into block first, and iterate each block
		my @startClusters = ();
		my $startClusterIndex = 0;
		my $clusterRepPos = $svs[0]->{SVPos1};
		my $clusterRepChrom1 = $svs[0]->{SVChrom1};
		my $clusterRepChrom2 = $svs[0]->{SVChrom2};
		# TODO: what's the appropriate merging distance for breakpoint???
		for(my $i=0; $i<scalar(@svs); ++$i) {
			if ($clusterRepChrom1 eq $svs[$i]->{SVChrom1} && $clusterRepChrom2 eq $svs[$i]->{SVChrom2} && abs($svs[$i]->{SVPos1}-$clusterRepPos)<=$G_XLStoVCF_MERGE_SPAN) {
				$startClusters[$startClusterIndex] = [] if (!defined $startClusters[$startClusterIndex]);
				push @{$startClusters[$startClusterIndex]}, $i;
				$clusterRepPos = $svs[$i]->{SVPos1};
			} else {
				$startClusterIndex++;
				$startClusters[$startClusterIndex] = [];
				push @{$startClusters[$startClusterIndex]}, $i;
				$clusterRepPos = $svs[$i]->{SVPos1};
				$clusterRepChrom1 = $svs[$i]->{SVChrom1};
				$clusterRepChrom2 = $svs[$i]->{SVChrom2};
			}
		}

		# next, subclusters within each block
		my @startEndClusters = ();
		for(my $j=0; $j<scalar(@startClusters); ++$j) {
			my $startGroupRef = $startClusters[$j];
			my @endClusters = ();
			my $endClusterIndex = 0;
			# need to sort the end
			@{$startGroupRef} = sort { $svs[$a]->{SVPos2} <=> $svs[$b]->{SVPos2} } @{$startGroupRef};
			my $clusterRepPos2 = $svs[$startGroupRef->[0]]->{SVPos2};
			# TODO: what's the appropriate merging distance for breakpoint???
			for(my $i=0; $i<scalar(@{$startGroupRef}); ++$i) {
				if (abs($svs[$startGroupRef->[$i]]->{SVPos2}-$clusterRepPos2)<=$G_XLStoVCF_MERGE_SPAN) {
					$endClusters[$endClusterIndex] = [] if (!defined $endClusters[$endClusterIndex]);
					push @{$endClusters[$endClusterIndex]}, $startGroupRef->[$i];
					$clusterRepPos2 = $svs[$startGroupRef->[$i]]->{SVPos2};
				} else {
					$endClusterIndex++;
					$endClusters[$endClusterIndex] = [];
					push @{$endClusters[$endClusterIndex]}, $startGroupRef->[$i];
					$clusterRepPos2 = $svs[$startGroupRef->[$i]]->{SVPos2};
				}
			}

			# we have the proper subclusters form, equivalent to condense to a structural variant
			foreach my $groupRef (@endClusters) {
				my @t = sort { $a<=>$b } @{$groupRef};
				push @startEndClusters, \@t;
			}
		}

		# let's output the vcf results set
		for(my $j=0; $j<scalar(@startEndClusters); ++$j) {
			my $groupRef = $startEndClusters[$j];
			next if (scalar(@{$groupRef})<$re);
			my $memberRef = $svs[$groupRef->[0]];
			my $chrom1 = $memberRef->{SVChrom1};
			my $chrom2 = $memberRef->{SVChrom2};
			# TODO: do the mixed strands matter?
			my $strand1 = $memberRef->{SVStrand1};
			my $strand2 = $memberRef->{SVStrand2};
			my @pos1s = ();
			my @pos2s = ();
			my @rnames = ();
			my %svtypes = ();
			for(my $i=0; $i<scalar(@{$groupRef}); ++$i) {
				my $memberRef = $svs[$groupRef->[$i]];
				push @pos1s, $memberRef->{SVPos1};
				push @pos2s, $memberRef->{SVPos2};
				push @rnames, $memberRef->{ReadId};
				$svtypes{$memberRef->{SVType}}++;
			}
			my $repPos1 = _median(\@pos1s);
			my $repPos2 = _median(\@pos2s);
			my $ciStartLow = List::Util::min(@pos1s) - $repPos1;
			my $ciStartHigh = List::Util::max(@pos1s) - $repPos1;
			my $ciEndLow = List::Util::min(@pos2s) - $repPos2;
			my $ciEndHigh = List::Util::max(@pos2s) - $repPos2;
			my $numRNames = scalar(@rnames);

			# print flanking#1 of breakpoint
			my @cols = ();
			push @cols, $chrom1; # CHROM
			push @cols, $repPos1; # POS
			$id++;
			my $mateId1 = $G_XLStoVCF_IDPREFIXES{$memberRef->{SVType}}.$id;
			$id++;
			my $mateId2 = $G_XLStoVCF_IDPREFIXES{$memberRef->{SVType}}.$id;
			push @cols, $mateId1; # ID
			push @cols, 'N'; # REF

			my $partP2 = '';
			if ('+' eq $strand2) { $partP2 = sprintf("[%s:%d[N", $chrom2, $repPos2); } else { $partP2 = sprintf("N]%s:%d]", $chrom2, $repPos2); }
			push @cols, $partP2; # format

			push @cols, '.'; # QUAL

			my $notes = '';
			if (1==$numRNames) {
				push @cols, 'singleton';
			} elsif ($ciStartLow<$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciStartHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciEndLow<$G_XLStoVCF_NON_CONVERGENT_SPAN_NEGATE || $ciEndHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN) {
				my @filters = ();
				if ($ciStartLow<$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciStartHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN) {
					push @filters, 'CIPOS';
				}
				if ($ciEndLow<$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciEndHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN) {
					push @filters, 'CIEND';
				}
				#push @cols, join(',', @filters);
				$notes = join('_', @filters);
				push @cols, 'PASS'; # TODO: >=2 ?
			} else {
				push @cols, 'PASS'; # TODO: >=2 ?
			}

			my @infos = ();
			if (0==$ciStartLow && 0==$ciStartHigh && 0==$ciEndLow && 0==$ciEndHigh && $numRNames>1) {
				push @infos, 'PRECISE';
			} else {
				push @infos, 'IMPRECISE';
			}
			push @infos, 'SVMETHOD=picky';
			push @infos, 'SVTYPE='.$svGroupType;
			push @infos, 'MATEID='.$mateId2;
			push @infos, 'RE='.$numRNames;
			push @infos, 'RNAMES='.join(",", @rnames);
			if ($numRNames>1) {
				push @infos, 'CIPOS='.$ciStartLow.','.$ciStartHigh;
				push @infos, 'CIEND='.$ciEndLow.','.$ciEndHigh;
			}
			push @infos, 'NOTE='.$notes;
			if (0!=$G_XLStoVCF_REPORT_ISVTYPES) {
				my @svtypes = (); grep { push @svtypes, sprintf("%s(%d)", $_, $svtypes{$_}); } (keys %svtypes);
				push @infos, 'ISVTYPE='.join(',', @svtypes);
			}
			push @infos, 'BERS='.$partP2;
			push @cols, join(';', @infos); # INFO

			my $format = 'GT';
			push @cols, $format; # INFO

			# TODO: handle each run in the future
			grep { push @cols, './.'; } @{$runsRef};

			# ready for output
			print join("\t", @cols), "\n";

			# END - flanking#1 of breakpoints


			# print flanking#2 of breakpoint
			@cols = ();
			push @cols, $chrom2; # CHROM
			push @cols, $repPos2; # POS
			#$id++;
			push @cols, $mateId2; # ID
			push @cols, 'N'; # REF

			my $partP1 = '';
			if ('+' eq $strand1) { $partP1 = sprintf("N]%s:%d]", $chrom1, $repPos1); } else { $partP1 = sprintf("[%s:%d[N", $chrom1, $repPos1); }
			push @cols, $partP1; # format

			push @cols, '.'; # QUAL

			$notes = '';
			if (1==$numRNames) {
				push @cols, 'singleton';
			} elsif ($ciStartLow<$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciStartHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciEndLow<$G_XLStoVCF_NON_CONVERGENT_SPAN_NEGATE || $ciEndHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN) {
				my @filters = ();
				if ($ciStartLow<$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciStartHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN) {
					push @filters, 'CIPOS';
				}
				if ($ciEndLow<$G_XLStoVCF_NON_CONVERGENT_SPAN || $ciEndHigh>$G_XLStoVCF_NON_CONVERGENT_SPAN) {
					push @filters, 'CIEND';
				}
				#push @cols, join(',', @filters);
				$notes = join('_', @filters);
				push @cols, 'PASS'; # TODO: >=2 ?
			} else {
				push @cols, 'PASS'; # TODO: >=2 ?
			}

			@infos = ();
			if (0==$ciStartLow && 0==$ciStartHigh && 0==$ciEndLow && 0==$ciEndHigh && $numRNames>1) {
				push @infos, 'PRECISE';
			} else {
				push @infos, 'IMPRECISE';
			}
			push @infos, 'SVMETHOD=picky';
			push @infos, 'SVTYPE='.$svGroupType;
			push @infos, 'MATEID='.$mateId1;
			push @infos, 'RE='.$numRNames;
			push @infos, 'RNAMES='.join(",", @rnames);
			if ($numRNames>1) {
				push @infos, 'CIPOS='.$ciStartLow.','.$ciStartHigh;
				push @infos, 'CIEND='.$ciEndLow.','.$ciEndHigh;
			}
			push @infos, 'NOTE='.$notes;
			if (0!=$G_XLStoVCF_REPORT_ISVTYPES) {
				my @svtypes = (); grep { push @svtypes, sprintf("%s(%d)", $_, $svtypes{$_}); } (keys %svtypes);
				push @infos, 'ISVTYPE='.join(',', @svtypes);
			}
			push @infos, 'BERS='.$partP1;
			push @cols, join(';', @infos); # INFO

			$format = 'GT';
			push @cols, $format; # INFO

			# TODO: handle each run in the future
			grep { push @cols, './.'; } @{$runsRef};

			# ready for output
			print join("\t", @cols), "\n";

			# END - flanking#2 of breakpoints
		}
	}

	return $id;
}

sub _reportClusteredSVs {
	my ($svGroupType, $svsRef, $runsRef, $re, $id) = @_;

	if ('BND' eq $svGroupType) {
		return _reportClusteredBreakpointSVs($svGroupType, $svsRef, $runsRef, $re, $id);
	} else {
		return _reportClusteredSpanSVs($svGroupType, $svsRef, $runsRef, $re, $id);
	}
}

sub _clusterSVsEvidences {
    my ($filesRef, $chromosomesRef, $svtypes, $svtypesRef, $re, $runsRef) = @_;

	my %svByTypes = ();
    foreach my $file (@{$filesRef}) {
        open INFILE, $file || die "Fail to open $file\n$!\n";
		print STDERR "Processing $file.." ;
		my ($block, $reportUnitCount, $reportUnit) = (5000, 1000, 'K');
        my $header = <INFILE>; chomp($header); $header =~ s/^#//;
        my @headers = split(/\t/, $header);
        my $numHeaders = scalar(@headers);
        while (<INFILE>) {
			next if ('#' eq substr($_,0,1));
			
			if (0==($. % $block)) {
				printf STDERR " %d%s..", int($. / $reportUnitCount), $reportUnit;
			}
            chomp();
            my @bits = split(/\t/);
            my %item = ();
            for(my $i=0; $i<$numHeaders; ++$i) {
                $item{$headers[$i]} = $bits[$i];
            }

            next if (defined $chromosomesRef && !exists $chromosomesRef->{$item{SVChrom}});

			# DEL-ok implemented .. <DEL>
			# INS-ok implemented .. <INS>
			# TDC-ok, TDSR-ok implemented .. <DUP>
			# INDEL-ok; implemented to report either insertion and deletion according to the net effect
			# TTLC-implemented .. BND
			# INV-implemented .. <INV>
			# INVB-implemented .. BND
            # adjustment for alignment resutls
			my $svGroupType = '';
            if ('DEL' eq $item{SVType}) {
				$svGroupType = 'DEL';
                my $net = $item{sDiff} - $item{qDiff};
                if ($net > 0) {
                    if ($item{sDiff}>0) {
                        if ($item{qDiff}>0) {
                            # A3: +ve, +ve
                            # no adjustment needed
                        } else {
                            # A2: +ve, -ve
                            $item{SVStart} = ('+' eq $item{refStrand_1}) ? $item{refEnd_1} : $item{refEnd_2};
                            $item{SVStart} += $item{qDiff};
                        }
                    } else {
                        if ($item{qDiff}>0) {
                            # insertion
                        } else {
                            # C2: -ve, -ve
                            $item{SVStart} = ('+' eq $item{refStrand_1}) ? $item{refStart_2} : $item{refStart_1};
                            $item{SVEnd} = $item{SVStart} + $net;
                        }
                    }
                }
                #if ($item{qDiff}<0) {
                #    $item{SVStart} = ('+' eq $item{refStrand_1}) ? $item{refEnd_1} : $item{refEnd_2};
                #    $item{SVStart} += $item{qDiff};
                #    $item{SVStart} -= $item{sDiff} if ($item{sDiff}<0);
                #    $item{SVEnd} = $item{SVStart} + $item{SVSpan};
                #}
				$svByTypes{$svGroupType} = [] if (!exists $svByTypes{$svGroupType});
				push @{$svByTypes{$svGroupType}}, \%item;
            } elsif ('INS' eq $item{SVType}) {
				$svGroupType = 'INS';
				$svByTypes{$svGroupType} = [] if (!exists $svByTypes{$svGroupType});
				push @{$svByTypes{$svGroupType}}, \%item;
            } elsif ('TDC' eq $item{SVType}) {
				$svGroupType = 'DUP';
				$svByTypes{$svGroupType} = [] if (!exists $svByTypes{$svGroupType});
				push @{$svByTypes{$svGroupType}}, \%item;
            } elsif ('TDSR' eq $item{SVType}) {
				$svGroupType = 'DUP';
				$svByTypes{$svGroupType} = [] if (!exists $svByTypes{$svGroupType});
				push @{$svByTypes{$svGroupType}}, \%item;
            } elsif ('INDEL' eq $item{SVType}) {
				# TODO: confirm how community interpret this class!
                if ($item{sDiff}>$item{qDiff}) {
                    # net deletion against reference
                    $item{SVEnd} -= ($item{sDiff}-$item{qDiff});
					$svGroupType = 'DEL';
					$svByTypes{$svGroupType} = [] if (!exists $svByTypes{$svGroupType});
					push @{$svByTypes{$svGroupType}}, \%item;
                } elsif ($item{sDiff}<$item{qDiff}) {
                    # net insertion against reference
                    $item{SVSpan} = ($item{qDiff}-$item{sDiff});
					$svGroupType = 'INS';
					$svByTypes{$svGroupType} = [] if (!exists $svByTypes{$svGroupType});
					push @{$svByTypes{$svGroupType}}, \%item;
                }
            } elsif ('INV' eq $item{SVType}) {
				$svGroupType = 'INV';
				$svByTypes{$svGroupType} = [] if (!exists $svByTypes{$svGroupType});
				push @{$svByTypes{$svGroupType}}, \%item;
            } elsif ('INVB' eq $item{SVType}) {
				$svGroupType = 'BND';
				$svByTypes{$svGroupType} = [] if (!exists $svByTypes{$svGroupType});
				push @{$svByTypes{$svGroupType}}, \%item;
            } elsif ('TTLC' eq $item{SVType}) {
				$svGroupType = 'BND';
				$svByTypes{$svGroupType} = [] if (!exists $svByTypes{$svGroupType});
				push @{$svByTypes{$svGroupType}}, \%item;
            } else {
				die "Yet to implement SVType for record $_\n";
			}
        }
		printf STDERR " %d.. done.\n", $.;
        close INFILE;
    }

    return if (0==scalar(keys %svByTypes));

    my $id = 0;
	foreach my $svGroupType (keys %svByTypes) {
		my $svsRef = $svByTypes{$svGroupType};
		my $nextId = _reportClusteredSVs($svGroupType, $svsRef, $runsRef, $re, $id);
		$id = $nextId;
	}
}

1;

__END__
