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

package PickerMT;

BEGIN {
	$VERSION = '0.1';
}

#####
# TODO:
#

use strict;
use warnings;

use Data::Dumper;
use Storable qw(dclone freeze thaw);
use Getopt::Long;
use base 'Exporter';

use utilities;


# constants
my @G_ReportCols = ('score', 'EG2', 'E', '=', '%=', 'X', '%X', 'D', '%D', 'I', '%I', 'qStart', 'qEnd', 'qStrand', 'qALen', 'q%', 'refId', 'refStart', 'refEnd', 'refStrand', 'refALen', 'cigar');

# TODO: needed
my $G_MAX_EG2_READLEN = 1.0e-49; # this is for reasonable read length
# was: EG2<1.74743e-12; # this is for 10x coverage
my $G_MIN_IDENTITY_VARIANCE = 65;  # was: %= >=55%

# extending overlap
my $G_MIN_FRACTION_UNIQUE = 0.5;
my $G_MIN_SPECIFIC_EXTENDED_BP = 170; #200; # rescue very specific alignment for extension
my $G_FLANKING_REMAIN_MIN_FRACTION_COV_ALIGN = 0.95;
my $G_FLANKING_MAX_OVERLAP_FRACTION_COV_ALIGN = 0.5;
my $G_EXTENSION_MAX_OVERLAP_FRACTION = 0.5;
my $G_EXTENSION_MAX_OVERLAP_BP = 2500;
#my $G_EXTENSION_MIN_LENGTH_FRACTION = 0.33;
my $G_EXTENSION_MIN_LENGTH_FRACTION = 0.70;
my $G_EXTENSION_MIN_LENGTH = 900;
# bound for proximity to best seed
my $G_IDENTITY_VARIANCE = 13; #10;

# stiched results threshold
my $G_MIN_FINAL_READ_COVERAGE = 0.65; #0.7;

# superset threshold
my $G_SIMILAR_IDENTITY_VARIANCE = 5;
my $G_SIMILAR_MIN_FRACTION_COV_ALIGN = 0.95;
my $G_MAX_SUBSET_DISTANCE_BP = 50; # collapse results further

# candidates stratification
my @G_EG2Bounds = (0.0, 1e-300, 1e-200, 1e-100, 1e-50);
my $G_MAX_CANDIDATE_SEEDS = 200;
my $G_MAX_CANDIDATE_EG2_SEEDS = 100;

# runtime bound control
my $G_MAX_REPORT_CANDIDATES = 10;
my $G_MAX_EXTEND_SEED_TIME = 5 * 60; # maximum 5 min on extension
# my $G_MAX_EXTEND_SEED_TIME = 10; # maximum 5 min on extension
my $G_MAX_EXTEND_DEPTH = 5; # do not try to stitch more than 5 fragments for each extension direction
my $G_MAX_EXTEND_CANDIDATES = 100;


# TODO: put in hash!
# reporting parameters
my $G_MAX_REPORT_SUBROW = 3;
my $G_REPORT_METRICS = 1;
my $G_REPORT_IOTIME = 1;
my $G_REPORT_MAX_IOTIME_SEC = 10;
my $G_REPORT_PROCESSTIME = 1;
my $G_REPORT_MAX_PROCESSTIME_SEC = 30;


#####

sub getMaxEG2 {
	return $G_MAX_EG2_READLEN;
}

sub getMinIdentityPercentage {
	return $G_MIN_IDENTITY_VARIANCE;
}

#####
# OUTPUT:
#####

sub printAlignments {
	my ($linesRef, $alignmentsRef, $prefix, $maxRows, $maxSubRows) = @_;
    
    my $pprefix = ''; $pprefix = $prefix if (defined $prefix);
    my $numAligns = scalar(@{$alignmentsRef});
    my $numRows = (!defined $maxRows || -1==$maxRows) ? $numAligns : $maxRows;
    $numRows = $numAligns if ($numRows>$numAligns);

    for(my $i=0; $i<$numRows; ++$i) {
        my $alignRef = $alignmentsRef->[$i];
        my @cols = ();
        
        grep { push @cols, $alignRef->{$_}; } ('score', 'EG2', 'E');
        
		my $denominator = 0; grep { $denominator += $alignRef->{profile}->{$_}; } @utilities::G_AlignmentTypes;
        
        grep { push @cols, $alignRef->{profile}->{$_}; push @cols, sprintf("%.2f", $alignRef->{profile}->{$_}*100.0/$denominator)} @utilities::G_AlignmentTypes;
        
        grep { push @cols, $alignRef->{read}->{$_}; } ('readstart', 'readend', 'readstrand', 'readalignlen');
        push @cols, sprintf("%.2f", $alignRef->{read}->{readalignlen}*100.0/$alignmentsRef->[0]->{read}->{readsize});
        grep { push @cols, $alignRef->{ref}->{$_}; } ('ref', 'refstart', 'refend', 'refstrand', 'refalignlen');
        
		push @cols, $alignRef->{'read'}->{'cigar'};
		
       
        # handle the prefix
        $cols[0] = $pprefix . (exists $alignRef->{category} ? $alignRef->{category} : ' ') . $cols[0];
		_appendLine(join("\t", @cols), $linesRef);
		
        # have to subgroup those similar alignments
        # further reduce the amount of reporting needed
        # e.g. if the is samereadspan, there is no need to report subset, or similar
        #      as it is already consider as multi-loci
        # similarly, if there is is similar alignment, it is also consider as multi-loci
        # for subset, this can be sub-optimal or close enough reads
        
        my $numSames = exists $alignRef->{'samereadspan'} ? scalar(@{$alignRef->{'samereadspan'}}) : 0;
        my $numSubsets = exists $alignRef->{'subset'} ? scalar(@{$alignRef->{'subset'}}) : 0;
        my $numSimilars = exists $alignRef->{'similar'} ? scalar(@{$alignRef->{'similar'}}) : 0;
        
        if ($numSames>0) {
			_appendLine(sprintf("=(%d)",$numSames), $linesRef) if ($numSames>$maxSubRows);
            printAlignments($linesRef, $alignRef->{'samereadspan'}, '=', $maxSubRows);
        } elsif ($numSubsets>0) {
			_appendLine(sprintf("~(%d)",$numSubsets), $linesRef) if ($numSubsets>$maxSubRows);
            printAlignments($linesRef, $alignRef->{'subset'}, '~', $maxSubRows);
            
            if ($numSimilars>0) {
				_appendLine(sprintf("-(%d)",$numSimilars), $linesRef) if ($numSimilars>$maxSubRows);
                printAlignments($linesRef, $alignRef->{'similar'}, '-', $maxSubRows);
            }
        } elsif ($numSimilars>0) {
			_appendLine(sprintf("-(%d)",$numSimilars), $linesRef) if ($numSimilars>$maxSubRows);
            printAlignments($linesRef, $alignRef->{'similar'}, '-', $maxSubRows);
        }
    }
}

#####
# EXTENSION:
#####

sub getSeedsRef {
    my ($alignmentsRef, $lowestIdentity, $resultsRef, $nonSeedIndexRef) = @_;
    
    my @seeds = ($alignmentsRef->[0]);
    my $seedEG2 = $seeds[0]->{'EG2'};
    my $nonSeedIndex = 1;
    for(; $nonSeedIndex < scalar(@{$alignmentsRef}); $nonSeedIndex++) {
        my $alignRef = $alignmentsRef->[$nonSeedIndex];
        last if ($alignRef->{'EG2'} > $seedEG2);
        if ($alignRef->{profile}->{'%='}>=$lowestIdentity) {
            push @seeds, $alignRef;
        } else {
            # we exclude those which do not have a sufficiently close %identity
        }
    }
    
    @{$resultsRef} = @seeds;
    $$nonSeedIndexRef = $nonSeedIndex;
}

sub getEndOverlapped {
    my ($seedRef, $qEndsOrderRef, $candidatesRef) = @_;
    
    # we should not go beyond half the length?
    my $maxOverlapBP = int ($G_EXTENSION_MAX_OVERLAP_FRACTION * ($seedRef->{read}->{readend}-$seedRef->{read}->{readstart}+1));
    $maxOverlapBP = $G_EXTENSION_MAX_OVERLAP_BP if ($maxOverlapBP>$G_EXTENSION_MAX_OVERLAP_BP);
    my $checkBoundary = int($seedRef->{read}->{readstart} + $maxOverlapBP);
    
    my $seedStartPos = $seedRef->{read}->{readstart};
    my $seedEndRef = $qEndsOrderRef->{$seedStartPos};
    while (-1!=$seedEndRef->{nextPos}) {
        if (scalar(@{$seedEndRef->{refs}})>0) {
            # extrapolate each of them
            foreach my $candidateRef (@{$seedEndRef->{refs}}) {
                if ($candidateRef->{read}->{readstart} < $seedStartPos) {
					my ($overlap, $seedLen, $alignLen) = utilities::computeIntervalCoverage ($seedRef->{read}->{readstart}, $seedRef->{read}->{readend}, $candidateRef->{read}->{readstart}, $candidateRef->{read}->{readend});
                    if (($alignLen - $overlap) >= $G_MIN_SPECIFIC_EXTENDED_BP && $overlap<$seedLen) {
                        my $overlapFractionPerSeed = $overlap / $seedLen;
                        my $overlapFractionPerAlign = $overlap / $alignLen;
                        my $uniqueFractionPerSeed = 1.0 - $overlapFractionPerSeed;
                        my $uniqueFractionPerAlign = 1.0 - $overlapFractionPerAlign;
                        if ($uniqueFractionPerAlign>=$G_MIN_FRACTION_UNIQUE && $uniqueFractionPerSeed>=$G_MIN_FRACTION_UNIQUE) {
                            # this alignment can extend our seed
                            push @{$candidatesRef}, [$candidateRef];
                        } else {
                            # if we can close the gap
                            my $remaining = $seedRef->{read}->{readstart};
                            if (($alignLen - $overlap)/$remaining>=$G_FLANKING_REMAIN_MIN_FRACTION_COV_ALIGN && $overlapFractionPerSeed<$G_FLANKING_MAX_OVERLAP_FRACTION_COV_ALIGN) {
                                push @{$candidatesRef}, [$candidateRef];
                            }
                        }
                    }
                }
            }
        }
        
        my $nextPos = $seedEndRef->{nextPos};
        last if ($nextPos>=$checkBoundary);
        $seedEndRef = $qEndsOrderRef->{$nextPos};
    }
}

sub getEndGapped {
    my ($seedRef, $qEndsOrderRef, $candidatesRef) = @_;
    
    my $seedStartPos = $seedRef->{read}->{readstart};
    my $seedEndRef = $qEndsOrderRef->{$seedStartPos};
    while (-1!=$seedEndRef->{prevPos}) {
        if (scalar(@{$seedEndRef->{refs}})>0) {
            # extrapolate each of them
            foreach my $candidateRef (@{$seedEndRef->{refs}}) {
                if ($candidateRef->{read}->{readstart} < $seedStartPos) {
                    my ($overlap, $seedLen, $alignLen) = utilities::computeIntervalCoverage ($seedRef->{read}->{readstart}, $seedRef->{read}->{readend}, $candidateRef->{read}->{readstart}, $candidateRef->{read}->{readend});
                    if (($alignLen - $overlap) >= $G_MIN_SPECIFIC_EXTENDED_BP) {
                        my $overlapFractionPerSeed = $overlap / $seedLen;
                        my $overlapFractionPerAlign = $overlap / $alignLen;
                        my $uniqueFractionPerSeed = 1.0 - $overlapFractionPerSeed;
                        my $uniqueFractionPerAlign = 1.0 - $overlapFractionPerAlign;
                        if ($uniqueFractionPerAlign>=$G_MIN_FRACTION_UNIQUE && $uniqueFractionPerSeed>=$G_MIN_FRACTION_UNIQUE) {
                            # this alignment can extend our seed
                            push @{$candidatesRef}, [$candidateRef];
                        } else {
                            # if we can close the gap
                            my $remaining = $seedRef->{read}->{readstart};
                            if (($alignLen - $overlap)/$remaining>=$G_FLANKING_REMAIN_MIN_FRACTION_COV_ALIGN && $overlapFractionPerSeed<$G_FLANKING_MAX_OVERLAP_FRACTION_COV_ALIGN) {
                                push @{$candidatesRef}, [$candidateRef];
                            }
                        }
                    }
                }
            }
        }
        
        my $prevPos = $seedEndRef->{prevPos};
        $seedEndRef = $qEndsOrderRef->{$prevPos};
    }
}

sub extendSeedLeft {
    my ($seedRef, $qEndsOrderRef, $startTime, $candidatesRef, $depth) = @_;
    
    return -1 if ((time-$startTime)>$G_MAX_EXTEND_SEED_TIME);
    
    $depth = 0 if (!defined $depth);
    return 0 if ($depth>=$G_MAX_EXTEND_DEPTH);
    
    # if we look to the 'prev', we introduce a 'gap', but we might hit 'end' (-1)
    getEndGapped($seedRef, $qEndsOrderRef, $candidatesRef);
    return -1 if ((time-$startTime)>$G_MAX_EXTEND_SEED_TIME);
    
    # if we look to the 'next', we introduce 'overlap', but we might hit 'end' (-1) if there are only anchor!!!
    getEndOverlapped($seedRef, $qEndsOrderRef, $candidatesRef);
    return -1 if ((time-$startTime)>$G_MAX_EXTEND_SEED_TIME);
    
    # nothing has been added
    return 0 if (0==scalar(@{$candidatesRef}));
    
    @{$candidatesRef} = sort { $a->[0]->{EG2}<=>$b->[0]->{EG2} || $b->[0]->{score}<=>$a->[0]->{score} } @{$candidatesRef};

    my $bestLength = $candidatesRef->[0]->[-1]->{read}->{readalignlen};
    @{$candidatesRef} = grep { $_->[-1]->{read}->{readalignlen}/$bestLength > $G_EXTENSION_MIN_LENGTH_FRACTION || $_->[-1]->{read}->{readalignlen}>$G_EXTENSION_MIN_LENGTH } @{$candidatesRef};
	
    my $extendTimeExceeded = 0;
    my @finalCandidates = ();
    foreach my $candidateRef (@{$candidatesRef}) {
        my $leftMostSeedRef = $candidateRef->[0];
        
        if ($seedRef->{read}->{readstart}>$G_MIN_SPECIFIC_EXTENDED_BP) {
            my @leftSeedCandidates = ();
            $extendTimeExceeded = extendSeedLeft($leftMostSeedRef, $qEndsOrderRef, $startTime, \@leftSeedCandidates,$depth+1);
            if (scalar(@leftSeedCandidates)>0) {
                # we have to incorporate the candidates into the current seed
                foreach my $lefttCandidateRef (@leftSeedCandidates) {
                    my @newCandidate = ();
                    push @newCandidate, @{$candidateRef};
                    push @newCandidate, @{$lefttCandidateRef};
                    if (0==existCandidate(\@newCandidate, \@finalCandidates)) { push @finalCandidates, \@newCandidate; }
                }
            } else {
                if (0==existCandidate($candidateRef, \@finalCandidates)) { push @finalCandidates, $candidateRef; }
            }
        } else {
            if (0==existCandidate($candidateRef, \@finalCandidates)) { push @finalCandidates, $candidateRef; }
        }
        
        last if (0!=$extendTimeExceeded);
    }
    
    @{$candidatesRef} = ();
    @{$candidatesRef} = @finalCandidates;
    
    return $extendTimeExceeded;
}

sub getStartOverlapped {
    my ($seedRef, $qStartsOrderRef, $candidatesRef) = @_;
    
    # we should not go beyond half the length?
    my $maxOverlapBP = int ($G_EXTENSION_MAX_OVERLAP_FRACTION * ($seedRef->{read}->{readend}-$seedRef->{read}->{readstart}+1));
    $maxOverlapBP = $G_EXTENSION_MAX_OVERLAP_BP if ($maxOverlapBP>$G_EXTENSION_MAX_OVERLAP_BP);
    my $checkBoundary = int($seedRef->{read}->{readend} - $maxOverlapBP);
    
    my $seedEndPos = $seedRef->{read}->{readend};
    my $seedStartRef = $qStartsOrderRef->{$seedEndPos};
    while (-1!=$seedStartRef->{prevPos}) {
        if (scalar(@{$seedStartRef->{refs}})>0) {
            # extrapolate each of them
            foreach my $candidateRef (@{$seedStartRef->{refs}}) {
                if ($candidateRef->{read}->{readend} > $seedEndPos) {
                    my ($overlap, $seedLen, $alignLen) = utilities::computeIntervalCoverage ($seedRef->{read}->{readstart}, $seedRef->{read}->{readend}, $candidateRef->{read}->{readstart}, $candidateRef->{read}->{readend});
                    if (($alignLen - $overlap) >= $G_MIN_SPECIFIC_EXTENDED_BP && $overlap<$seedLen) {
                        my $overlapFractionPerSeed = $overlap / $seedLen;
                        my $overlapFractionPerAlign = $overlap / $alignLen;
                        my $uniqueFractionPerSeed = 1.0 - $overlapFractionPerSeed;
                        my $uniqueFractionPerAlign = 1.0 - $overlapFractionPerAlign;
                        if ($uniqueFractionPerAlign>=$G_MIN_FRACTION_UNIQUE && $uniqueFractionPerSeed>=$G_MIN_FRACTION_UNIQUE) {
                            # this alignment can extend our seed
                            push @{$candidatesRef}, [$candidateRef];
                        } else {
                            # if we can close the gap
                            my $remaining = $seedRef->{read}->{readsize} - $seedRef->{read}->{readend} + 1;
                            if (($alignLen - $overlap)/$remaining>=$G_FLANKING_REMAIN_MIN_FRACTION_COV_ALIGN && $overlapFractionPerSeed<$G_FLANKING_MAX_OVERLAP_FRACTION_COV_ALIGN) {
                                push @{$candidatesRef}, [$candidateRef];
                            }
                        }
                    }
                }
            }
        }
        
        my $prevPos = $seedStartRef->{prevPos};
        last if ($prevPos<=$checkBoundary);
        $seedStartRef = $qStartsOrderRef->{$prevPos};
    }
}

sub getStartGapped {
    my ($seedRef, $qStartsOrderRef, $candidatesRef) = @_;
    
    my $seedEndPos = $seedRef->{read}->{readend};
    my $seedStartRef = $qStartsOrderRef->{$seedEndPos};
    while (-1!=$seedStartRef->{nextPos}) {
        if (scalar(@{$seedStartRef->{refs}})>0) {
            # extrapolate each of them
            foreach my $candidateRef (@{$seedStartRef->{refs}}) {
                if ($candidateRef->{read}->{readend} > $seedEndPos) {
                    my ($overlap, $seedLen, $alignLen) = utilities::computeIntervalCoverage ($seedRef->{read}->{readstart}, $seedRef->{read}->{readend}, $candidateRef->{read}->{readstart}, $candidateRef->{read}->{readend});
                    if (($alignLen - $overlap) >= $G_MIN_SPECIFIC_EXTENDED_BP) {
                        my $overlapFractionPerSeed = $overlap / $seedLen;
                        my $overlapFractionPerAlign = $overlap / $alignLen;
                        my $uniqueFractionPerSeed = 1.0 - $overlapFractionPerSeed;
                        my $uniqueFractionPerAlign = 1.0 - $overlapFractionPerAlign;
                        if ($uniqueFractionPerAlign>=$G_MIN_FRACTION_UNIQUE && $uniqueFractionPerSeed>=$G_MIN_FRACTION_UNIQUE) {
                            # this alignment can extend our seed
                            push @{$candidatesRef}, [$candidateRef];
                        } else {
                            # if we can close the gap
                            my $remaining = $seedRef->{read}->{readsize} - $seedRef->{read}->{readend} + 1;
                            if (($alignLen - $overlap)/$remaining>=$G_FLANKING_REMAIN_MIN_FRACTION_COV_ALIGN && $overlapFractionPerSeed<$G_FLANKING_MAX_OVERLAP_FRACTION_COV_ALIGN) {
                                push @{$candidatesRef}, [$candidateRef];
                            }
                        }
                    }
                }
            }
        }
        
        my $nextPos = $seedStartRef->{nextPos};
        $seedStartRef = $qStartsOrderRef->{$nextPos};
    }
}

sub existCandidate {
    my ($candidateRef, $candidatesRef) = @_;
    
    return 1 if (0==scalar(@{$candidateRef}));
    
    my %candidateRowIds = ();
    grep { $candidateRowIds{$_->{rowId}} = 1 } @{$candidateRef};
    my $isUnique = 1;
    foreach my $existCandidateRef (@{$candidatesRef}) {
		my $checkRowIdsRef = dclone(\%candidateRowIds);
        grep { delete $checkRowIdsRef->{$_->{rowId}} } @{$existCandidateRef};
        $isUnique = scalar(keys %{$checkRowIdsRef}) > 0 ? 1 : 0;
        last if (0==$isUnique);
    }
    return (0==$isUnique) ? 1 : 0;
}

sub extendSeedRight {
    my ($seedRef, $qStartsOrderRef, $startTime, $candidatesRef, $depth) = @_;
    
    return -1 if ((time-$startTime)>$G_MAX_EXTEND_SEED_TIME);
    
    $depth = 0 if (!defined $depth);
    return 0 if ($depth>=$G_MAX_EXTEND_DEPTH);
    
    # if we look to the 'prev', we introduce 'overlap', but we might hit 'end' (-1) if there are only anchor!!!
    getStartOverlapped($seedRef, $qStartsOrderRef, $candidatesRef);
    return -1 if ((time-$startTime)>$G_MAX_EXTEND_SEED_TIME);
    
    # if we look to the 'next', we introduce a 'gap', but we might hit 'end' (-1)
    getStartGapped($seedRef, $qStartsOrderRef, $candidatesRef);
    return -1 if ((time-$startTime)>$G_MAX_EXTEND_SEED_TIME);
    
    # nothing has been added
    return 0 if (0==scalar(@{$candidatesRef}));

    @{$candidatesRef} = sort { $a->[0]->{EG2}<=>$b->[0]->{EG2} || $b->[0]->{score}<=>$a->[0]->{score} } @{$candidatesRef};
    
    my $bestLength = $candidatesRef->[0]->[-1]->{read}->{readalignlen};
    @{$candidatesRef} = grep { $_->[-1]->{read}->{readalignlen}/$bestLength > $G_EXTENSION_MIN_LENGTH_FRACTION || $_->[-1]->{read}->{readalignlen}>$G_EXTENSION_MIN_LENGTH} @{$candidatesRef};

    my $extendTimeExceeded = 0;
    my @finalCandidates = ();
    foreach my $candidateRef (@{$candidatesRef}) {
        my $rightMostSeedRef = $candidateRef->[-1];
        
        if ($rightMostSeedRef->{read}->{readsize}-$rightMostSeedRef->{read}->{readend} > $G_MIN_SPECIFIC_EXTENDED_BP) {
            my @rightSeedCandidates = ();
            $extendTimeExceeded = extendSeedRight($rightMostSeedRef, $qStartsOrderRef, $startTime, \@rightSeedCandidates, $depth+1);
            if (scalar(@rightSeedCandidates)>0) {
                # we have to incorporate the candidates into the current seed
                foreach my $rightCandidateRef (@rightSeedCandidates) {
                    my @newCandidate = ();
                    push @newCandidate, @{$candidateRef};
                    push @newCandidate, @{$rightCandidateRef};
                    if (0==existCandidate(\@newCandidate, \@finalCandidates)) { push @finalCandidates, \@newCandidate; }
                }
            } else {
                if (0==existCandidate($candidateRef, \@finalCandidates)) { push @finalCandidates, $candidateRef; }
            }
        } else {
            if (0==existCandidate($candidateRef, \@finalCandidates)) { push @finalCandidates, $candidateRef; }
        }
        
        last if (0!=$extendTimeExceeded);
    }
    
    @{$candidatesRef} = ();
    @{$candidatesRef} = @finalCandidates;
    
    return $extendTimeExceeded;
}

sub extendSeed {
    my ($seedRef, $qStartsOrderRef, $qEndsOrderRef, $startTime, $candidatesRef, $numRightCandidatesRef, $numLeftCandidatesRef, $reductionRef) = @_;
    
    my @rightCandidates = (); my $numRightCandidates = 0;
    my @leftCandidates = (); my $numLeftCandidates = 0;
    
    my $extendTimeExceeded = 0;
    if (($seedRef->{read}->{readsize}-$seedRef->{read}->{readend})>$G_MIN_SPECIFIC_EXTENDED_BP) {
        # need to extend 3' end ... $qStartsOrderRef
        $extendTimeExceeded = extendSeedRight($seedRef, $qStartsOrderRef, $startTime, \@rightCandidates);
        $numRightCandidates = scalar(@rightCandidates);
        #print "seedRef = ", Dumper($seedRef), "\n";
        #print "rightCandidates = ", Dumper(\@rightCandidates), "\n";
    }
	
    if ($seedRef->{read}->{readstart}>$G_MIN_SPECIFIC_EXTENDED_BP && 0==$extendTimeExceeded) {
        # need to extend 5' end ... $qEndsOrderRef
        $extendTimeExceeded = extendSeedLeft($seedRef, $qEndsOrderRef, $startTime, \@leftCandidates);
        $numLeftCandidates = scalar(@leftCandidates);
        #print "seedRef = ", Dumper($seedRef), "\n";
        #print "leftCandidates = ", Dumper(\@leftCandidates), "\n";
    }
	
    $$numRightCandidatesRef = $numRightCandidates;
    $$numLeftCandidatesRef = $numLeftCandidates;
    $$reductionRef = 0;
    if ($numRightCandidates>$G_MAX_EXTEND_CANDIDATES) {
        $numRightCandidates = $G_MAX_EXTEND_CANDIDATES;
        $$reductionRef = 1;
    }
    if ($numLeftCandidates>$G_MAX_EXTEND_CANDIDATES) {
        $numLeftCandidates = $G_MAX_EXTEND_CANDIDATES;
        $$reductionRef = 1;
    }
    # let's generate the combination
    if ($numRightCandidates>0) {
        if ($numLeftCandidates>0) {
            # both left and right extension are possible
            foreach my $leftCandidateRef (@leftCandidates[0..($numLeftCandidates-1)]) {
                foreach my $rightCandidateRef (@rightCandidates[0..($numRightCandidates-1)]) {
                    my @extendedSolution = @{$leftCandidateRef};
                    push @extendedSolution, $seedRef;
                    push @extendedSolution, @{$rightCandidateRef};
                    push @{$candidatesRef}, \@extendedSolution;
                }
            }
            foreach my $leftCandidateRef (@leftCandidates[0..($numLeftCandidates-1)]) {
                grep {
                    if (exists $_->{category}) {
                        if ('S' eq $_->{category} && (exists $seedRef->{category} && 'e' ne $seedRef->{category})) { $_->{category} = 'E'; $_->{extended} = 1; }
                    } else {
                        $_->{category} = 'e';
                    }
                } @{$leftCandidateRef};
            }
            foreach my $rightCandidateRef (@rightCandidates[0..($numRightCandidates-1)]) {
                grep {
                    if (exists $_->{category}) {
                        if ('S' eq $_->{category} && (exists $seedRef->{category} && 'e' ne $seedRef->{category})) { $_->{category} = 'E'; $_->{extended} = 1; }
                    } else {
                        $_->{category} = 'e';
                    }
                } @{$rightCandidateRef};
            }
        } else {
            # there is only right extension possible
            foreach my $rightCandidateRef (@rightCandidates[0..($numRightCandidates-1)]) {
                my @extendedSolution = @{$rightCandidateRef};
                push @extendedSolution, $seedRef;
                push @{$candidatesRef}, \@extendedSolution;
                
                grep {
                    if (exists $_->{category}) {
                        if ('S' eq $_->{category} && (exists $seedRef->{category} && 'e' ne $seedRef->{category})) { $_->{category} = 'E'; $_->{extended} = 1; }
                    } else {
                        $_->{category} = 'e';
                    }
                } @{$rightCandidateRef};
            }
        }
    } else {
        if ($numLeftCandidates>0) {
            # there is only left extension possible
            foreach my $leftCandidateRef (@leftCandidates[0..($numLeftCandidates-1)]) {
                my @extendedSolution = @{$leftCandidateRef};
                push @extendedSolution, $seedRef;
                push @{$candidatesRef}, \@extendedSolution;
                
                grep {
                    if (exists $_->{category}) {
                        if ('S' eq $_->{category} && (exists $seedRef->{category} && 'e' ne $seedRef->{category})) { $_->{category} = 'E'; $_->{extended} = 1; }
                    } else {
                        $_->{category} = 'e';
                    }
                } @{$leftCandidateRef};
            }
        } else {
            my @onlySolution = ($seedRef);
            push @{$candidatesRef}, \@onlySolution;
        }
    }
    
    return $extendTimeExceeded;
}

sub getCandidateMetrics {
    my ($candidateRef) = @_;
    
    my $cscore = 0;
    grep { $cscore+= $_->{score} } @{$candidateRef};
    
    # compute the coverage?!
    my $numCandidate = scalar(@{$candidateRef});
    my $coverageFraction = 0;
    my $startPos = 0;
    my $endPos = 0;
    if (0==$numCandidate) {
        # do nothing
    } elsif (1==$numCandidate) {
        $coverageFraction = $candidateRef->[0]->{read}->{readalignlen} / $candidateRef->[0]->{read}->{readsize};
        $startPos = $candidateRef->[0]->{read}->{readstart};
        $endPos = $candidateRef->[0]->{read}->{readend};
    } else {
        my $x1s = $candidateRef->[0]->{read}->{readstart};
        my $x1e = $candidateRef->[0]->{read}->{readend};
        $coverageFraction = $x1e - $x1s + 1;
        for(my $i=1; $i<$numCandidate; ++$i) {
            my $candidateRef = $candidateRef->[$i];
            my ($overlap, $seedLen, $alignLen) = utilities::computeIntervalCoverage ($x1s, $x1e, $candidateRef->{read}->{readstart}, $candidateRef->{read}->{readend});
            if (0==$overlap) {
                $coverageFraction += $alignLen;
                
                # consider if this is necessary: penalize the gap between region
                #my $gap = ($candidateRef->{read}->{readstart} - $x1e + 1);
                my $gap = ($x1e < $candidateRef->{read}->{readstart}) ? ($candidateRef->{read}->{readstart} - $x1e - 1) : ($x1s - $candidateRef->{read}->{readend} - 1);
                $cscore -= $gap;
            } else {
                # there are some overlap
                $coverageFraction += ($alignLen - $overlap);
                
                # consider if this is necessary, penalize the overlapping region
                $cscore -= $overlap;
            }
            $x1s = $candidateRef->{read}->{readstart} if ($candidateRef->{read}->{readstart} < $x1s);
            $x1e = $candidateRef->{read}->{readend} if ($candidateRef->{read}->{readend} > $x1e);
        }
        $coverageFraction /= $candidateRef->[0]->{read}->{readsize};
        $startPos = $x1s;
        $endPos = $x1e;
    }
    
    return ($cscore, $coverageFraction, $startPos, $endPos);
}

sub getBestCandidateMetrics {
    my ($candidatesRef) = @_;
    
    my $bestScore = 0;
    my $bestCoverageFraction = 0.0;
    my $bestStartPos = 0;
    my $bestEndPos = 0;
    foreach my $candidateRef (@{$candidatesRef}) {
        my ($bcscore, $coverageFraction, $startPos, $endPos) = getCandidateMetrics($candidateRef);
        # do we track the score and coverage separately?
        $bestScore = $bcscore if ($bcscore>$bestScore);
        if ($coverageFraction>$bestCoverageFraction) {
            $bestCoverageFraction = $coverageFraction;
            $bestStartPos = $startPos;
            $bestEndPos = $endPos;
        }
        
    }
    
    return ($bestScore, $bestCoverageFraction, $bestStartPos, $bestEndPos);
}

sub collapsedSameSpans {
    my ($candidatesRef, $resultsRef) = @_;
    
    @{$resultsRef} = ();
    my $numCandidates = scalar(@{$candidatesRef});
    if ($numCandidates>0) {
        push @{$resultsRef}, $candidatesRef->[0];
        for(my $i=1; $i<$numCandidates; ++$i) {
            my $alignRef = $candidatesRef->[$i];
            my $isExact = 0;
            foreach my $resultRef (@{$resultsRef}) {
                if ($resultRef->{read}->{readstart} == $alignRef->{read}->{readstart} && $resultRef->{read}->{readend} == $alignRef->{read}->{readend}) {
                    $resultRef->{'samereadspan'} = [] if (!exists $resultRef->{'samereadspan'});
                    push @{$resultRef->{'samereadspan'}}, $alignRef;
                    $isExact = 1; last;
                }
            }
            push @{$resultsRef}, $alignRef if (0==$isExact);
        }
    }
}

sub processSupersetSpans {
    my ($candidatesRef, $bestIdentity, $resultsRef) = @_;
    
    @{$resultsRef} = ();
    my $numCandidates = scalar(@{$candidatesRef});
    if ($numCandidates>0) {
        push @{$resultsRef}, $candidatesRef->[0];
        for(my $i=1; $i<$numCandidates; ++$i) {
            my $alignRef = $candidatesRef->[$i];
            my $isSimilar = 0;
            foreach my $resultRef (@{$resultsRef}) {
				my ($overlap, $seedLen, $alignLen) = utilities::computeIntervalCoverage ($resultRef->{read}->{readstart}, $resultRef->{read}->{readend}, $alignRef->{read}->{readstart}, $alignRef->{read}->{readend});
                # if it is a subset, it might be very similar or just really a subset
                # else, it should be left alone
                my $overlapFractionPerSeed = $overlap / $seedLen;
                my $overlapFractionPerAlign = $overlap / $alignLen;
                
                if (1==$overlapFractionPerAlign) {
                    # alignRef is a subset of resultRef
                    # we have to decide how large a subset can be for us to discard!
                    if (($alignRef->{read}->{readstart}-$resultRef->{read}->{readstart}+1<=$G_MAX_SUBSET_DISTANCE_BP)
                        && ($resultRef->{read}->{readend}-$alignRef->{read}->{readend}-+1<=$G_MAX_SUBSET_DISTANCE_BP)) {
                            $resultRef->{'subset'} = [] if (!exists $resultRef->{'subset'});
                            push @{$resultRef->{'subset'}}, $alignRef;
                            if (exists $alignRef->{samereadspan}) {
                                grep { push @{$resultRef->{'subset'}}, $_; } @{$alignRef->{samereadspan}};
                                delete $alignRef->{samereadspan};
                            }
                            $isSimilar = 1; last;
                    }
                }
                
                if ($overlapFractionPerAlign>=$G_SIMILAR_MIN_FRACTION_COV_ALIGN && $overlapFractionPerSeed>=$G_SIMILAR_MIN_FRACTION_COV_ALIGN) {
                    # no need to check '%=' as these are prefiltered already
                    if (abs($alignRef->{profile}->{'%='}-$bestIdentity)<=$G_SIMILAR_IDENTITY_VARIANCE) {
                        # possible to be similar
                        $resultRef->{'similar'} = [] if (!exists $resultRef->{'similar'});
                        push @{$resultRef->{'similar'}}, $alignRef;
                        # if alignment X is similar to alignment Y,
                        #   then all alignment Z which has the same span as alignment X will be similar to alignment Y by transitivity
                        # let's transfer that so that the reported results is not confusing
                        if (exists $alignRef->{samereadspan}) {
                            grep { push @{$resultRef->{'similar'}}, $_; } @{$alignRef->{samereadspan}};
                            delete $alignRef->{samereadspan};
                        }
                        $isSimilar = 1; last;
                    }
                }
            }
            push @{$resultsRef}, $alignRef if (0==$isSimilar);
        }
    }
}

sub isMultiLoci {
    my ($candidatesRef) = @_;
    
    return 1 if (scalar(@{$candidatesRef})>1);
    
    foreach my $candidateRef (@{$candidatesRef}) {
        foreach my $alignRef (@{$candidateRef}) {
            if (exists $alignRef->{'samereadspan'}) {
                return 1;
            }
            if (exists $alignRef->{'similar'}) {
                return 1;
            }
        }
        
    }
    
    return 0;
}

#####
# pre-processing for faster exploration
#####

sub makeStartsTraversalList {
    my ($alignsRef, $qStartsOrderRef) = @_;
    
    # let's order by qstart (asc), qend (desc)
    my @qStartsOrder = sort { $a->{read}->{readstart}<=>$b->{read}->{readstart} || $b->{read}->{readend}<=>$a->{read}->{readend} } @{$alignsRef};
    %{$qStartsOrderRef} = ();
    #my %qStartsOrder = ();
    foreach my $alignRef (@qStartsOrder) {
        my $startPos = $alignRef->{read}->{readstart};
        if (!exists $qStartsOrderRef->{$startPos}) { $qStartsOrderRef->{$startPos} = {anchors=>[], refs=>[$alignRef], prevPos=>-1, nextPos=>-1, currPos=>$startPos}; }
        else { push @{$qStartsOrderRef->{$startPos}->{refs}}, $alignRef; }
        
        # we inject the end position of the alignment here so that we can quickly home in on locus where our next extension alignment can be
        my $endPos = $alignRef->{read}->{readend};
        if (!exists $qStartsOrderRef->{$endPos}) { $qStartsOrderRef->{$endPos} = {anchors=>[$alignRef], refs=>[], prevPos=>-1, nextPos=>-1, currPos=>$endPos}; }
        else { push @{$qStartsOrderRef->{$endPos}->{anchors}}, $alignRef;}
    }
    my @qStartsLookup = sort { int($a) <=> int($b) } keys %{$qStartsOrderRef};
    my $numStartOrders = scalar(@qStartsLookup);
    for(my $i=0; $i<$numStartOrders; ++$i) {
        my $startPos = $qStartsLookup[$i];
        if (0==$i) { $qStartsOrderRef->{$startPos}->{nextPos} = $qStartsLookup[1]; }
        elsif ($i==($numStartOrders-1)) { $qStartsOrderRef->{$startPos}->{prevPos} = $qStartsLookup[$numStartOrders-2] if (($numStartOrders-2)>=0); }
        else {
            $qStartsOrderRef->{$startPos}->{prevPos} = $qStartsLookup[$i-1];
            $qStartsOrderRef->{$startPos}->{nextPos} = $qStartsLookup[$i+1];
        }
    }
}

sub makeEndsTraversalList {
    my ($alignsRef, $qEndsOrderRef) = @_;
    
    # let's order by qend (desc), qstart (asc)
    my @qEndsOrder = sort { $b->{read}->{readend}<=>$a->{read}->{readend} || $a->{read}->{readstart}<=>$b->{read}->{readstart} } @{$alignsRef};
    %{$qEndsOrderRef} = ();
    foreach my $alignRef (@qEndsOrder) {
        my $endPos = $alignRef->{read}->{readend};
        if (!exists $qEndsOrderRef->{$endPos}) { $qEndsOrderRef->{$endPos} = {anchors=>[], refs=>[$alignRef], prevPos=>-1, nextPos=>-1, currPos=>$endPos}; }
        else { push @{$qEndsOrderRef->{$endPos}->{refs}}, $alignRef; }
        
        # we inject the start position of the alignment here so that we can quickly home in on locus where our next extension alignment can be
        my $startPos = $alignRef->{read}->{readstart};
        if (!exists $qEndsOrderRef->{$startPos}) { $qEndsOrderRef->{$startPos} = {anchors=>[$alignRef], refs=>[], prevPos=>-1, nextPos=>-1, currPos=>$startPos}; }
        else { push @{$qEndsOrderRef->{$startPos}->{anchors}}, $alignRef;}
    }
    my @qEndsLookup = sort { int($a) <=> int($b) } keys %{$qEndsOrderRef};
    my $numEndOrders = scalar(@qEndsLookup);
    for(my $i=0; $i<$numEndOrders; ++$i) {
        my $endPos = $qEndsLookup[$i];
        if (0==$i) { $qEndsOrderRef->{$endPos}->{nextPos} = $qEndsLookup[1]; }
        elsif ($i==($numEndOrders-1)) { $qEndsOrderRef->{$endPos}->{prevPos} = $qEndsLookup[$numEndOrders-2] if (($numEndOrders-2)>=0); }
        else {
            $qEndsOrderRef->{$endPos}->{prevPos} = $qEndsLookup[$i-1];
            $qEndsOrderRef->{$endPos}->{nextPos} = $qEndsLookup[$i+1];
        }
    }
}

#####
# collapsing alignments for managibility
#####

sub collapseIfSameSpan {
    my ($seedsRef, $collapsedRef) = @_;
    
    @{$collapsedRef} = ();
    my $numSeeds = scalar(@{$seedsRef});
    return if (0==$numSeeds);
    
    push @{$collapsedRef}, $seedsRef->[0];
    for(my $i=1; $i < $numSeeds; $i++) {
        my $alignRef = $seedsRef->[$i];
        my $isExact = 0;
        foreach my $seedRef (@{$collapsedRef}) {
            if ($seedRef->{read}->{readstart} == $alignRef->{read}->{readstart} && $seedRef->{read}->{readend} == $alignRef->{read}->{readend}) {
                $seedRef->{'samereadspan'} = [] if (!exists $seedRef->{'samereadspan'});
                push @{$seedRef->{'samereadspan'}}, $alignRef;
                $isExact = 1; last;
            }
        }
        push @{$collapsedRef}, $alignRef if (0==$isExact);
    }
}

#####
# processing alignments collection
#####

sub _appendLine {
	my ($line, $linesRef) = @_;
	push @{$linesRef}, $line;
}

sub _extendLine {
	my ($line, $linesRef) = @_;
	$linesRef->[$#$linesRef] .= $line;
}

sub processAlignments {
	my ($currReadRef, $alignmentsRef, $unfilteredCounts, $ioTime, $alignLinesRef, $logLinesRef) = @_;
    
    my $startTime = time();
    
    my $numOfAligns = scalar(@{$alignmentsRef});
    my %reportMetrics = ('#uf_aligns'=>$unfilteredCounts, '#aligns'=>$numOfAligns, '#seeds'=>0, '#nonseeds'=>0);
    # need local time for knowing how long the processing has been
    _appendLine(sprintf("%s Processing read %s (%d,%d)", utilities::getTimeStamp(), $currReadRef->{read}, $unfilteredCounts, $numOfAligns), $logLinesRef);
	
    # cases that we are handling:
    # 1. there is no alignment
    # 2. there is only a single alignment
    # 3. multiple alignments to handle
    
    if (0==$numOfAligns) {
        # case 1: there is no alignment
		
		my @b = ("# ".$currReadRef->{readsize}, $currReadRef->{read}, "{!!!}");
		push @b, sprintf("align(%d,%d)", $reportMetrics{'#uf_aligns'}, $reportMetrics{'#aligns'}),
			sprintf("seed(%d)", $reportMetrics{'#seeds'}),
			sprintf("nonseed(%d)", $reportMetrics{'#nonseeds'}) if (0!=$G_REPORT_METRICS);
		_appendLine(join("\t", @b), $alignLinesRef);
		_appendLine("# ".join("\t", @G_ReportCols), $alignLinesRef);
		_appendLine("", $alignLinesRef);
		
    } elsif (1==$numOfAligns) {
        # case 2: there is only a single alignment

        grep { $_->{category} = 'S' } @{$alignmentsRef};
        $reportMetrics{'#seeds'} = 1;
		
		my @b = ("# ".$currReadRef->{readsize}, $currReadRef->{read}, "(1)");
		push @b, sprintf("align(%d,%d)", $reportMetrics{'#uf_aligns'}, $reportMetrics{'#aligns'}),
		sprintf("seed(%d)", $reportMetrics{'#seeds'}),
		sprintf("nonseed(%d)", $reportMetrics{'#nonseeds'}) if (0!=$G_REPORT_METRICS);
		_appendLine(join("\t", @b), $alignLinesRef);
		_appendLine("# ".join("\t", @G_ReportCols), $alignLinesRef);
		# TODO:
        printAlignments ($alignLinesRef, $alignmentsRef, undef, -1, $G_MAX_REPORT_SUBROW);
		_appendLine("", $alignLinesRef);
		
    } else {
        # case 3: multiple alignments to handle
        
		my $localAlignsRef = dclone($alignmentsRef);
		
        # sort the alignments first with smaller EG2, then score
        @{$localAlignsRef} = sort { $a->{EG2}<=>$b->{EG2} || $b->{score}<=>$a->{score} } @{$localAlignsRef};
        
        # get the average of the seed's percentage identity
		my $bestIdentity = utilities::getBestAlignmentIdentity($localAlignsRef);
        my $lowestIdentity = $bestIdentity - $G_IDENTITY_VARIANCE;
        
        # we will be stricter here
        # we exclude those which do not have a sufficiently close %identity
        @{$localAlignsRef} = grep { $_->{'EG2'}<$G_MAX_EG2_READLEN && $_->{profile}->{'%='}>=$lowestIdentity } @{$localAlignsRef};

        my @candidates = ();
        my $extendTimeExceeded = 0;
        if (scalar(@{$localAlignsRef})>0) {
            my $rowId = 0; foreach my $alignRef (@{$localAlignsRef}) { $alignRef->{rowId} = $rowId; $rowId++; }
            
            # we tag the seed first before collapse same and similar span
            my @seeds = ();
            my $nonSeedIndex = 1;
            getSeedsRef ($localAlignsRef, $lowestIdentity, \@seeds, \$nonSeedIndex);
            # update the alignment classification
            grep { $_->{category} = 'S' } @seeds;
            
            # let's collapse all alignments by read coordinates
            my @results = ();
            collapsedSameSpans ($localAlignsRef, \@results);
            @{$localAlignsRef} = @results; @results = ();
            
            # let's collapse all alignments by similar read coordinates
            processSupersetSpans ($localAlignsRef, $bestIdentity, \@results);
            @{$localAlignsRef} = @results; @results = ();
            
            # we now get the condensed seeds
            @seeds = ();
            $nonSeedIndex = 1;
            getSeedsRef ($localAlignsRef, $lowestIdentity, \@seeds, \$nonSeedIndex);
            
            my $numSeeds = scalar(@seeds);
            my $numNonSeeds = scalar(@{$localAlignsRef}) - $numSeeds;
            $reportMetrics{'#seeds'} = $numSeeds;
            $reportMetrics{'#nonseeds'} = $numNonSeeds;
			_extendLine(sprintf(" seed(%d) nonseed(%d)", $numSeeds, $numNonSeeds), $logLinesRef);
			
            my $extensionRef = $localAlignsRef; # by default, we use all alignment
            # possible alternative logics:
			#    test out seeds when #seeds>=2, if number of candidates is not sufficient, extend to non-seed alignments
            #    in additional, if there is only a single seed, we have to include non-seed alignments
            #    if there are too many seeds, we have to restrict the number of seeds tested
            if ($numSeeds>=5) {
                $extensionRef = \@seeds;
                if ($numSeeds>$G_MAX_CANDIDATE_SEEDS) {
                    my @slice = @seeds[0..($G_MAX_CANDIDATE_SEEDS-1)];
                    $extensionRef = \@slice;
                }
            } else {
                my %EG2BoundsCount = (); grep { $EG2BoundsCount{$_} = 0 } @G_EG2Bounds;
                foreach my $alignRef (@{$extensionRef}) {
                    foreach my $eg2 (@G_EG2Bounds) { if ($alignRef->{'EG2'}<=$eg2) { $EG2BoundsCount{$eg2}++; last; } }
                }
                my $bound = 0;
                my $total = 0;
                foreach my $eg2 (@G_EG2Bounds) {
                    if (($total + $EG2BoundsCount{$eg2})>$G_MAX_CANDIDATE_EG2_SEEDS) {
                        last;
                    } else {
                        $total += $EG2BoundsCount{$eg2};
                        $bound = $eg2;
                    }
                }
                my @subsets = grep { $_->{'EG2'}<=$bound } @{$extensionRef};
                $extensionRef = \@subsets;
            }
            
            # let's check that there is a good cause to try
            my ($ignoreScore, $maxCoverageFraction, $minStartPos, $maxEndPos) = getCandidateMetrics($extensionRef);
            if ($maxCoverageFraction>=$G_MIN_FINAL_READ_COVERAGE) {
                
                # let's order by qstart (asc), qend (desc)
                my %qStartsOrder = ();
                makeStartsTraversalList ($extensionRef, \%qStartsOrder);
                
                # let's order by qend (desc), qstart (asc)
                my %qEndsOrder = ();
                makeEndsTraversalList ($extensionRef, \%qEndsOrder);
                
                # let's collapse those which will be consider multiple-loci
                my @collapsedSeeds = ();
                collapseIfSameSpan($extensionRef, \@collapsedSeeds);
                
                # by now, we have two look up tables "%qStartsOrder" and "%qEndsOrder"
                # we can now iterate thru' the seeds and work out possible candidates
                # each candidates seeded by a seed should be scored
                # finally, the score should determine which are the candidates!!!
                my $bestCandidateRef = undef;
                my $bestScore = 0;
                my $bestCoverageFraction = 0;
                my $bestStartPos = 0;
                my $bestEndPos = 0;
                my $minScore = 0;
                my $minCoverageFraction = 0;
                my $extendStartTime = time;
                my $numRightCandidates = 0;
                my $numLeftCandidates = 0;
                my $candidatesReduction = 0;
                my $numExtensionDone = 0;
                foreach my $seedRef (@collapsedSeeds) {
                    # let's check each of the seed
                    next if (exists $seedRef->{extended} && 0!=$seedRef->{extended});
                    if (defined $bestCandidateRef && $bestCoverageFraction>0.9) {
                        next if (($seedRef->{read}->{readalignlen}/$seedRef->{read}->{readsize})<$minCoverageFraction);
                    }
                    
                    my @seedCandidates = ();
                    $extendTimeExceeded = extendSeed ($seedRef, \%qStartsOrder, \%qEndsOrder, $extendStartTime, \@seedCandidates, \$numRightCandidates, \$numLeftCandidates, \$candidatesReduction);
                    next if (0==scalar(@seedCandidates));
                    
                    if (0!=$candidatesReduction) {
						_extendLine(sprintf(" reduced(R,%d,%d)..", $numRightCandidates, $G_MAX_EXTEND_CANDIDATES), $logLinesRef) if ($G_MAX_EXTEND_CANDIDATES<$numRightCandidates);
						_extendLine(sprintf(" reduced(L,%d,%d)..", $numLeftCandidates, $G_MAX_EXTEND_CANDIDATES), $logLinesRef) if ($G_MAX_EXTEND_CANDIDATES<$numLeftCandidates);
                    }
                    
                    if (!defined $bestCandidateRef) {
                        ($bestScore, $bestCoverageFraction, $bestStartPos, $bestEndPos) = getBestCandidateMetrics(\@seedCandidates);
                        if ($bestCoverageFraction>=$G_MIN_FINAL_READ_COVERAGE) {
                            $bestCandidateRef = \@seedCandidates;
                            $minScore = 0.9*$bestScore; # TODO: 90% constant
                            $minCoverageFraction = $bestCoverageFraction - 0.1;
                        } else {
                            # if the best did not meet the criteria, we skip all candidates
                            # next;
                            @seedCandidates = ();
                        }
                    } else {
                        ($bestScore, $bestCoverageFraction, $bestStartPos, $bestEndPos) = getBestCandidateMetrics(\@seedCandidates);
                        if (0.9*$bestScore>$minScore) {
                            $bestCandidateRef = \@seedCandidates;
                            $minScore = 0.9*$bestScore; # TODO: 90% constant
                            $minCoverageFraction = $bestCoverageFraction - 0.1; # TODO: adjust coverage?
                        }
                    }
                    
                    foreach my $candidateRef (@seedCandidates) {
                        my ($score, $coverageFraction, $startPos, $endPos) = getCandidateMetrics($candidateRef);
                        next if ($coverageFraction<$G_MIN_FINAL_READ_COVERAGE);
                        # if we have good coverage but not as great a score, we will still like to know about this candidates
                        # thus, the 'or' expression
                        if ($score>=$minScore || $coverageFraction>=$minCoverageFraction) {
                            if ($coverageFraction/$bestCoverageFraction>=$G_SIMILAR_MIN_FRACTION_COV_ALIGN) {
                                if (0==existCandidate($candidateRef, \@candidates)) {
                                    push @candidates, $candidateRef;
                                }
                            }
                        }
                    }

                    # remove candidates that has fall out of range!
                    my @sieves = ();
                    foreach my $candidateRef (@candidates) {
                        my ($score, $coverageFraction, $startPos, $endPos) = getCandidateMetrics($candidateRef);
                        if ($score>=$minScore) {
                        #if ($score>=$minScore || $coverageFraction>=$minCoverageFraction) {
                            if ($coverageFraction/$bestCoverageFraction>=$G_SIMILAR_MIN_FRACTION_COV_ALIGN) {
                                if (0==existCandidate($candidateRef, \@sieves)) {
                                    push @sieves, $candidateRef;
                                }
                            }
                        }
                    }
                    @candidates = (); push @candidates, @sieves;
                    
                    # should have a flag for this?
                    my $currCount = scalar(@candidates);
                    if ($currCount>0) {
                        if (0!=isMultiLoci(\@candidates)) {
							_extendLine(" multi-loci detected..", $logLinesRef);
							_extendLine(" exceeded $G_MAX_EXTEND_SEED_TIME extension time..", $logLinesRef) if (0!=$extendTimeExceeded);
                            last;
                        } elsif ($currCount>$G_MAX_REPORT_CANDIDATES) {
							_extendLine(" exceeded $G_MAX_REPORT_CANDIDATES candidates..", $logLinesRef);
							_extendLine(" exceeded $G_MAX_EXTEND_SEED_TIME extension time..", $logLinesRef) if (0!=$extendTimeExceeded);
                            last;
                        }
                    } elsif (0!=$extendTimeExceeded) {
                        # time exceeded, we should stop and report the situation
                        # and any possible results?
						_extendLine(" exceeded $G_MAX_EXTEND_SEED_TIME extension time..", $logLinesRef);
                        last;
                    }
                }
                
            } else {
                # the alignments are all skewed to one-end of the read!
                # impossible to get sufficient read coverage
                # DO NOTHING as there is no canddiate possible!!!
				_extendLine(sprintf(" mpc(%.2f)", $maxCoverageFraction*100.0), $logLinesRef);
            }
        } else {
            $reportMetrics{'#seeds'} = 0;
            $reportMetrics{'#nonseeds'} = 0;
        }
		_extendLine(sprintf(" candidates(%d)", scalar(@candidates)), $logLinesRef);
		_extendLine(sprintf(" %sio(%d)", ($ioTime>$G_REPORT_MAX_IOTIME_SEC?"*":""), $ioTime), $logLinesRef) if (0!=$G_REPORT_IOTIME);
        if (0!=$G_REPORT_PROCESSTIME) {
            my $endTime = time();
            my $processTime = $endTime - $startTime;
			_extendLine(sprintf(" %sprocess(%d)", ($processTime>$G_REPORT_MAX_PROCESSTIME_SEC?"*":""), $processTime), $logLinesRef)
        }

		# let's report the candidate(s)
        my $numCandidates = scalar(@candidates);
		my @b = ("# ".$currReadRef->{readsize}, $currReadRef->{read}, ($numCandidates<=1 ? "(X)":"[ ]"));
		push @b, sprintf("align(%d,%d)", $reportMetrics{'#uf_aligns'}, $reportMetrics{'#aligns'}),
		sprintf("seed(%d)", $reportMetrics{'#seeds'}),
		sprintf("nonseed(%d)", $reportMetrics{'#nonseeds'}) if (0!=$G_REPORT_METRICS);
		push @b, "extExcess" if (0!=$extendTimeExceeded);
		_appendLine(join("\t", @b), $alignLinesRef);
		_appendLine("# ".join("\t", @G_ReportCols), $alignLinesRef);
		
		
        # note that sorting by the read coordinate only make sense when all the seeds are equal, else still better by EG2
        if ($numCandidates>1) {
            my $candidateId = 0;
            my $numCandidates = scalar(@candidates);
            foreach my $candidateRef (@candidates) {
                $candidateId++;
                @{$candidateRef} = sort { $a->{read}->{readstart}<=>$b->{read}->{readstart} } @{$candidateRef};
				_appendLine(sprintf("### candidate#%d/%d", $candidateId, $numCandidates), $alignLinesRef);
                printAlignments ($alignLinesRef, $candidateRef, undef, -1, $G_MAX_REPORT_SUBROW);
            }
        } elsif (1==$numCandidates) {
            @{$candidates[0]} = sort { $a->{read}->{readstart}<=>$b->{read}->{readstart} } @{$candidates[0]};
			_appendLine("### candidate#1/1", $alignLinesRef);
            printAlignments ($alignLinesRef, $candidates[0], undef, -1, $G_MAX_REPORT_SUBROW);
        } else {
			_appendLine("### NO candidate left", $alignLinesRef);
        }
        
		_appendLine("", $alignLinesRef);
    }
}

1;

__END__
