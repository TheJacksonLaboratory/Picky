#####
#
# Picky - Structural Variants Pipeline for (ONT) long read
#
# Created Aug 16, 2016
# Copyright (c) 2016-2017  Chee-Hong WONG
#                          Genome Technologies
#                          The Jackson Laboratory
#
#####

package utilities;

use strict;
use warnings;

our $VERSION = '0.1';

use base 'Exporter';

#####

sub makeRange {
	my ($rangerRef) = @_;
	
	my $MAX_end = 1e12;
	if ($rangerRef->{start}>0) {
		if ($rangerRef->{block}>0) {
			if ($rangerRef->{end}>0) {
				# start defined, end defined, block defined : marking
				# ONLY <start>-<end>
				$rangerRef->{specified} = 1;
			} else {
				# start defined, end undef, block defined : marking
				$rangerRef->{end} = $rangerRef->{start} + $rangerRef->{block} - 1; $rangerRef->{specified} = 1;
			}
		} else {
			if ($rangerRef->{end}>0) {
				# start defined, end defined, block undef : marking
				$rangerRef->{specified} = 1;
			} else {
				# start defined, end undef, block undef : marking
				$rangerRef->{end} = $MAX_end; $rangerRef->{specified} = 1;
			}
		}
	} else {
		if ($rangerRef->{block}>0) {
			if ($rangerRef->{end}>0) {
				# start undef, end defined, block defined : marking
				$rangerRef->{start} = $rangerRef->{end} - $rangerRef->{block} + 1; $rangerRef->{specified} = 1;
			} else {
				# start undef, end undef, block defined : marking
				$rangerRef->{start} = 1; $rangerRef->{end} = $rangerRef->{block}; $rangerRef->{specified} = 1;
			}
		} else {
			if ($rangerRef->{end}>0) {
				# start undef, end defined, block undef : marking
				$rangerRef->{start} = 1; $rangerRef->{specified} = 1;
			} else {
				# start undef, end undef, block undef : everything
				$rangerRef->{start} = 1; $rangerRef->{end} = $MAX_end; $rangerRef->{specified} = 0;
			}
		}
	}
	
	#'.'.$startId.'-'.$endId
	$rangerRef->{label} = sprintf(".%d-%d", $rangerRef->{start}, $rangerRef->{end}) if (0!=$rangerRef->{specified});
	return $rangerRef->{specified};
}

sub getQueryIdValue {
	my ($qline) = @_;
	
	chomp($qline);
	my @bits = split(/\s+/, $qline);
	#WCH: gave up internal serial run id to work with PacBio for now
	#@bits = split(/:/, $bits[1]);
	#return int($bits[1]);
	return $bits[1];
}

sub getTimeStamp {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	return sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
}

#####

our @G_AlignmentTypes = ('=', 'X', 'D', 'I');

sub getBestAlignmentIdentity {
	my ($alignmentsRef) = @_;
	my $averageIdentity = $alignmentsRef->[0]->{profile}->{'='}*100.0/$alignmentsRef->[0]->{profile}->{'TOTAL'};
	return $averageIdentity;
}

sub _unused_getAverageIdentity {
	my ($alignmentsRef) = @_;
	
	my $numAlignments = scalar(@{$alignmentsRef});
	return 0.0 if (0==$numAlignments);
	
	my $averageIdentity = 0;
	foreach my $alignRef (@{$alignmentsRef}) {
		$averageIdentity += $alignRef->{profile}->{'='}*100.0/$alignRef->{profile}->{'TOTAL'};
	}
	$averageIdentity /= $numAlignments;
	
	return $averageIdentity;
}

sub _profileTwoSeqs {
	my ($refSeq, $querySeq, $profileRef) = @_;
	
	# counting is w.r.t. Reference
	# ref: A, query: C, X
	# ref: A, query: A, =
	# ref: -, query: G, I
	# ref: T, query: -, D
	$profileRef->{'X'} = 0;
	$profileRef->{'='} = 0;
	$profileRef->{'I'} = 0;
	$profileRef->{'D'} = 0;
	#$profileRef->{'Ifrag'} = 0;
	#$profileRef->{'Dfrag'} = 0;
	die "Reference length (",length($refSeq),") != Query length (",length($querySeq),")\n" if (length($refSeq)!=length($querySeq));
	for(my $i=0; $i<length($refSeq); $i++) {
		my $refbase = substr($refSeq, $i, 1);
		my $querybase = substr($querySeq, $i, 1);
		
		if ('-' eq $refbase) {
			if ('-' eq $querybase) {
				# R:-,Q:-
				die "Both reference and query contain '-' @ ", $i, "\n";
			} else {
				# R:-,Q:[ACGT]
				$profileRef->{'I'}++;
			}
		} else {
			if ('-' eq $querybase) {
				# R:[ACGT],Q:-
				$profileRef->{'D'}++;
			} else {
				# R:[ACGT],Q:[ACGT]
				if ((lc $refbase) eq (lc $querybase)) {
					$profileRef->{'='}++;
				} else {
					$profileRef->{'X'}++;
				}
			}
		}
	}
	
	# pre-compute common required value
	my $denominator = 0; grep { $denominator += $profileRef->{$_}; } @G_AlignmentTypes;
	$profileRef->{'TOTAL'} = $denominator;
	$profileRef->{'%='} = $profileRef->{'='} * 100.0 / $profileRef->{'TOTAL'};
}

sub parseMAFRecord {
	my ($aline, $sline, $qline, $qualityLine, $alignRef, $keepline) = @_;
	
	# let's process the block
	# a score=1897 EG2=0 E=0
	# s chr6                        160881684 3801 + 171115067 CCCAAGAAAAC
	# s WTD01:50183:2D:P:3:M:L13738        32 3791 +     13738 CCCAAGAAAAT
	# q WTD01:50183:2D:P:3:M:L13738                            4/0825;7<=3
	
	my $recordLine = (defined $keepline) ? $keepline : 0;
	
	%{$alignRef} = ();
	
	# work on ^a line
	chomp($aline);
	$alignRef->{lines} = [$aline] if (0!=$recordLine);
	my ($line) = $aline =~ /^a\s+(.*)/;
	foreach my $keyValue (split(/\s+/, $line)) {
		my ($key, $value) = split(/\=/, $keyValue);
		$alignRef->{$key} = $value;
	}
	
	# work on ^s line
	chomp($sline);
	push @{$alignRef->{lines}}, $sline if (0!=$recordLine);
	my @bits = split(/\s+/, $sline);
	my %ref=();
	$ref{'ref'} = $bits[1];
	$ref{'refstart'} = int($bits[2]); # it is 0-based
	$ref{'refalignlen'} = int($bits[3]);
	$ref{'refstrand'} = $bits[4];
	$ref{'refsize'} = int($bits[5]);
	$ref{'refseq'} = $bits[6];
	#my $sbases = $bits[6]; $sbases =~ s/\-//g; $ref{'refend'} = $ref{'refstart'} + length($sbases); # start is 0-based
	$ref{'refend'} = $ref{'refstart'} + $ref{'refalignlen'}; # start is 0-based
	$alignRef->{'ref'} = \%ref;
	
	# work on ^s line
	chomp($qline);
	push @{$alignRef->{lines}}, $qline if (0!=$recordLine);
	@bits = split(/\s+/, $qline);
	my %read=();
	$read{'read'} = $bits[1];
	$read{'readstart'} = int($bits[2]); # it is 0-based
	$read{'readalignlen'} = int($bits[3]);
	$read{'readstrand'} = $bits[4];
	$read{'readsize'} = int($bits[5]);
	$read{'readseq'} = $bits[6];
	#my $qbases = $bits[6]; $qbases =~ s/\-//g; $read{'readend'} = $read{'readstart'} + length($qbases); # start is 0-based
	$read{'readend'} = $read{'readstart'} + $read{'readalignlen'}; # start is 0-based
	$alignRef->{'read'} = \%read;
	
	#####
	# we need to make read as the reference '+' coordinates for easier assembly later
	if ('-' eq $read{'readstrand'}) {
		$ref{'refstrand'} = '-';
		$read{'readstrand'} = '+';
		
		$read{'readstart'} = $read{'readsize'} - $read{'readstart'} - $read{'readalignlen'} - 1; # start is 0-based
		$read{'readend'} = $read{'readstart'} + $read{'readalignlen'} + 1; # +1 'cos start is 0-based
		
		# NOTE: we do not reverse complement the sequences as the profile computation is the same
		#       however, the CIGAR generation should be taken care of
	}
	
	my %profile = ();
	_profileTwoSeqs($ref{'refseq'}, $read{'readseq'}, \%profile);
	$alignRef->{'profile'} = \%profile;
	
	# work on ^q line
	chomp($qualityLine);
	push @{$alignRef->{lines}}, $qualityLine if (0!=$recordLine);
	
	
	# clear sequences to reduce memory consumption and easier debugging
	delete $ref{'refseq'};
	delete $read{'readseq'};
}

sub generateCIGAR {
	my ($refSeq, $querySeq, $qStart, $qSize, $qStrand, $cigarRef) = @_;
	
	die "Reference length (",length($refSeq),") != Query length (",length($querySeq),")\n" if (length($refSeq)!=length($querySeq));
	
	# WCH: 20170215, qStart is passed as 0-based
	# $qStart;
	
	my @cigarBits = ();
	if ($qStart>0) {
		# there is soft clipping needed
		push @cigarBits, {type=>'S', count=>$qStart};
	}
	my $type = ''; my $count = 0;
	my $currType = '';
	my $endPos = $qStart;
	for(my $i=0; $i<length($refSeq); $i++) {
		my $refbase = substr($refSeq, $i, 1);
		my $querybase = substr($querySeq, $i, 1);
		
		if ('-' eq $refbase) {
			if ('-' eq $querybase) {
				# R:-,Q:-
				die "Both reference and query contain '-' @ ", $i, "\n";
			} else {
				# R:-,Q:[ACGT]
				$currType = 'I';
				$endPos++;
			}
		} else {
			if ('-' eq $querybase) {
				# R:[ACGT],Q:-
				$currType = 'D';
			} else {
				# R:[ACGT],Q:[ACGT]
				if ((lc $refbase) eq (lc $querybase)) {
					$currType = '=';
				} else {
					$currType = 'X';
				}
				$endPos++;
			}
		}
		
		if ($currType eq $type) {
			$count++;
		} else {
			# need to flush
			push @cigarBits, {type=>$type, count=>$count} if ($type ne '');
			# reset to new type
			$type = $currType; $count = 1;
		}
	}
	# last information is not flushed yet
	push @cigarBits, {type=>$type, count=>$count} if ($type ne '');
	
	# just in case we have clipping at the right end
	my $remaining = $qSize - $endPos;
	push @cigarBits, {type=>'S', count=>$remaining} if ($remaining>0);
	
	@{$cigarRef} = (); push @{$cigarRef}, @cigarBits;
}

sub getCIGARFromBits {
	my ($cigarRef, $reverse) = @_;
	my $cigar = '';
	if (defined $reverse && 0!=$reverse) {
		grep { $cigar .= $_->{count}.$_->{type} } reverse @{$cigarRef};
	} else {
		grep { $cigar .= $_->{count}.$_->{type} } @{$cigarRef};
	}
	return $cigar;
}

sub parseMAFRecordGenCigar {
	my ($aline, $sline, $qline, $qualityLine, $alignRef, $keepline) = @_;
	
	# let's process the block
	# a score=1897 EG2=0 E=0
	# s chr6                        160881684 3801 + 171115067 CCCAAGAAAAC
	# s WTD01:50183:2D:P:3:M:L13738        32 3791 +     13738 CCCAAGAAAAT
	# q WTD01:50183:2D:P:3:M:L13738                            4/0825;7<=3
	
	my $recordLine = (defined $keepline) ? $keepline : 0;
	
	%{$alignRef} = ();
	
	# work on ^a line
	chomp($aline);
	$alignRef->{lines} = [$aline] if (0!=$recordLine);
	my ($line) = $aline =~ /^a\s+(.*)/;
	foreach my $keyValue (split(/\s+/, $line)) {
		my ($key, $value) = split(/\=/, $keyValue);
		$alignRef->{$key} = $value;
	}
	
	# work on ^s line
	chomp($sline);
	push @{$alignRef->{lines}}, $sline if (0!=$recordLine);
	my @bits = split(/\s+/, $sline);
	my %ref=();
	$ref{'ref'} = $bits[1];
	$ref{'refstart'} = int($bits[2]); # it is 0-based
	$ref{'refalignlen'} = int($bits[3]);
	$ref{'refstrand'} = $bits[4];
	$ref{'refsize'} = int($bits[5]);
	$ref{'refseq'} = $bits[6];
	#my $sbases = $bits[6]; $sbases =~ s/\-//g; $ref{'refend'} = $ref{'refstart'} + length($sbases); # start is 0-based
	$ref{'refend'} = $ref{'refstart'} + $ref{'refalignlen'}; # start is 0-based
	$alignRef->{'ref'} = \%ref;
	
	# work on ^s line
	chomp($qline);
	push @{$alignRef->{lines}}, $qline if (0!=$recordLine);
	@bits = split(/\s+/, $qline);
	my %read=();
	$read{'read'} = $bits[1];
	$read{'readstart'} = int($bits[2]); # it is 0-based
	$read{'readalignlen'} = int($bits[3]);
	$read{'readstrand'} = $bits[4];
	$read{'readsize'} = int($bits[5]);
	$read{'readseq'} = $bits[6];
	#my $qbases = $bits[6]; $qbases =~ s/\-//g; $read{'readend'} = $read{'readstart'} + length($qbases); # start is 0-based
	$read{'readend'} = $read{'readstart'} + $read{'readalignlen'}; # start is 0-based
	$alignRef->{'read'} = \%read;
	
	#####
	# we need to make read as the reference '+' coordinates for easier assembly later
	if ('-' eq $read{'readstrand'}) {
		$ref{'refstrand'} = '-';
		$read{'readstrand'} = '+';
		
		# from .maf file
		# Coordinates are 0-based.  For - strand matches, coordinates
		# in the reverse complement of the 2nd sequence are used.
		$read{'readstart'} = $read{'readsize'} - $read{'readstart'} - $read{'readalignlen'}; # start is 0-based
		$read{'readend'} = $read{'readstart'} + $read{'readalignlen'}; # +1 'cos start is 0-based
		
		# NOTE: we do not reverse complement the sequences as the profile computation is the same
		#       however, the CIGAR generation should be taken care of
		$ref{'refseq'} = reverse $ref{'refseq'};
		$read{'readseq'} = reverse $read{'readseq'};
	}
	
	my %profile = ();
	_profileTwoSeqs($ref{'refseq'}, $read{'readseq'}, \%profile);
	$alignRef->{'profile'} = \%profile;
	
	# generate the CIGAR string
	my @cigar = ();
	generateCIGAR($ref{'refseq'}, $read{'readseq'}, $read{'readstart'}, $read{'readsize'}, $read{'readstrand'}, \@cigar);
	my $cigar = getCIGARFromBits(\@cigar, '-' eq $ref{refstrand} ?1:0);
	$read{'cigar'} = $cigar;
	
	# work on ^q line
	chomp($qualityLine);
	push @{$alignRef->{lines}}, $qualityLine if (0!=$recordLine);
	
	# clear sequences to reduce memory consumption and easier debugging
	delete $ref{'refseq'};
	delete $read{'readseq'};
}

#####

sub computeIntervalCoverage {
    my ($x1s, $x1e, $x2s, $x2e) = @_;

    if ($x1e < $x2s) {
        # x1 is smaller than x2
        return (0, $x1e-$x1s+1, $x2e-$x2s+1)
    } elsif ($x2e < $x1s) {
        # x2 is smaller than x1
        return (0, $x1e-$x1s+1, $x2e-$x2s+1)
    } else {
        # possible overlap
        my @xs = ($x1s, $x1e, $x2s, $x2e);
        @xs = sort {$a<=>$b} @xs;
        my $interval = $xs[3] - $xs[0] + 1;
        my $x1 = ($x1e-$x1s+1);
        my $x2 = ($x2e-$x2s+1);
        my $overlap = $x1 + $x2 - $interval;
        return ($overlap, $x1, $x2);
    }
}

#####

sub loadReadFastqFile{
	my ($file, $seqsRef) = @_;
	open FQFILE, $file || die "Fail to open $file\n$!\n";
	while (<FQFILE>) {
		chomp();
		$_ =~ s/^@//;
		my @bits = split(/\s+/);
		my $id = $bits[0];
		die "Read id $id already exists.\n" if (exists $seqsRef->{$id});
		my $seq = <FQFILE>; chomp($seq); $seqsRef->{$id} = {seq=>$seq};
		my $ignore = <FQFILE>;
		my $quality = <FQFILE>; chomp($quality); $seqsRef->{$id}->{qual} = $quality;
	}
	close FQFILE;
}

#####

sub loadGenomeFastaFile {
	my ($file, $seqsRef) = @_;
	
	print STDERR "Loading $file..\n";
	
	%{$seqsRef} = ();
	open GFILE, $file || die "Fail to open $file\n$!\n";
	my $id = undef;
	while (<GFILE>) {
		chomp();
		if (/^>/) {
			if (defined $id) { print STDERR " ", length($seqsRef->{$id})," read.\n"; }
			my $line = $_;
			$line =~ s/^>//;
			my @bits = split(/\s+/, $line);
			$id = lc($bits[0]);
			die "Duplicate entry $id!\n" if (exists $seqsRef->{$id});
			$seqsRef->{$id} = '';
			print STDERR "Loading ", $id, "..";
		} else {
			$seqsRef->{$id} .= $_;
		}
	}
	if (defined $id) { print STDERR " ", length($seqsRef->{$id})," read.\n"; }
	close GFILE;
	print STDERR "$file loaded.\n";
}

#####

our $G_BASECALLER_KMER=5;
# $cutoffsRef->{'_del_min_sdiff'} = $G_DEL_MIN_SDIFF
# $cutoffsRef->{'_del_max_qdiff'} = $G_DEL_MAX_QDIFF
# $cutoffsRef->{'_inv_max_sdiff'} = $G_INTRA_TRANSLOCATION_MIN_SDIFF
sub isAffectedByHomopolymer {
	my ($readId, $readSeq, $readQual, $readSize, $alignmentIRef, $alignmentJRef, $genomeSequencesRef, $cutoffsRef) = @_;
	
	if ('+' eq $alignmentIRef->{refStrand}) {
		if ($alignmentIRef->{refEnd}<$alignmentJRef->{refStart}) {
			my $span = $alignmentJRef->{refStart} - $alignmentIRef->{refEnd}; #0-based; remove + 1
			if ($span>$cutoffsRef->{'_del_min_sdiff'}) {
				my $readGap = $alignmentJRef->{qStart} - $alignmentIRef->{qEnd}; #0-based; remove + 1
				if ($readGap<$cutoffsRef->{'_del_max_qdiff'}) {
					if ($readGap<=$cutoffsRef->{'_inv_max_sdiff'}) {
						my $leftFlankStart = $alignmentIRef->{qEnd} - $G_BASECALLER_KMER; # TODO: need - 1?
						my $leftFlankSeq = uc(substr($readSeq, $leftFlankStart, $G_BASECALLER_KMER));
						my $leftNT = substr($leftFlankSeq, -1, 1);
						
						my $rightFlankStart = $alignmentJRef->{qStart};
						my $rightFlankSeq = uc(substr($readSeq, $rightFlankStart, $G_BASECALLER_KMER));
						my $rightNT = substr($rightFlankSeq, 0, 1);
						
						my $deleteLen = $span;
						my $deleteSequence = substr($genomeSequencesRef->{lc($alignmentIRef->{refId})}, $alignmentIRef->{refEnd}, $span); # check, no need -1
						
						# to determine the possible poly-5-mer
						my $leftCheckRef = undef; my $rightCheckRef = undef;
						if ($leftFlankSeq eq ($leftNT x $G_BASECALLER_KMER)) {
							$leftCheckRef = {where=>'L', type=>'F', polymer=>$leftFlankSeq, nt=>$leftNT};
							if ($rightFlankSeq eq ($rightNT x $G_BASECALLER_KMER)) {
								$rightCheckRef = {where=>'R', type=>'F', polymer=>$rightFlankSeq, nt=>$rightNT};
							}
						} else {
							if ($rightFlankSeq eq ($rightNT x $G_BASECALLER_KMER)) {
								$rightCheckRef = {where=>'R', type=>'F', polymer=>$rightFlankSeq, nt=>$rightNT};
							}
						}
						# we may have to attempt to look for split polymer caused by suboptimal alignment
						if (!defined $leftCheckRef) {
							my $leftFlankSeqSkip1 = uc(substr($readSeq, $leftFlankStart+1, $G_BASECALLER_KMER));
							if ($leftFlankSeqSkip1 eq ($leftNT x $G_BASECALLER_KMER)) {
								$leftCheckRef = {where=>'L', type=>'S', polymer=>$leftFlankSeqSkip1, nt=>$leftNT};
							}
						}
						if (!defined $rightCheckRef) {
							my $rightFlankSeqSkip1 = uc(substr($readSeq, $rightFlankStart-1, $G_BASECALLER_KMER));
							if ($rightFlankSeqSkip1 eq ($rightNT x $G_BASECALLER_KMER)) {
								$rightCheckRef = {where=>'R', type=>'S', polymer=>$rightFlankSeqSkip1, nt=>$rightNT};
							}
						}
						
						my %results = ();
						_simulatePolymer($deleteSequence, $leftCheckRef, $rightCheckRef, $leftFlankSeq, $rightFlankSeq, \%results);
						if ('' eq $results{deletedSeq}) {
							# detect homopolymer
							print STDERR 'RMV', "\t", 'HOMOPOLYMER', "\t", '+', "\t", $readId, "\t", defined $leftCheckRef ? $leftCheckRef->{type}.'-'.$leftCheckRef->{polymer} : $leftFlankSeq, "\t", $deleteSequence, "\t", defined $rightCheckRef ? $rightCheckRef->{type}.'-'.$rightCheckRef->{polymer} : $rightFlankSeq, "\n";
							return 1;
						} else {
							if ($results{leftPoly}>0 || $results{rightPoly}>0) {
								print STDERR (length($results{deletedSeq})<=$cutoffsRef->{'_del_min_sdiff'})?"RMV":"STY", "\t", 'PARTIAL', "\t", '+', "\t", $readId, "\t", length($deleteSequence), "\t", length($results{deletedSeq}), "\t", defined $leftCheckRef ? $leftCheckRef->{type}.'-'.$leftCheckRef->{polymer} : $leftFlankSeq, "\t", $results{deletedSeq}, "\t", defined $rightCheckRef ? $rightCheckRef->{type}.'-'.$rightCheckRef->{polymer} : $rightFlankSeq, "\n";
								if (length($results{deletedSeq})<=$cutoffsRef->{'_del_min_sdiff'}) {
									# if after adjusting for homopolymer the deleted sequence is too short, it is affected by homopolymer
									return 1;
								}
							}
						}
					}
				}
			}
		}
	} else {
		if ($alignmentJRef->{refEnd}<$alignmentIRef->{refStart}) {
			my $span = $alignmentIRef->{refStart} - $alignmentJRef->{refEnd}; #0-based; remove + 1
			if ($span>$cutoffsRef->{'_del_min_sdiff'}) {
				my $readGap = $alignmentJRef->{qStart} - $alignmentIRef->{qEnd}; #0-based; remove + 1
				if ($readGap<$cutoffsRef->{'_del_max_qdiff'}) {
					if ($readGap<=$cutoffsRef->{'_inv_max_sdiff'}) {
						# TODO: decide which is left and which is right, and revcomp?
						my $leftFlankStart = $alignmentJRef->{qStart}; # TODO: needs -1?
						my $leftFlankSeq = uc(substr($readSeq, $leftFlankStart, $G_BASECALLER_KMER));
						# have to reverse complement 'cos
						$leftFlankSeq =~ tr/ACGT/TGCA/; $leftFlankSeq = reverse $leftFlankSeq;
						my $leftNT = substr($leftFlankSeq, -1, 1);
						
						my $rightFlankStart = $alignmentIRef->{qEnd} - $G_BASECALLER_KMER;
						my $rightFlankSeq = uc(substr($readSeq, $rightFlankStart, $G_BASECALLER_KMER));
						$rightFlankSeq =~ tr/ACGT/TGCA/; $rightFlankSeq = reverse $rightFlankSeq;
						my $rightNT = substr($rightFlankSeq, 0, 1);
						
						my $deleteLen = $span;
						my $deleteSequence = substr($genomeSequencesRef->{lc($alignmentIRef->{refId})}, $alignmentJRef->{refEnd}, $span); # check, no need -1
						
						# to determine the possible poly-5-mer
						my $leftCheckRef = undef; my $rightCheckRef = undef;
						if ($leftFlankSeq eq ($leftNT x $G_BASECALLER_KMER)) {
							$leftCheckRef = {where=>'L', type=>'F', polymer=>$leftFlankSeq, nt=>$leftNT};
							if ($rightFlankSeq eq ($rightNT x $G_BASECALLER_KMER)) {
								$rightCheckRef = {where=>'R', type=>'F', polymer=>$rightFlankSeq, nt=>$rightNT};
							}
						} else {
							if ($rightFlankSeq eq ($rightNT x $G_BASECALLER_KMER)) {
								$rightCheckRef = {where=>'R', type=>'F', polymer=>$rightFlankSeq, nt=>$rightNT};
							}
						}
						# we may have to attempt to look for split polymer caused by suboptimal alignment
						if (!defined $leftCheckRef) {
							my $leftFlankSeqSkip1 = uc(substr($readSeq, $leftFlankStart+1, $G_BASECALLER_KMER));
							if ($leftFlankSeqSkip1 eq ($leftNT x $G_BASECALLER_KMER)) {
								$leftCheckRef = {where=>'L', type=>'S', polymer=>$leftFlankSeqSkip1, nt=>$leftNT};
							}
						}
						if (!defined $rightCheckRef) {
							my $rightFlankSeqSkip1 = uc(substr($readSeq, $rightFlankStart-1, $G_BASECALLER_KMER));
							if ($rightFlankSeqSkip1 eq ($rightNT x $G_BASECALLER_KMER)) {
								$rightCheckRef = {where=>'R', type=>'S', polymer=>$rightFlankSeqSkip1, nt=>$rightNT};
							}
						}
						
						my %results = ();
						_simulatePolymer($deleteSequence, $leftCheckRef, $rightCheckRef, $leftFlankSeq, $rightFlankSeq, \%results);
						if ('' eq $results{deletedSeq}) {
							# detect homopolymer
							print STDERR 'RMV', "\t", 'HOMOPOLYMER', "\t", '-', "\t", $readId, "\t", defined $leftCheckRef ? $leftCheckRef->{type}.'-'.$leftCheckRef->{polymer} : $leftFlankSeq, "\t", $deleteSequence, "\t", defined $rightCheckRef ? $rightCheckRef->{type}.'-'.$rightCheckRef->{polymer} : $rightFlankSeq, "\n";
							return 1;
						} else {
							if ($results{leftPoly}>0 || $results{rightPoly}>0) {
								print STDERR (length($results{deletedSeq})<=$cutoffsRef->{'_del_min_sdiff'})?"RMV":"STY", "\t", 'PARTIAL', "\t", '-', "\t", $readId, "\t", length($deleteSequence), "\t", length($results{deletedSeq}), "\t", defined $leftCheckRef ? $leftCheckRef->{type}.'-'.$leftCheckRef->{polymer} : $leftFlankSeq, "\t", $results{deletedSeq}, "\t", defined $rightCheckRef ? $rightCheckRef->{type}.'-'.$rightCheckRef->{polymer} : $rightFlankSeq, "\n";
								if (length($results{deletedSeq})<=$cutoffsRef->{'_del_min_sdiff'}) {
									# if after adjusting for homopolymer the deleted sequence is too short, it is affected by homopolymer
									return 1;
								}
							}
						}
					}
				}
			}
		}
	}
	
	return 0;
}

##### internal functions

sub _removePolymer {
	my ($sequence, $ntL, $ntR, $flankL, $flankR, $kmer, $threshold) = @_;
	
	my $seqLen = length($sequence);
	my @bits = ();
	for(my $i=0; $i<$seqLen; ++$i) {
		$bits[$i] = {nt=>substr($sequence, $i, 1), lscore=>0, rscore=>0};
	}
	
	my $poly = 0;
	my $leftPoly = 0;
	my $rightPoly = 0;
	
	# let's scan from Left
	my $lscore = $kmer;
	my $l=0;
	if (defined $ntL) {
		for($l=0; $l<$seqLen; ++$l) {
			my $bitRef = $bits[$l];
			if ($l<$kmer) {
				$lscore--;
			} else {
				if ($ntL eq $bits[$l-$kmer]->{nt}) {
					$lscore--;
				}
			}
			if ($ntL eq $bitRef->{nt}) {
				$lscore++
			}
			$bitRef->{lscore} = $lscore;
			last if ($lscore<$threshold);
		}
	}
	
	# let's scan from right
	my $rscore = $kmer;
	my $r = $seqLen;
	if (defined $ntR) {
		for(my $i=0; $i<$seqLen; ++$i) {
			$r--;
			my $bitRef = $bits[$r];
			if ($i<$kmer) {
				$rscore--;
			} else {
				if ($ntR eq $bits[$r+$kmer]->{nt}) {
					$rscore--;
				}
			}
			if ($ntR eq $bitRef->{nt}) {
				$rscore++
			}
			$bitRef->{rscore} = $rscore;
			last if ($rscore<$threshold);
		}
	}
	
	my $remaining = '';
	if ($l<=$r) {
		my $remainLen = $r-$l+1;
		$remaining = substr($sequence, $l, $remainLen);
		
		# WE remove known base caller effect
		$remaining =~ s/T{6,}/TTTTT/g;
		$remaining =~ s/A{6,}/AAAAA/g;
		
		# something remain, set up the appropriate count
		$leftPoly = $l;
		$rightPoly = $seqLen - $r;
	} else {
		# nothing left
		$leftPoly = $l;
		$rightPoly = $seqLen - $r;
	}
	$poly |= 1 if ($leftPoly>0);
	$poly |= 2 if ($rightPoly>0);
	
	return ($poly, $leftPoly, $rightPoly, $remaining);
}


sub _simulatePolymer {
	my ($deleteSequence, $leftCheckRef, $rightCheckRef, $flankL, $flankR, $resultsRef) = @_;
	$resultsRef->{poly} = $resultsRef->{leftPoly} = $resultsRef->{rightPoly} = 0;
	if (defined $leftCheckRef || defined $rightCheckRef) {
		# let's check
		my $nonHomopolymer = uc($deleteSequence);
		
		my ($poly, $leftPoly, $rightPoly, $remaining) =
		_removePolymer($nonHomopolymer,
		(defined $leftCheckRef) ? $leftCheckRef->{nt} : undef,
		(defined $rightCheckRef) ? $rightCheckRef->{nt} : undef,
		$flankL, $flankR, 5, 4);
		
		$resultsRef->{poly} = $poly;
		$resultsRef->{leftPoly} = $leftPoly;
		$resultsRef->{rightPoly} = $rightPoly;
		$resultsRef->{deletedSeq} = $remaining;
	} else {
		$resultsRef->{deletedSeq} = $deleteSequence;
	}
}

##### END - internal funcitons
1;

__END__
