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
# Created Feb 16, 2018
# Copyright (c) 2018       Chee-Hong WONG
#                          Genome Technologies
#                          The Jackson Laboratory
#
#####

package samToAlign;

BEGIN {
	$VERSION = '0.1';
}

use strict;
use warnings;

use Data::Dumper;
use base 'Exporter';
our @EXPORT = qw(runSamToAlign);

use utilities;

#my $G_SAM_CIGAR = 0;
my $G_SAM_CIGAR = 1;
my @G_SAM_ReportCols = ('score', 'EG2', 'E', '=', '%=', 'X', '%X', 'D', '%D', 'I', '%I', 'qStart', 'qEnd', 'qStrand', 'qALen', 'q%', 'refId', 'refStart', 'refEnd', 'refStrand', 'refALen', 'cigar');

sub runSamToAlign {
    # convert .sam to .align format for Picky
	my $startTime = time;

    # we omit the headers ; no .sam output in callSV if we already have .sam file as input
    # thus, we are just processing .sam record; but need to be name sorted
    my $seenHeader = 0;
    my $sortedByName = 1;
    my $line = undef;
    while ($line = <STDIN>) {
        if ($line =~ /^\@/) {
            $seenHeader = 1;
            if ($line =~/^\@HD\t/) {
                if ($line =~/\tSO:queryname/i) {
                    $sortedByName = 1;
                } else {
                    $sortedByName = 0;
                }
            }
        } else {
            last;
        }
    }
    die "# ERROR: You sam file has no header line(s). Consider: samtools view -h\n" if (0==$seenHeader);
    die "# ERROR: You sam file is not sorted by queryname. Consider: samtools sort -n\n" if (0==$sortedByName);

    my $rowCount = 0; my $readCount = 0; my $block = 1000; my $triggerBlock = $block;
    printf STDERR "%s Processing sam-->align..", utilities::getTimeStamp();

    my $prevReadId = '';
    my @records = ();
    while ($line) {
        $rowCount++;
        chomp($line);
        my @bits = split(/\t/, $line);
        if ('' eq $prevReadId || $bits[0] ne $prevReadId) {
            if ('' ne $prevReadId) {
                # transform collection to .align format
                _writeAlign(\@records);
            }

            @records = ();
            $prevReadId = $bits[0];

            $readCount++;
            if ($readCount>=$triggerBlock) {
                $triggerBlock+=$block;

                printf STDERR " %0.3fM(%0.3fM rows)..", $readCount/1e6, $rowCount/1e6;
            }
        }
        # correct all sam records of the same read id
        push @records, $line;

        $line = <STDIN>;
    }

    # process the last buffered read
    if ('' ne $prevReadId) {
        # transform collection to .align format
        _writeAlign(\@records);

        $readCount++;
    }
    printf STDERR " %0.3fM(%0.3fM rows).. done.", $readCount/1e6, $rowCount/1e6;

	my $totalTime = time - $startTime;
	printf STDERR "%s Total run time = %d secs\n", utilities::getTimeStamp(), $totalTime;
}

sub _writeAlign {
    my ($recsRef) = @_;

    my $numOfRecords = scalar(@{$recsRef});
    return if (0==$numOfRecords);

    ## 2214  MinION3_20161013_FNFAB42260_MN20093_sequencing_run_Chip98_Genomic_R9_4_480bps_17755_ch271_read1365_strand1.fast5        {!!!}   align(10,0)     seed(0) nonseed(0)
    ## score EG2     E       =       %=      X       %X      D       %D      I       %I      qStart  qEnd    qStrand qALen   q%      refId   refStart        refEnd  refStrand       refALen cigar

    ## 1035  llssbzms2p35x_20161128_FNFAB45271_MN17073_sequencing_run_Hu_Bir_R94_1Dlig_fc5_92573_ch84_read182_strand.fast5   (1)     align(47,1)     seed(1) nonseed(0)
    ## score EG2     E       =       %=      X       %X      D       %D      I       %I      qStart  qEnd    qStrand qALen   q%      refId   refStart        refEnd  refStrand       refALen cigar
    #S55     7.3e-15 4e-20   95      85.59   4       3.60    9       8.11    3       2.70    835     937     +       102     9.86    chr8    26406193        26406301        +       108     835S6=1I2=1I15=1X3=1X6=2D16=3D1=2D2=1X6=1X14=1I12=1D8=1D4=98S

    ## 3064  MinION2_20161013_FNFAB42706_MN16454_sequencing_run_Chip99_Genomic_R9_4_480bps_44361_ch254_read57_strand1.fast5  (X)     align(36,6)     seed(0) nonseed(0)
    ## score EG2     E       =       %=      X       %X      D       %D      I       %I      qStart  qEnd    qStrand qALen   q%      refId   refStart        refEnd  refStrand       refALen cigar
    #### NO candidate left

    ## 2965  LomanLabz_PC_20161128_FNFAF01253_MN17024_sequencing_run_20161128_Human_Qiagen_1D_R9_4_50330_ch84_read1083_strand.fast5  [ ]     align(67,65)    seed(1) nonseed(11)
    ## score EG2     E       =       %=      X       %X      D       %D      I       %I      qStart  qEnd    qStrand qALen   q%      refId   refStart        refEnd  refStrand       refALen cigar
    #### candidate#1/2
    # :
    #### candidate#2/2
    # :

    # there is no selection here as the aligner has taken care of it
    # thus, there is either a representative alignment or there isn't

    my @bits = split(/\t/, $recsRef->[0]);
    return if (4==($bits[1]&4)); # unmapped

    # print header
    my $readLength = 0; # NOT used, but informative for checking
    $readLength = _computeReadLength ($bits[5]);
    my $readId = $bits[0];
    print '# ', $readLength, "\t", $readId, "\n"; # we omit the optional 3rd-6th cells

    # print columns header
    print '# ', join("\t", @G_SAM_ReportCols), "\n";

    # print alignemtns
    # TODO: just write a lot of '.' for values and see what are mandatory
    # [ ] score (opt)
    # [ ] EG2 (opt)
    # [ ] E (opt)
    # [ ] = (opt)
    # [ ] %= (opt)
    # [ ] X (opt)
    # [ ] %X (opt)
    # [ ] D (opt)
    # [ ] %D (opt)
    # [ ] I (opt)
    # [ ] %I (opt)
    # [X] qStart (TODO: need to compute)
    # [X] qEnd (TODO: need to compute)
    # [X] qStrand (fixed:+)
    # [ ] qALen (opt; sam)
    # [ ] q% (opt)
    # [X] refId (col2)
    # [X] refStart (col3)
    # [X] refEnd (TODO: need to compute)
    # [X] refStrand (col1) 
    # [ ] refALen (opt; sam)
    # [ ] cigar (opt; sam)
    if (1==$numOfRecords) {
        # no candidate block for single-fragment
        foreach my $line (@{$recsRef}) {
            @bits = split(/\t/, $line);
            my %item = ();
            grep { $item{$_} = '.'; } @G_SAM_ReportCols;

            # fill up the mandatory information
            $item{qStrand} = '+';
            $item{refId} = $bits[2];
            $item{refStrand} = (16==(16&$bits[1])) ? '-' : '+';
            my $qStart = 0;
            my $qEnd = 0;
            my $refStart = int($bits[3]); #POS is 1-based
            my $refEnd = 0;
            ($qStart, $qEnd, $refStart, $refEnd) = _computeCoordinates($item{refStrand}, $bits[5], $refStart);
            $item{qStart} = $qStart;
            $item{qEnd} = $qEnd;
            $item{refStart} = $refStart;
            $item{refEnd} = $refEnd;

            $item{cigar} = $bits[5] if (0!=$G_SAM_CIGAR);

            # report
            my @cols = ();
            grep { push @cols, $item{$_}; } @G_SAM_ReportCols;
            print join("\t", @cols), "\n";
            print "\n";
        }
    } else {
        print '### candidate#1/1', "\n";
        # TODO:
        foreach my $line (@{$recsRef}) {
            @bits = split(/\t/, $line);
            my %item = ();
            grep { $item{$_} = '.'; } @G_SAM_ReportCols;

            # fill up the mandatory information
            $item{qStrand} = '+';
            $item{refId} = $bits[2];
            $item{refStrand} = (16==(16&$bits[1])) ? '-' : '+';
            my $qStart = 0;
            my $qEnd = 0;
            my $refStart = int($bits[3]); #POS is 1-based
            my $refEnd = 0;
            ($qStart, $qEnd, $refStart, $refEnd) = _computeCoordinates($item{refStrand}, $bits[5], $refStart);
            $item{qStart} = $qStart;
            $item{qEnd} = $qEnd;
            $item{refStart} = $refStart;
            $item{refEnd} = $refEnd;

            $item{cigar} = $bits[5] if (0!=$G_SAM_CIGAR);

            # report
            my @cols = ();
            grep { push @cols, $item{$_}; } @G_SAM_ReportCols;
            print join("\t", @cols), "\n";
        }
        print "\n";
    }
}

sub _computeCoordinates {
    my ($refStrand, $cigar, $refStart) = @_;
    my ($cqStart, $cqEnd, $crefStart, $crefEnd) = (0, 0, $refStart-1, $refStart-1); # 0-based

    my @cigarBits = split(/([MIDNSHP=X])/, $cigar);
    my $numberOfCigarBits = scalar(@cigarBits);
    my @cigarOps = ();
    for(my $i=1; $i<$numberOfCigarBits; $i+=2) {
        push @cigarOps, {ops=>$cigarBits[$i], bases=>$cigarBits[$i-1]};
    }
    @cigarBits = ();
    @cigarOps = reverse @cigarOps if ('-' eq $refStrand);
    my $numberOfCigarOps = scalar(@cigarOps);
    if ($numberOfCigarOps>=0) {
        my $startIdx = 0;
        if ('S' eq $cigarOps[0]->{ops} || 'H' eq $cigarOps[0]->{ops}) {
            $cqStart += int($cigarOps[0]->{bases}); # 0-based
            $cqEnd += int($cigarOps[0]->{bases}); # 0-based
            $startIdx = 1;
        }
        for(my $i=$startIdx; $i<$numberOfCigarOps; ++$i) {
            if ('M' eq $cigarOps[$i]->{ops} || '=' eq $cigarOps[$i]->{ops} || 'X' eq $cigarOps[$i]->{ops}) {
                # increase query, and reference
                $cqEnd += int($cigarOps[$i]->{bases});
                $crefEnd += int($cigarOps[$i]->{bases});
            } elsif ('I' eq $cigarOps[$i]->{ops}) {
                # increase query
                $cqEnd += int($cigarOps[$i]->{bases});
            } elsif ('D' eq $cigarOps[$i]->{ops} || 'N' eq $cigarOps[$i]->{ops}) {
                # increase reference
                $crefEnd += int($cigarOps[$i]->{bases});
            } elsif ('S' eq $cigarOps[$i]->{ops} || 'H' eq $cigarOps[$i]->{ops}) {
                # should only happen at the end of cigar!!
                # TODO: bait if this is not the case!
                die "CIGAR op=$cigarOps[$i]->{ops} is only allowed at start and end, but found at $i\n" if ($i!=($numberOfCigarOps-1));
            }
        }
    }

    return ($cqStart, $cqEnd, $crefStart, $crefEnd); # return 0-based start, 0-based end
}

sub _computeReadLength {
    my ($cigar) = @_;

    # read length does not need to consider the reference strand
    my $readLength = 0;
    my @cigarBits = split(/([MIDNSHP=X])/, $cigar);
    my $numberOfCigarBits = scalar(@cigarBits);
    if ($numberOfCigarBits>2) {
        for(my $i=1; $i<$numberOfCigarBits; $i+=2) {
            if ('M' eq $cigarBits[$i] || '=' eq $cigarBits[$i] || 'X' eq $cigarBits[$i]) {
                $readLength += int($cigarBits[$i-1]);
            } elsif ('I' eq $cigarBits[$i]) {
                $readLength += int($cigarBits[$i-1]);
            } elsif ('D' eq $cigarBits[$i] || 'N' eq $cigarBits[$i]) {
                # do nothing
            } elsif ('S' eq $cigarBits[$i] || 'H' eq $cigarBits[$i]) {
                $readLength += int($cigarBits[$i-1]);
            }
        }
    }

    return $readLength;
}
