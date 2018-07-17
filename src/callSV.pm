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

package callSV;

BEGIN {
	$VERSION = '0.1';
}

#####
# TODO:
# 1) abstract the sam file header
# 2) bed file writing
# 3) parameters to control SV calling
#

use strict;
use warnings;

use Data::Dumper;
use Storable qw(dclone);
use Getopt::Long;
use base 'Exporter';
our @EXPORT = qw(runCallStructuralVariants);

use utilities;
use SVCaller;
use SVTDSRCaller;
use SVTDCCaller;
use SVCTLCCaller;
use SVINVCaller;
use SVDELCaller;
use SVINSCaller;
use SVINDELCaller;
use SVTTLCCaller;

#####
my %G_LabelsIndicator = ('TDSR'=>0x01, 'TDC'=>0x02, 'INV'=>0x04, 'DEL'=>0x08, 'INS'=>0x10, 'INDEL'=>0x20, 'CTLC'=>0x40, 'TTLC'=>0x80, 'LAST'=>0x100);
my @G_SVLabels = grep { 'LAST' ne $_} sort keys %G_LabelsIndicator;
my $G_NumIndicatorBits = scalar(grep { $_ ne 'LAST' } keys %G_LabelsIndicator);
my %G_IndicatorsLabels = (0=>'None'); grep { $G_IndicatorsLabels{$G_LabelsIndicator{$_}}=$_ ;} keys %G_LabelsIndicator;

my @G_ReportCols = ('score', 'EG2', 'E', '=', '%=', 'X', '%X', 'D', '%D', 'I', '%I', 'qStart', 'qEnd', 'qStrand', 'qALen', 'q%', 'refId', 'refStart', 'refEnd', 'refStrand', 'refALen');
my @G_ReportLabels = ('read', 'rlen', 'fid', 'ftotal', 'type', 'FragSV');
push @G_ReportLabels, sort { $G_LabelsIndicator{$a}<=>$G_LabelsIndicator{$b} } @G_SVLabels;
push @G_ReportLabels, 'score', 'EG2', 'E', 'Identity', 'IdentityP', 'Mismatch', 'MismatchP', 'Deletion', 'DeletionP', 'Insertion', 'InsertionP', 'qStart', 'qEnd', 'qStrand', 'qALen', 'qP', 'refId', 'refStart', 'refEnd', 'refStrand', 'refALen', 'notes', 'sdiff', 'qdiff', 'diffclass';
#####

my $G_USAGE = "
$0 callSV --in <alignFile> --fastq <fqFile> --lastpara <last parameters> [--genome <genomeFastaFile> --removehomopolymerdeletion] [--sam] [--exlucde <chromosomeToExeclude> [--exlucde <anotherChromosomeToExeclude>]]

  --oprefix STR   prefix for output files
  --fastq STR     .fastq file
  --lastpara STR  lastal parameters used
  --removehomopolymerdeletion
                  exclude DEL and INDEL possibly affected by homopolymer
  --genome STR    genome sequence in .fasta file
  --sam           flag to output .sam file
  --exclude STR   exclude SV invovling specified chromosome
                  (specify each chromosome with --exclude individually)
  --multiloci     report SVs on best alignment of multi-loci aligments
";

sub runCallStructuralVariants {
	my $oprefix = undef;
	my $fqfile = undef;
	my $genome = undef;
	my %genomeSequences = ();
	my $G_removeHomopolymerDeletion = 0;
	my $outputSam = 0;
	my $lastpara = '-r1 -q1 -a0 -b2 -Q1'; # optional parameters -v -P28
	my $G_readseq = 0;
	my @G_excludeChrs = ();
	my $G_excludeChrsRef = undef;
	my $G_multiloci = 0;
	my $help = 0;
	
	GetOptions (
	"oprefix=s"    => \$oprefix,
	"fastq=s"      => \$fqfile,
	"genome=s"     => \$genome,
	"removehomopolymerdeletion" => \$G_removeHomopolymerDeletion,
	"sam"          => \$outputSam,
	"lastpara"     => \$lastpara,
	"exclude=s"    => \@G_excludeChrs,
	"multiloci"    => \$G_multiloci,
	"help!"        => \$help)
	or die("Error in command line arguments\n$G_USAGE");
	
	die "$G_USAGE" if ($help);
	
	die "Please specify the prefix for output files\n$G_USAGE" if (!defined $oprefix || '' eq $oprefix);
	die "Please specify the fastq files\n$G_USAGE" if ((!defined $fqfile || '' eq $fqfile) && (0!=$outputSam));
	
	my %svParameters = ();
	$svParameters{'_ctlc_min_sdiff'} = 1000000000;  # exclude from this study
	$svParameters{'_inv_max_sdiff'} = 100000;
	$svParameters{'_tdc_min_span'} = 100;
	$svParameters{'_tdc_min_sqdelta'} = 100;
	$svParameters{'_del_min_sdiff'} = 20;
	$svParameters{'_del_max_qdiff'} = 20;
	if (0!=$G_removeHomopolymerDeletion) {
		$svParameters{'_remove_homopolymer_artifact'} = $G_removeHomopolymerDeletion;
		$svParameters{'_genome_sequence_hash'} = \%genomeSequences;
	} else {
		$svParameters{'_remove_homopolymer_artifact'} = 0;
		$svParameters{'_genome_sequence_hash'} = undef;
	}
	$svParameters{'_ins_min_qdiff'} = 20;
	$svParameters{'_ins_max_sdiff'} = 20;
	$svParameters{'_indel_min_qdiff'} = 20;
	$svParameters{'_indel_min_sdiff'} = 20;
	
	if (scalar(@G_excludeChrs)>0) {
		my %excludeChrs = ();
		grep { $excludeChrs{lc($_)}=1; } @G_excludeChrs;
		$G_excludeChrsRef = \%excludeChrs;
	}
	
	my $outfile = sprintf("%s.profile.xls", $oprefix);
	open OUTP, ">$outfile" || die "Fail to open $outfile\n$!\n";
	# for R, we have the first line as header
	print OUTP join("\t", @G_ReportLabels), "\n";
	
	my $excludeFile = sprintf("%s.profile.exclude", $oprefix);
	open my $fhExclude, ">$excludeFile" || die "Fail to open $excludeFile\n$!\n";
	
	# read the sequence
	my %readSequences = ();
	$G_readseq = 1 if (0!=$outputSam); # if we are writing SAM, we need the read sequence
	utilities::loadReadFastqFile($fqfile, \%readSequences) if (0!=$G_readseq);

	# load the genome sequences; need for homopolymer checking
	if ($G_removeHomopolymerDeletion) {
		die "Please specify a genome sequence fasta file for homopolymer detection\n$G_USAGE" if (!defined $genome || '' eq $genome);
		utilities::loadGenomeFastaFile($genome, \%genomeSequences);
	} else {
		if (defined $genome && '' ne $genome) {
			print STDERR "Genome sequence fasta file is needed for homopolymer detection only.\n";
		}
	}
	
	my %svCallers = ();
	my @svCallers = ();
	my $svCallerRef; my $SVFile; my $SVType;
	$SVType = 'TDSR'; $SVFile = sprintf("%s.profile.%s.xls", $oprefix, $SVType); $svCallerRef = new SVTDSRCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'TDC'; $SVFile = sprintf("%s.profile.%s.xls", $oprefix, $SVType); $svCallerRef = new SVTDCCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'CTLC'; $SVFile = sprintf("%s.profile.%s.xls", $oprefix, $SVType); $svCallerRef = new SVCTLCCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'INV'; $SVFile = sprintf("%s.profile.%s.xls", $oprefix, $SVType); $svCallerRef = new SVINVCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'DEL'; $SVFile = sprintf("%s.profile.%s.xls", $oprefix, $SVType); $svCallerRef = new SVDELCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'INS'; $SVFile = sprintf("%s.profile.%s.xls", $oprefix, $SVType); $svCallerRef = new SVINSCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'INDEL'; $SVFile = sprintf("%s.profile.%s.xls", $oprefix, $SVType); $svCallerRef = new SVINDELCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'TTLC'; $SVFile = sprintf("%s.profile.%s.xls", $oprefix, $SVType); $svCallerRef = new SVTTLCCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	my $noneHandlerRef;
	$SVType = 'NONE'; $SVFile = sprintf("%s.profile.%s.xls", $oprefix, $SVType); $noneHandlerRef = new SVCaller(\%svParameters, $SVFile);
	
	my $fhSam = undef;
	my $fhBed = undef;
	my $pg_id = 'lastal';
	my $pg_pn = 'lastal';
	my $pg_vn = '755';
	my $pg_db = 'hg19.lastdb';
	if (0!=$outputSam) {
		my $samfile = sprintf("%s.profile.sam", $oprefix);
		open $fhSam, ">$samfile" || die "Fail to open $samfile\n$!\n";
		print $fhSam "\@HD	VN:1.3\n"; # "	SO:coordinate\n"
		
		my $bedfile = sprintf("%s.profile.bed", $oprefix);
		open $fhBed, ">$bedfile" || die "Fail to open $bedfile\n$!\n";
	}
	
	
	my $SRCFILE = undef;
	my $readRef = undef;
	my $candidateId = ();
	my $candidateRef = undef;
	my $candidateLastRowRef = undef;
	my @candidates = ();
	my @colHeaders = ();
	my $numCols = 0;
	my %seenReadIds = ();
	my %statistics = ('TOTAL'=>0, 'POTENTIAL'=>0, '0-candidate'=>0, 'Single-candidate 0-fragment'=>0, 'Multi-candidates'=>0, 'Single-candidate single-fragment'=>0, 'Single-candidate multi-fragments multi-loci'=>0, 'Single-candidate multi-fragments single-locus'=>0);
	while (<STDIN>) {
		chomp();
		next if ('' eq $_);
		
		if ('# ' eq substr($_, 0, 2)) {
			## 3196  WTD01:3:2D:P:3:M:L3196  (X)
			## score EG2     E       =       %=      X       %X      D       %D      I       %I      qStart  qEnd    qStrand qALen   q%      refId   refStart        refEnd  refStrand       refALen
			
			# parse for some information
			my $afterComment = substr($_, 2);
			if ($afterComment =~ /^score\s*/) {
				# '# score '
				$_ =~ s/^\#\s+//;
				@colHeaders = split(/\t/, $_);
				$numCols = scalar(@colHeaders);
				$readRef->{cols} = dclone(\@colHeaders);
				$readRef->{numCols} = $numCols;
			} elsif ($afterComment =~ /^\d+/) {
				# '# 3196'
				### check the results
				my $line = $_; profileResults ($readRef, $G_excludeChrsRef, \%svCallers, \@svCallers, $noneHandlerRef, $outputSam, \%seenReadIds, \%statistics, $fhSam, $fhExclude, $fhBed, $G_multiloci); $_ = $line; # 'cos classifyResult perform file read operation
				
				$_ =~ s/^#\s+//;
				my @bits = split(/\t/, $_);
				$readRef = {read=>$bits[1], readsize=>$bits[0], candidates=>\@candidates};
				if (0!=$G_readseq) {
					if (exists $readSequences{$readRef->{read}}) {
						my $seqQualRef = $readSequences{$readRef->{read}};
						$readRef->{readseq} = $seqQualRef->{seq};
						$readRef->{readqual} = $seqQualRef->{qual};
					} else {
						$readRef->{readseq} = '';
						$readRef->{readqual} = '';
					}
					delete $readSequences{$readRef->{read}};
				}
				@candidates = ();
			} elsif ($afterComment =~ /^\@PG\_/) {
				# FILE=WTD09.2D.subset.v755.hg19.maf
				# @PG_ID        lastal
				# @PG_PN        lastal
				# @PG_VN        755
				# @PG_DB        hg19.lastdb
				# @PG_END
				do {
					my @pgs = split(/\s+/, $afterComment);
					if ('@PG_ID' eq $pgs[0]) {
						$pg_id = $pgs[1];
					} elsif ('@PG_PN' eq $pgs[0]) {
						$pg_pn = $pgs[1];
					} elsif ('@PG_VN' eq $pgs[0]) {
						$pg_vn = $pgs[1];
					} elsif ('@PG_DB' eq $pgs[0]) {
						$pg_db = $pgs[1];
					} else {
						print $fhExclude "Unrecognized \@PG tag\n\t", $_, "\n";
					}
					$_ = <STDIN>; chomp();
					$afterComment = '';
					$afterComment = substr($_, 2) if ('# ' eq substr($_, 0, 2));
				} while ($afterComment!~/^\@PG_END/);
				
				if (0!=$outputSam) {
					die "Please specify the lastal database used for the alignment process." if (!defined $pg_db);
					my $lastdbPrefix = $pg_db;
					$lastdbPrefix =~ s/\.lastdb$// if ($lastdbPrefix =~ /\.lastdb$/);
					my $lastdbSeqDic = $lastdbPrefix . '.seq.dict';
					die "Lastal database dictionary '$lastdbSeqDic' does not exist. $G_USAGE" if (! -f $lastdbSeqDic);
					open SEQDICT, $lastdbSeqDic || die "Fail to open $lastdbSeqDic\n$!\n";
					while (<SEQDICT>) {
						print $fhSam $_;
					}
					close SEQDICT;
					# TODO: have to allow user to specify this information as aligner may be other than last
					printf $fhSam "\@PG\tID:%s\tPN:%s\tVN:%s\tCL:lastal %s %s %s\n", $pg_id, $pg_pn, $pg_vn, $lastpara, $pg_db, $fqfile;
				}
			} else {
				print $fhExclude "Unrecognized header1:\n\t$_";
			}
			
		} elsif ('## ' eq substr($_, 0, 3)) {
			# ignore, section can be ignore
			print $fhExclude "Unrecognized header2:\n\t$_";
		} elsif ('### ' eq substr($_, 0, 4)) {
			#### candidate#1/5
			my $afterComment = substr($_, 4);
			if ('candidate#' eq substr($afterComment, 0, 10) && substr($afterComment, 10) =~ /^(\d+)\/\d+/) {
				$candidateId = $1;
				my @candidate = ();
				$candidateRef = \@candidate;
				$candidateLastRowRef = undef;
				
				push @candidates, $candidateRef;
			} elsif ($afterComment =~ /^NO candidate left/) {
				my @candidate = ();
				$candidateRef = \@candidate;
				$candidateLastRowRef = undef;
				
				push @candidates, $candidateRef;
			} else {
				print $fhExclude "Unrecognized header3:\n\t$_";
			}
		} else {
			# keep the record
			#S2214   0       0       2889    86.94   152     4.57    201     6.05    81      2.44    42      3164    +       3122    97.68   chr8    12336432        12339674        +       3242
			#=(5)
			#=S2030  0       0       2796    83.99   244     7.33    207     6.22    82      2.46    42      3164    +       3122    97.68   chr12   8428734 8431981 +       3247
			#=S1989  0       0       2775    83.61   252     7.59    197     5.94    95      2.86    42      3164    +       3122    97.68   chr11   67606513        67609737        +       3224
			#=S1972  0       0       2774    83.35   256     7.69    206     6.19    92      2.76    42      3164    +       3122    97.68   chr11   3488611 3491847 +       3236
			
			if (0==scalar(@{$readRef->{candidates}})) {
				# address the single seed case
				my @candidate = ();
				$candidateRef = \@candidate;
				$candidateLastRowRef = undef;
				
				push @candidates, $candidateRef;
			}
			
			my @bits = split(/\t/, $_);
			if ($bits[0]=~ /\-\(\d+\)/) {
				$bits[0] =~ /\-\((\d+)\)/ ;
				$candidateLastRowRef->{'#similar'} = $1;
				next;
			} elsif ($bits[0]=~ /\=\(\d+\)/) {
				$bits[0] =~ /\=\((\d+)\)/ ;
				$candidateLastRowRef->{'#samereadspan'} = $1;
				next;
			} elsif ($bits[0]=~ /\~\(\d+\)/) {
				$bits[0] =~ /\~\((\d+)\)/ ;
				$candidateLastRowRef->{'#subset'} = $1;
				next;
			}
			my %align = ();
			for(my $i=0; $i<$readRef->{numCols}; ++$i) {
				$align{$readRef->{cols}->[$i]} = $bits[$i];
			}
			if ($bits[0]=~/^\-/) {
				if (defined $candidateLastRowRef) {
					$align{$readRef->{cols}->[0]} =~ s/^\-//;
					push @{$candidateLastRowRef->{similar}}, \%align;
				} else {
					print $fhExclude "No associated previous alignment record:\n\t$_";
				}
			} elsif ($bits[0]=~/^\=/) {
				if (defined $candidateLastRowRef) {
					$align{$readRef->{cols}->[0]} =~ s/^\=//;
					push @{$candidateLastRowRef->{samereadspan}}, \%align;
				} else {
					print $fhExclude "No associated previous alignment record:\n\t$_";
				}
			} elsif ($bits[0]=~/^\~/) {
				if (defined $candidateLastRowRef) {
					$align{$readRef->{cols}->[0]} =~ s/^\~//;
					push @{$candidateLastRowRef->{subset}}, \%align;
				} else {
					print $fhExclude "No associated previous alignment record:\n\t$_";
				}
			} else {
				push @{$candidateRef}, \%align;
				$candidateLastRowRef = \%align;
				
				$align{samereadspan} = [];
				$align{similar} = [];
				$align{subset} = [];
				$align{'#samereadspan'} = 0;
				$align{'#similar'} = 0;
				$align{'#subset'} = 0;
			}
		}
	}
	# classify the candidates before proceeding
	###
	profileResults ($readRef, $G_excludeChrsRef, \%svCallers, \@svCallers, $noneHandlerRef, $outputSam, \%seenReadIds, \%statistics, $fhSam, $fhExclude, $fhBed, $G_multiloci);
	
	# release resources no longer in use
	%genomeSequences = (); # free up large memory
	%svCallers = (); @svCallers = (); $noneHandlerRef = undef;
	close OUTP;
	if (0!=$outputSam) { close $fhSam; close $fhBed; }
	
	####
	# let's report
	_summarizeSVCallMetrics(\%statistics, \%G_LabelsIndicator, $fhExclude);
	close $fhExclude;
}

sub printProfile {
    my ($fh, $readRef, $colsRef, $alignmentsRef, $readType, $marker, $indicators) = @_;

    my $numAligns = scalar(@{$alignmentsRef});
    for(my $i=0, my $j=1; $i<$numAligns; ++$i,++$j) {
        my $alignRef = $alignmentsRef->[$i];
        my @cols = ($readRef->{read}, $readRef->{readsize}, $i+1, $numAligns, $readType);
        # output FragSV!!!
        if (-1==$indicators) {
            push @cols, '.'; #FragSV
            push @cols, ('.') x $G_NumIndicatorBits;
        } elsif (0==$indicators) {
            push @cols, 'N'; #FragSV
            push @cols, ('N') x $G_NumIndicatorBits;
        } else {
            # FragSV
            if (!exists $alignRef->{SVs}) {
                push @cols, 'N';
            } else {
                push @cols, join(",", sort keys %{$alignRef->{SVs}});
            }
            foreach my $bit (0..($G_NumIndicatorBits-1)) {
                my $value = 1 << $bit;
                push @cols, (0!=($value & $indicators)) ? 'Y' : 'N';
            }
        }
        grep { push @cols, $alignRef->{$_}; } @G_ReportCols;
        push @cols, $marker; # 'notes'
        if ($i==($numAligns-1)) {
            push @cols, ('n.a.') x 3;
        } elsif ($alignmentsRef->[$i]->{refStrand} ne $alignmentsRef->[$j]->{refStrand}) {
            # different strand, cannot compute distance appropriately
            push @cols, ('n.a.') x 3;
        } elsif ($alignmentsRef->[$i]->{refId} ne $alignmentsRef->[$j]->{refId}) {
            # different chromosome, cannot compute distance appropriately
            push @cols, ('n.a.') x 3;
        } else {
            if ('+' eq $alignmentsRef->[$i]->{refStrand}) {
                # subject-diff
                my $sdiff = $alignmentsRef->[$j]->{refStart} - $alignmentsRef->[$i]->{refEnd};
                push @cols, $sdiff;
                my $sdifflabel = '=';
                $sdifflabel = '<' if ($sdiff<0);
                $sdifflabel = '>' if ($sdiff>0);
                
                # query-diff
                my $qdiff = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
                push @cols, $qdiff;
                my $qdifflabel = '=';
                $qdifflabel = '<' if ($qdiff<0);
                $qdifflabel = '>' if ($qdiff>0);
                
                # diff-class
                push @cols, $sdifflabel.$qdifflabel;
            } else {
                # subject-diff
                my $sdiff = $alignmentsRef->[$i]->{refStart} - $alignmentsRef->[$j]->{refEnd};
                push @cols, $sdiff;
                my $sdifflabel = '=';
                $sdifflabel = '<' if ($sdiff<0);
                $sdifflabel = '>' if ($sdiff>0);
                
                # query-diff
                my $qdiff = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
                push @cols, $qdiff;
                my $qdifflabel = '=';
                $qdifflabel = '<' if ($qdiff<0);
                $qdifflabel = '>' if ($qdiff>0);
                
                # diff-class
                push @cols, $sdifflabel.$qdifflabel;
            }
        }
        
        print $fh join("\t", @cols), "\n";
    }
}

sub printSAM {
    my ($fh, $readid, $readseq, $readqual, $readsize, $marker, $candidateRef, $category, $numCandidates, $fhBed) = @_;
    
    # $category = {., MC, SCSF, SCMFML, SCMFSL}
    # SAM record need to be ordered by coordinates
    my $alignmentsRef = undef;
    my $numRows = 0;
    if (defined $candidateRef) {
        $alignmentsRef = dclone($candidateRef);
        $numRows = scalar(@{$alignmentsRef});
    }
    
    # SAM record
    if (0==$numRows) {
        # record as unmapped with mapq=0
        
        my @cols = ();
        push @cols, $readid;
        push @cols, 4; #f lag
        push @cols, '*'; # rname
        push @cols, 0; # pos
        push @cols, 0; # mapq
        push @cols, '*'; # cigar
        push @cols, '*'; # rnext
        push @cols, 0; # pnext
        push @cols, $readsize; # tlen
        push @cols, ('' eq $readseq) ? '*' : $readseq; # seq
        push @cols, ('' eq $readqual) ? '*' : $readqual; # quality
        push @cols, 'zk:Z:'.$category;
        push @cols, 'CO:Z:'.$marker if ('' ne $marker);
        print $fh join("\t", @cols), "\n";
        return;
    }
    
    ### SA:Z:chr7,94283306,+,8000S7107M1664S,0,0;chr7,94289306,+,14195S1664M,0,0
    my $repSA = '';
    my $repCigar = '';
    #if ($numRows>1) {
    if ($numRows>=1) {
        my $alignRef = $alignmentsRef->[0];
        # TODO: mapq, nm
        my $mapq = 60;
        my $nm = 0;
		$repCigar = $alignRef->{'cigar'};
        my @bits = ($alignRef->{'refId'}, $alignRef->{'refStart'}+1, $alignRef->{refStrand}, $repCigar, $mapq, $nm);
        $repSA = join(',', @bits);
    }
	my $mapq = 60; # TODO: how do we determine the mapq w.r.t. EG2?
    $mapq = 0 if ('MC' eq $category || 'SCMFML' eq $category);
    for(my $i=0; $i<$numRows; ++$i) {
        my $alignRef = $alignmentsRef->[$i];
        
        my @cols = ();
        push @cols, $readid;
        #flag
        my $flag = 0;
        if ('-' eq $alignRef->{refStrand}) {
            $flag |= 0x10; # reverse
        }
        if ($i>0) {
            $flag |= 0x0800; # supplementary
        }
        push @cols, $flag;
        push @cols, $alignRef->{'refId'};
        push @cols, $alignRef->{'refStart'}+1;
        #TODO: mapq
        push @cols, $mapq;
        my $cigar = '';
        if (0==$i) {
            $cigar = $repCigar;
        } else {
			$cigar = $alignRef->{'cigar'};
        }
        push @cols, $cigar;
        #rnext
        my $rnext = '*';
        push @cols, $rnext;
        #pnext
        my $pnext = 0;
        push @cols, $pnext;
        #tlen
        my $tlen = $readsize;
        push @cols, $tlen;
        #seq
        my $seq = '*';
        if ('' ne $readseq) {
            if ('-' eq $alignRef->{refStrand}) {
                $seq = reverse $readseq;
                $seq =~ tr/ACGTNacgtn/TGCANtgcan/;
            } else {
                $seq = $readseq;
            }
        }
        push @cols, $seq;
        #quality
        my $quality = ('' eq $readqual) ? '*' : $readqual;
        push @cols, $quality;
		
        # TODO: to provide more information in SAM for these subgroups of similar alignments
        if (exists $alignRef->{'samereadspan'}) {
            my $numSubRows = scalar(@{$alignRef->{'samereadspan'}});
            $numSubRows = $alignRef->{'#samereadspan'} if ($numSubRows<$alignRef->{'#samereadspan'});
        }
        if (exists $alignRef->{'subset'}) {
            my $numSubRows = scalar(@{$alignRef->{'subset'}});
            $numSubRows = $alignRef->{'#subset'} if ($numSubRows<$alignRef->{'#subset'});
        }
        if (exists $alignRef->{'similar'}) {
            my $numSubRows = scalar(@{$alignRef->{'similar'}});
            $numSubRows = $alignRef->{'#similar'} if ($numSubRows<$alignRef->{'#similar'});
        }
        
        # need to report the evalue for making sense of results
        # %Id, q%, score, EG2
        # zi:f:%Id
        # zq:f:q%
        # zl:i:qALen
        # zs:i:score
        # ze:f:EG2
        # zt:c:[SEe] ##type
        # zn:Z:<notes> ## notes
        push @cols, 'zi:f:'.$alignRef->{'%='};
        push @cols, 'zq:f:'.$alignRef->{'q%'};
        push @cols, 'zl:i:'.$alignRef->{'qALen'};
        push @cols, 'zs:i:'.substr($alignRef->{'score'}, 1);
        push @cols, 'ze:f:'.$alignRef->{'EG2'};
        push @cols, 'zt:c:'.substr($alignRef->{'score'}, 0, 1);
        push @cols, 'zc:i:'.$numCandidates;
        push @cols, 'zk:Z:'.$category;
        push @cols, 'zn:Z:f'.($i+1).'/'.$numRows;
        
        # let's inject the "SA" tag
		if (0==0) {
			# WCH: coded per ngmlr
			if ($numRows>1) {
				my @fragments = ();
				for(my $j=0; $j<$numRows; ++$j) {
					next if ($j==$i);
					my $fragAlignRef = $alignmentsRef->[$j];
					# TODO: mapq, nm
					my $mapq = 60;
					my $nm = 0;
					my $cigar = $fragAlignRef->{'cigar'};
					my @bits = ($fragAlignRef->{'refId'}, $fragAlignRef->{'refStart'}+1, $fragAlignRef->{refStrand}, $cigar, $mapq, $nm);
					push @fragments, join(',', @bits);
				}
				push @cols, 'SA:Z:'.join(";", @fragments).';';
			}
		} else {
			# WCH: old code
			if (0==$i) {
				# combine all cigar
				if ($numRows>1) {
					my @fragments = ();
					for(my $i=1; $i<$numRows; ++$i) {
						my $fragAlignRef = $alignmentsRef->[$i];
						# TODO: mapq, nm
						my $mapq = 60;
						my $nm = 0;
						my $cigar = $cigar = $fragAlignRef->{'cigar'};
						my @bits = ($fragAlignRef->{'refId'}, $fragAlignRef->{'refStart'}+1, $fragAlignRef->{refStrand}, $cigar, $mapq, $nm);
						push @fragments, join(',', @bits);
					}
					push @cols, 'SA:Z:'.join(";", @fragments);
				}
			} else {
				# only the first
				### SA:Z:chr7,94283306,+,8000S7107M1664S,0,0;chr7,94289306,+,14195S1664M,0,0
				push @cols, 'SA:Z:'.$repSA;
			}
		}
        
        push @cols, 'CO:Z:'.$marker if ('' ne $marker);
        
        print $fh join("\t", @cols), "\n";
        
        # TODO: refactor
        # overload writing output to bed file for easier navigation
		my $mapq = 60;
		# bed format is (start, end]
		@cols = ($alignRef->{'refId'}, $alignRef->{'refStart'}, $alignRef->{'refEnd'}, $readid.'_f'.($i+1), $mapq, , $alignRef->{'refStrand'}, $marker);
		print $fhBed join("\t", @cols), "\n";
    }
}

sub profileResults {
    my ($readRef, $excludedChrsRef, $namedSVCallersRef, $orderedSVCallersRef, $noneHandlerRef, $outputSam, $seenReadIdsRef, $statisticsRef, $fhSam, $fhExclude, $fhBed, $G_multiloci) = @_;
    
    return if (!defined $readRef);
    
    if (exists $seenReadIdsRef->{$readRef->{read}}) {
        print $fhExclude "# ", $readRef->{read}, " DUPLICATED ENTRY!\n";
        return ;
    }
    $seenReadIdsRef->{$readRef->{read}}=1;
    
    $statisticsRef->{'TOTAL'}++;
    
    # let's attempt to classify the read alignments
	my $candidateType = '';
    my $numCandidates = scalar(@{$readRef->{candidates}});
    if (0==$numCandidates) {
        $statisticsRef->{'0-candidate'}++;
        print $fhExclude "# ", $readRef->{read}, " ", $numCandidates, " candidate\n";
		printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, '', undef, '.', $numCandidates, $fhBed) if (0!=$outputSam);
        return;
    }
    if ($numCandidates>1) {
        $statisticsRef->{'Multi-candidates'}++;
		$candidateType = 'MC';
		if (0==$G_multiloci) {
	        printProfile (*OUTP, $readRef, $readRef->{cols}, $readRef->{candidates}->[0], $candidateType, '', -1);
			printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, '', $readRef->{candidates}->[0], $candidateType, $numCandidates, $fhBed) if (0!=$outputSam);
        	return;
		}
    }

    # there is only a single candidate!!!
    my $numFragments = scalar(@{$readRef->{candidates}->[0]});
    if (0==$numFragments) {
        $statisticsRef->{'Single-candidate 0-fragment'}++ if ('' eq $candidateType);
        print $fhExclude "# ", $readRef->{read}, " ", 0, " candidate\n";
		printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, '', undef, $candidateType, $numCandidates, $fhBed) if (0!=$outputSam);
        return;
    }
    if (1==$numFragments) {
        # there is no SV here
		if ('' eq $candidateType) {
	        $statisticsRef->{'Single-candidate single-fragment'}++;
			$candidateType = 'SCSF';
		}
        printProfile (*OUTP, $readRef, $readRef->{cols}, $readRef->{candidates}->[0], $candidateType, '', -1);
		# mapq = 60 given it is single-fragment
		printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, '', $readRef->{candidates}->[0], $candidateType, $numCandidates, $fhBed) if (0!=$outputSam);
        return;
    }
    
    # single candidate, multiple fragments!!!
    # can be same region, or SVs
    my $multiFragments = 0;
    my @countCigar = ();
    foreach my $alignRef (@{$readRef->{candidates}->[0]}) {
        my $numSameReadSpan = scalar(@{$alignRef->{'samereadspan'}});
        $numSameReadSpan = $alignRef->{'#samereadspan'} if ($numSameReadSpan<$alignRef->{'#samereadspan'});
        my $numSubset = scalar(@{$alignRef->{'subset'}});
        $numSubset = $alignRef->{'#subset'} if ($numSubset<$alignRef->{'#subset'});
        my $numSimilar = scalar(@{$alignRef->{'similar'}});
        $numSimilar = $alignRef->{'#similar'} if ($numSimilar<$alignRef->{'#similar'});
        
        my $cigar = '1';
        if ($numSameReadSpan>0) {
            $cigar .= '='.$numSameReadSpan;
            $multiFragments = 1;
        }
        if ($numSubset>0) {
            $cigar .= '~'.$numSubset;
            $multiFragments = 1;
        }
        if ($numSimilar>0) {
            $cigar .= '-'.$numSimilar;
            $multiFragments = 1;
        }
        
        push @countCigar, $cigar;
    }

    if (0!=$multiFragments) {
        # multi-loci
		if ('' eq $candidateType) {
	        $statisticsRef->{'Single-candidate multi-fragments multi-loci'}++;
			$candidateType = 'SCMFML';
		}
		if (0==$G_multiloci) {
			# write sam file as multi-mapped (first candidate)
			printProfile (*OUTP, $readRef, $readRef->{cols}, $readRef->{candidates}->[0], $candidateType, join(",", @countCigar), -1);
			printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, '', $readRef->{candidates}->[0], $candidateType, $numCandidates, $fhBed) if (0!=$outputSam);
			return;
		}
    }

    # question: can be inversion and TD at the same time?
	if ('' eq $candidateType) {
    	$statisticsRef->{'Single-candidate multi-fragments single-locus'}++;
		$candidateType = 'SCMFSL';
	}
    my $labelIndicator = 0;
    my $marker = '';
    my @SVs = ();
    if (0==_containExcludeChrs($readRef->{candidates}->[0], $excludedChrsRef)) {
		# $namedSVCallersRef, $orderedSVCallersRef
		foreach my $svCallerRef (@{$orderedSVCallersRef}) {
			if ($svCallerRef->containsSV($readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, $readRef->{candidates}->[0])>0) {
				push @SVs, $svCallerRef->toCO();
				$labelIndicator |= $G_LabelsIndicator{$svCallerRef->getSVType()};
				
				$svCallerRef->writeSVList($readRef);
			}
		}
		
        $marker = join(";", @SVs);
        $statisticsRef->{'Label-'.$labelIndicator}++;
    }
	$noneHandlerRef->writeSVList($readRef, $readRef->{candidates}->[0]) if (0==scalar(@SVs));
	
    # write potential SV candidates sam file
    printProfile (*OUTP, $readRef, $readRef->{cols}, $readRef->{candidates}->[0], $candidateType, $marker, $labelIndicator);
	printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, $marker, $readRef->{candidates}->[0], $candidateType, $numCandidates, $fhBed) if (0!=$outputSam);
}

##### internal functions

sub _containExcludeChrs {
	my ($alignmentsRef, $excludedChrsRef) = @_;
	
	return 0 if (!defined $excludedChrsRef);
	
	my $numFrags = 0;
	foreach my $alignmentRef (@{$alignmentsRef}) {
		$numFrags++ if (exists $excludedChrsRef->{lc($alignmentRef->{refId})});
	}
	
	return $numFrags;
}

sub _summarizeSVCallMetrics {
	my ($statisticsRef, $labelsIndicatorRef, $fh) = @_;
	
	$fh = *STDOUT if (!defined $fh);
	
	print $fh "\n";
	print $fh "\n";
	
	print $fh "# SUMMARY\n";
	print $fh "# ---\n";
	printf $fh "# TOTAL\t\t%d\t100.0%%\n", $statisticsRef->{'TOTAL'};
	my $denominator = $statisticsRef->{'TOTAL'}; $denominator = 1 if (0==$denominator);
	my $numerator = undef;
	$numerator = $statisticsRef->{'0-candidate'}; printf $fh "# without good alignment\t\t%d\t%.2f%%\n", $numerator, $numerator*100.0/$denominator; #if ($numerator>0);
	$numerator = $statisticsRef->{'Single-candidate 0-fragment'}; printf $fh "# without candidate\t\t%d\t%.2f%%\n", $numerator, $numerator*100.0/$denominator;
	$numerator = $statisticsRef->{'Multi-candidates'}; printf $fh "# Multi-candidates\t\t%d\t%.2f%%\n", $numerator, $numerator*100.0/$denominator;
	$numerator = $statisticsRef->{'Single-candidate single-fragment'}; printf $fh "# Single-candidate single-fragment\t\t%d\t%.2f%%\n", $numerator, $numerator*100.0/$denominator;
	$numerator = $statisticsRef->{'Single-candidate multi-fragments multi-loci'}; printf $fh "# Single-candidate multi-fragments multi-loci\t\t%d\t%.2f%%\n", $numerator, $numerator*100.0/$denominator;
	$numerator = $statisticsRef->{'Single-candidate multi-fragments single-locus'}; printf $fh "# Single-candidate multi-fragments single-locus\t\t%d\t%.2f%%\n", $numerator, $numerator*100.0/$denominator;
	
	my %SVTotals = ();
	# compute once before reporting the collated summary followed by details
	for(my $i=0; $i<$labelsIndicatorRef->{'LAST'}; ++$i) {
		my $label = sprintf("Label-%d", $i);
		if (exists $statisticsRef->{$label}) {
			if (0==$i) {
				$SVTotals{$i} += $statisticsRef->{$label};
			} else {
				$numerator = $statisticsRef->{$label};
				my @SVs = ();
				while (my ($key, $value) = each %G_LabelsIndicator) {
					if ($value==($i & $value)) {
						push @SVs, $key;
						$SVTotals{$value} += $numerator;
					}
				}
			}
		}
	}
	
	# report collated summary
	printf $fh "# SV Types\n";
	foreach my $labelIndicator (sort { $a<=>$b } keys %G_IndicatorsLabels) {
		next if ($labelIndicator==$labelsIndicatorRef->{'LAST'});
		$numerator = (exists $SVTotals{$labelIndicator}) ? $SVTotals{$labelIndicator} : 0;
		printf $fh "# \tSV=%s\t%d\t%.2f%%\n", $G_IndicatorsLabels{$labelIndicator}, $numerator, $numerator*100.0/$denominator;
	}
	
	# report details
	printf $fh "# SV Grid Cells\n";
	for(my $i=0; $i<$labelsIndicatorRef->{'LAST'}; ++$i) {
		my $label = sprintf("Label-%d", $i);
		if (exists $statisticsRef->{$label}) {
			if (0==$i) {
				$numerator = $statisticsRef->{$label};
				printf $fh "# \tSVCell=%s\t%d\t%.2f%%\n", $G_IndicatorsLabels{$i}, $numerator, $numerator*100.0/$denominator;
			} else {
				$numerator = $statisticsRef->{$label};
				my @SVs = ();
				while (my ($key, $value) = each %G_LabelsIndicator) {
					if ($value==($i & $value)) {
						push @SVs, $key;
					}
				}
				printf $fh "# \tSVCell=%s\t%d\t%.2f%%\n", join(",", @SVs), $numerator, $numerator*100.0/$denominator;
			}
		}
	}
	
	print $fh "# ---\n";
}

##### END - internal functions



1;

__END__
