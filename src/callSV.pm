#####
#
# Picky - Structural Variants Pipeline for (ONT) long read
#
# Copyright (c) 2016-2017  Chee-Hong WONG  The Jackson Laboratory
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
$0 callSV --in <alignFile> --fastq <fqFile> --genome <genomeFastaFile>

  --in STR        .align file
  --fastq STR     .fastq file
  --genome STR    genome sequence in .fasta file
  --removehomopolymerdeletion
                  exclude DEL and INDEL possibly affected by homopolymer
  --sam           flag to output .sam file
  --readseq       read sequenece to be output in .sam file
  --exclude STR   exclude SV invovling specified chromosome
  --correctminus  flag to adjust for off-by-one error in .align file
";

sub runCallStructuralVariants {
	my $file = undef;
	my $fqfile = undef;
	my $outfile = undef;
	my $genome = undef;
	my %genomeSequences = ();
	my $G_removeHomopolymerDeletion = 0;
	my $outputSam = 0;
	my $G_readseq = 0;
	my @G_excludeChrs = ();
	my $G_excludeChrsRef = undef;
	my $G_correctMinus = 0; # WCH: 20170216 - fix the off by 1 error in .align output by pick-maf.pl
	
	GetOptions (
	"in=s"    => \$file,
	"fastq=s" => \$fqfile,
	"genome=s" => \$genome,
	"removehomopolymerdeletion" => \$G_removeHomopolymerDeletion,
	"sam"     => \$outputSam,
	"readseq"     => \$G_readseq,
	"correctminus" => \$G_correctMinus,
	"exclude=s" => \@G_excludeChrs)
	or die("Error in command line arguments\n$G_USAGE");
	
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
	
	open INFILE, "$file" || die "Fail to open $file\n$!\n";
	
	$outfile = $file; $outfile =~ s/.align$// if ($outfile =~ /.align$/);
	
	my $prefixFile = $outfile;
	$outfile .= '.profile.xls';
	open OUTP, ">$outfile" || die "Fail to open $outfile\n$!\n";
	# for R, we have the first line as header
	print OUTP join("\t", @G_ReportLabels), "\n";
	
	my $excludeFile = $prefixFile . '.profile.exclude';
	open my $fhExclude, ">$excludeFile" || die "Fail to open $excludeFile\n$!\n";
	
	if (!defined $fqfile) {
		my @bits = split(/\./,$file);
		if (scalar(@bits)>2) {
			$fqfile = sprintf("%s.%s.fastq", $bits[0], $bits[1]);
		} elsif (scalar(@bits)>1) {
			$fqfile = sprintf("%s.fastq", $bits[0]);
		} else {
			$fqfile = $file . '.fastq';
		}
	}
	
	# read the sequence
	my %readSequences = ();
	utilities::loadReadFastqFile($fqfile, \%readSequences) if (0!=$G_readseq);

	# load the genome sequences; need for homopolymer checking
	utilities::loadGenomeFastaFile($genome, \%genomeSequences) if (defined $genome && '' ne $genome);
	
	my %svCallers = ();
	my @svCallers = ();
	my $svCallerRef; my $SVFile; my $SVType;
	$SVType = 'TDSR'; $SVFile = sprintf("%s.profile.%s.xls", $prefixFile, $SVType); $svCallerRef = new SVTDSRCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'TDC'; $SVFile = sprintf("%s.profile.%s.xls", $prefixFile, $SVType); $svCallerRef = new SVTDCCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'CTLC'; $SVFile = sprintf("%s.profile.%s.xls", $prefixFile, $SVType); $svCallerRef = new SVCTLCCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'INV'; $SVFile = sprintf("%s.profile.%s.xls", $prefixFile, $SVType); $svCallerRef = new SVINVCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'DEL'; $SVFile = sprintf("%s.profile.%s.xls", $prefixFile, $SVType); $svCallerRef = new SVDELCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'INS'; $SVFile = sprintf("%s.profile.%s.xls", $prefixFile, $SVType); $svCallerRef = new SVINSCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'INDEL'; $SVFile = sprintf("%s.profile.%s.xls", $prefixFile, $SVType); $svCallerRef = new SVINDELCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	$SVType = 'TTLC'; $SVFile = sprintf("%s.profile.%s.xls", $prefixFile, $SVType); $svCallerRef = new SVTTLCCaller(\%svParameters, $SVFile);
	push @svCallers, $svCallerRef; $svCallers{$SVType} = $svCallerRef;
	my $noneHandlerRef;
	$SVType = 'NONE'; $SVFile = sprintf("%s.profile.%s.xls", $prefixFile, $SVType); $noneHandlerRef = new SVCaller(\%svParameters, $SVFile);
	
	if (0!=$outputSam) {
		my $samfile = $file;
		$samfile .= '.profile.sam';
		open OUTS, ">$samfile" || die "Fail to open $samfile\n$!\n";
		# TODO: print the same file header!!!
		# print OUTS "\@HD	VN:1.3	SO:coordinate\n";
		print OUTS "\@HD	VN:1.3\n";
		print OUTS "\@SQ	SN:chr1	LN:249250621\n";
		print OUTS "\@SQ	SN:chr2	LN:243199373\n";
		print OUTS "\@SQ	SN:chr3	LN:198022430\n";
		print OUTS "\@SQ	SN:chr4	LN:191154276\n";
		print OUTS "\@SQ	SN:chr5	LN:180915260\n";
		print OUTS "\@SQ	SN:chr6	LN:171115067\n";
		print OUTS "\@SQ	SN:chr7	LN:159138663\n";
		print OUTS "\@SQ	SN:chr8	LN:146364022\n";
		print OUTS "\@SQ	SN:chr9	LN:141213431\n";
		print OUTS "\@SQ	SN:chr10	LN:135534747\n";
		print OUTS "\@SQ	SN:chr11	LN:135006516\n";
		print OUTS "\@SQ	SN:chr12	LN:133851895\n";
		print OUTS "\@SQ	SN:chr13	LN:115169878\n";
		print OUTS "\@SQ	SN:chr14	LN:107349540\n";
		print OUTS "\@SQ	SN:chr15	LN:102531392\n";
		print OUTS "\@SQ	SN:chr16	LN:90354753\n";
		print OUTS "\@SQ	SN:chr17	LN:81195210\n";
		print OUTS "\@SQ	SN:chr18	LN:78077248\n";
		print OUTS "\@SQ	SN:chr19	LN:59128983\n";
		print OUTS "\@SQ	SN:chr20	LN:63025520\n";
		print OUTS "\@SQ	SN:chr21	LN:48129895\n";
		print OUTS "\@SQ	SN:chr22	LN:51304566\n";
		print OUTS "\@SQ	SN:chrX	LN:155270560\n";
		print OUTS "\@SQ	SN:chrY	LN:59373566\n";
		print OUTS "\@SQ	SN:chrM	LN:16571\n";
		# TODO: we need to report the program info correctly!
		print OUTS "\@PG	ID:lastal	PN:lastal	VN:755	CL:last-755/src/lastal -r1 -q1 -a0 -b2 -v -P 28 -Q 1 hg19.lastdb ",$fqfile,"\n";
		
		my $bedfile = $file;
		$bedfile .= '.profile.bed';
		open OUTB, ">$bedfile" || die "Fail to open $bedfile\n$!\n";
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
	while (<INFILE>) {
		chomp();
		next if ('' eq $_);
		
		if ('# ' eq substr($_, 0, 2)) {
			## 3196  WTD01:3:2D:P:3:M:L3196  (X)
			## score EG2     E       =       %=      X       %X      D       %D      I       %I      qStart  qEnd    qStrand qALen   q%      refId   refStart        refEnd  refStrand       refALen
			
			# TODO: parse for some information
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
				my $line = $_; profileResults ($readRef, $G_excludeChrsRef, \%svCallers, \@svCallers, $noneHandlerRef, $outputSam, \%seenReadIds, \%statistics, *OUTS, $fhExclude); $_ = $line; # 'cos classifyResult perform file read operation
				
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
			} elsif ($afterComment =~ /^FILE=/) {
				my $newSrcFile = $';
				if (!defined $SRCFILE) {
					$SRCFILE = $newSrcFile;
					open INSRC, "$SRCFILE" || die "Fail to open $SRCFILE\n$!\n";
				} else {
					if ($SRCFILE ne $newSrcFile) {
						$SRCFILE = $newSrcFile;
						open INSRC, "$SRCFILE" || die "Fail to open $SRCFILE\n$!\n";
					}
				}
				print $fhExclude $_, "\n";
			} else {
				print $fhExclude "Unrecognized header1:\n\t$_";
			}
			
		} elsif ('## ' eq substr($_, 0, 3)) {
			# TODO: ignore, section can be ignore
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
			if (0!=$G_correctMinus) {
				# need to correct the off-by-1 problem
				$align{qStart}++ if ('-' eq $align{refStrand});
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
	profileResults ($readRef, $G_excludeChrsRef, \%svCallers, \@svCallers, $noneHandlerRef, $outputSam, \%seenReadIds, \%statistics, *OUTS, $fhExclude);
	
	# release resources no longer in use
	%genomeSequences = (); # free up large memory
	close INFILE;
	%svCallers = (); @svCallers = (); $noneHandlerRef = undef;
	close INSRC if (defined $SRCFILE);
	close OUTP;
	if (0!=$outputSam) { close OUTS; close OUTB; }
	
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
    my ($fh, $readid, $readseq, $readqual, $readsize, $marker, $candidateRef, $category, $numCandidates) = @_;
    
    # $category = {., MC, SCSF, SCMFML, SCMFSL}
    # SAM record need to be orderedby coordinates
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
    #my $KEY_STRAND = 'qStrand';
    my $KEY_STRAND = 'refStrand';
    my $repSA = '';
    my $repStrand = '';
    my $repCigar = '';
    #if ($numRows>1) {
    if ($numRows>=1) {
        my $alignRef = $alignmentsRef->[0];
        # TODO: mapq, nm
        my $mapq = 60;
        my $nm = 0;
        $repStrand = $alignRef->{$KEY_STRAND};
        #$repCigar = getCIGARFromBits($alignRef->{'readcigar'}, '-' eq $repStrand ?1:0);
        $repCigar = getCIGARFromBits($alignRef->{'readcigar'}, '-' eq $repStrand ?1:0);
        my @bits = ($alignRef->{'refId'}, $alignRef->{'refStart'}+1, $repStrand, $repCigar, $mapq, $nm);
        $repSA = join(',', @bits);
    }
    my $mapq = 60;
    $mapq = 0 if ('MC' eq $category || 'SCMFML' eq $category);
    for(my $i=0; $i<$numRows; ++$i) {
        my $alignRef = $alignmentsRef->[$i];
        
        my @cols = ();
        push @cols, $readid;
        #TODO: flag
        my $flag = 0;
        if ('-' eq $alignRef->{$KEY_STRAND}) {
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
            #$cigar = getCIGARFromBits($alignRef->{'readcigar'}, '-' eq $alignRef->{'qStrand'}?1:0);
            $cigar = getCIGARFromBits($alignRef->{'readcigar'}, '-' eq $alignRef->{$KEY_STRAND}?1:0);
        }
        push @cols, $cigar;
        #TODO: rnext
        my $rnext = '*';
        push @cols, $rnext;
        #TODO: pnext
        my $pnext = 0;
        push @cols, $pnext;
        #TODO: tlen
        my $tlen = $readsize;
        push @cols, $tlen;
        #TODO: seq
        my $seq = '*';
        if ('' ne $readseq) {
            if ('-' eq $alignRef->{$KEY_STRAND}) {
                ### TODO: check if we need to reverse complement for '-v' strand mapping
                $seq = reverse $readseq;
                $seq =~ tr/ACGTNacgtn/TGCANtgcan/;
            } else {
                $seq = $readseq;
            }
        }
        push @cols, $seq;
        #TODO: quality
        my $quality = ('' eq $readqual) ? '*' : $readqual;
        push @cols, $quality;
        
        
        # TODO: WCH: FIXME:
        # have to subgroup those similar alignments
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
        
        # TODO: need to report the evalue for making sense of results
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
        if (0==$i) {
            # combine all cigar
            if ($numRows>1) {
                my @fragments = ();
                for(my $i=1; $i<$numRows; ++$i) {
                    my $fragAlignRef = $alignmentsRef->[$i];
                    # TODO: mapq, nm
                    my $mapq = 60;
                    my $nm = 0;
                    #my $cigar = getCIGARFromBits($fragAlignRef->{'readcigar'}, '-' eq $fragAlignRef->{'qStrand'}?1:0);
                    my $cigar = getCIGARFromBits($fragAlignRef->{'readcigar'}, '-' eq $fragAlignRef->{$KEY_STRAND}?1:0);
                    my @bits = ($fragAlignRef->{'refId'}, $fragAlignRef->{'refStart'}+1, $fragAlignRef->{$KEY_STRAND}, $cigar, $mapq, $nm);
                    push @fragments, join(',', @bits);
                }
                push @cols, 'SA:Z:'.join(";", @fragments);
            }
        } else {
            # only the first
            ### SA:Z:chr7,94283306,+,8000S7107M1664S,0,0;chr7,94289306,+,14195S1664M,0,0
            push @cols, 'SA:Z:'.$repSA;
        }
        
        push @cols, 'CO:Z:'.$marker if ('' ne $marker);
        
        print $fh join("\t", @cols), "\n";
        
        # TODO: refactor
        # overload writing output to bed file for easier navigation
		my $mapq = 60;
		# bed format is (start, end]
		@cols = ($alignRef->{'refId'}, $alignRef->{'refStart'}, $alignRef->{'refEnd'}, $readid.'_f'.($i+1), $mapq, , $alignRef->{'refStrand'}, $marker);
		print OUTB join("\t", @cols), "\n";
    }
}

sub parseMAFRecordForCigar {
    my ($aline, $sline, $qline, $qualityLine, $alignRef) = @_;
    
    # let's process the block
    # a score=1897 EG2=0 E=0
    # s chr6                        160881684 3801 + 171115067 CCCAAGAAAAC
    # s WTD01:50183:2D:P:3:M:L13738        32 3791 +     13738 CCCAAGAAAAT
    # q WTD01:50183:2D:P:3:M:L13738                            4/0825;7<=3
    
    %{$alignRef} = ();
    
    # work on ^a line
    chomp($aline);
    my ($line) = $aline =~ /^a\s+(.*)/;
    foreach my $keyValue (split(/\s+/, $line)) {
        my ($key, $value) = split(/\=/, $keyValue);
        $alignRef->{$key} = $value;
    }
    
    # work on ^s line
    chomp($sline);
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
        
        $read{'readstart'} = $read{'readsize'} - ($read{'readstart'} + $read{'readalignlen'}); # start is 0-based
        $read{'readend'} = $read{'readstart'} + $read{'readalignlen'}; # start is 0-based
        
        # NOTE: we do not reverse complement the sequences as the profile computation is the same
        #       however, the CIGAR generation should be taken care of
        $ref{'refseq'} = reverse $ref{'refseq'};
        $read{'readseq'} = reverse $read{'readseq'};
    }
    
    # TODO: generate the CIGAR string
    # TODO: confirm the recording in MAF format for -ve strand!!!
    my @cigar = ();
    generateCIGAR($ref{'refseq'}, $read{'readseq'}, $read{'readstart'}, $read{'readsize'}, $read{'readstrand'}, \@cigar);
    $read{'readcigar'} = \@cigar;
    
    # work on ^q line
    chomp($qualityLine);
    @bits = split(/\s+/, $qualityLine);
    $read{'readquality'} = $bits[2];
    
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

sub ReadMAFRecord {
    my ($alignmentsRef) = @_;
    
    # work on a single record
    # TODO: logic to start from the beginning of the file if necessary
    my $wanted = $alignmentsRef->{'line#'};
    while (<INSRC>) {
        my $linenumber = $.;
        if ($linenumber==$wanted) {
            # let's process the record!!!
            my $aline = $_;
            my $sline = <INSRC>;
            my $qline = <INSRC>;
            my $qualityLine = <INSRC>;
            
			my $queryId = utilities::getQueryIdValue ($qline);
            my %alignment = ();
            parseMAFRecordForCigar ($aline, $sline, $qline, $qualityLine, \%alignment);
        
            $alignmentsRef->{refseq} = $alignment{ref}->{refseq};
            $alignmentsRef->{readseq} = $alignment{read}->{readseq};
            $alignmentsRef->{readcigar} = $alignment{read}->{readcigar};
            $alignmentsRef->{readquality} = $alignment{read}->{readquality};
            
            last;
        } elsif ($linenumber>$wanted) {
            # past what we want!
            die "Wanted line#", $wanted, " but already at line#", $linenumber, "!\n!\t$_\n", Dumper($alignmentsRef);
        }
    }
}

sub LoadMAFRecord {
    my ($alignmentsRef) = @_;
    
    # order the line number of the records
    # read each of them the sequences
    # compute the CIGAR base on that
    my @lineOrders = ();
    for(my $i=0; $i<scalar(@{$alignmentsRef}); ++$i) {
        push @lineOrders, {order=>$i, ref=>$alignmentsRef->[$i]};
    }
    @lineOrders = sort { int($a->{ref}->{'line#'}) <=> int($b->{ref}->{'line#'})} @lineOrders;
    
    foreach my $lineOrderRef (@lineOrders) {
        ReadMAFRecord($lineOrderRef->{ref});
    }
}

sub profileResults {
    my ($readRef, $excludedChrsRef, $namedSVCallersRef, $orderedSVCallersRef, $noneHandlerRef, $outputSam, $seenReadIdsRef, $statisticsRef, $fhSam, $fhExclude) = @_;
    
    return if (!defined $readRef);
    
    if (exists $seenReadIdsRef->{$readRef->{read}}) {
        print $fhExclude "# ", $readRef->{read}, " DUPLICATED ENTRY!\n";
        return ;
    }
    $seenReadIdsRef->{$readRef->{read}}=1;
    
    $statisticsRef->{'TOTAL'}++;
    
    # let's attempt to classify the read alignments
    my $numCandidates = scalar(@{$readRef->{candidates}});
    if (0==$numCandidates) {
        $statisticsRef->{'0-candidate'}++;
        print $fhExclude "# ", $readRef->{read}, " ", $numCandidates, " candidate\n";
        # TODO: write sam file as unmapped
        if (0!=$outputSam) {
            # there is nothing to load, it unmapped
            #LoadMAFRecord($readRef->{candidates}->[0]);
            printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, '', undef, '.', $numCandidates);
        }
        return;
    }
    if ($numCandidates>1) {
        $statisticsRef->{'Multi-candidates'}++;
        # TODO: DO NOT?? write sam file as multi-mapped
        # TODO: DO NOT?? write the first candidate
        printProfile (*OUTP, $readRef, $readRef->{cols}, $readRef->{candidates}->[0], 'MC', '', -1);
        if (0!=$outputSam) {
            # Need to load for cigar string
            LoadMAFRecord($readRef->{candidates}->[0]);
            # TODO: write multi-mapped as mapq=0
            printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, '', $readRef->{candidates}->[0], 'MC', $numCandidates);
        }
        return;
    }

    # there is only a single candidate!!!
    my $numFragments = scalar(@{$readRef->{candidates}->[0]});
    if (0==$numFragments) {
        $statisticsRef->{'Single-candidate 0-fragment'}++;
        print $fhExclude "# ", $readRef->{read}, " ", 0, " candidate\n";
        # write sam file as unmapped
        if (0!=$outputSam) {
            # there is nothing to load, it unmapped
            #LoadMAFRecord($readRef->{candidates}->[0]);
            printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, '', undef, '.', $numCandidates);
        }
        return;
    }
    if (1==$numFragments) {
        # there is no SV here
        $statisticsRef->{'Single-candidate single-fragment'}++;
        # write sam file
        # write the candidate
        printProfile (*OUTP, $readRef, $readRef->{cols}, $readRef->{candidates}->[0], 'SCSF', '', -1);
        if (0!=$outputSam) {
            # Need to load for cigar string
            LoadMAFRecord($readRef->{candidates}->[0]);
            # mapq = 60 given it is single-fragment
            printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, '', $readRef->{candidates}->[0], 'SCSF', $numCandidates);
        }
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
        $statisticsRef->{'Single-candidate multi-fragments multi-loci'}++;
        # write sam file as multi-mapped
        # write the first candidate
        printProfile (*OUTP, $readRef, $readRef->{cols}, $readRef->{candidates}->[0], 'SCMFML', join(",", @countCigar), -1);
        if (0!=$outputSam) {
            # Need to load for cigar string
            LoadMAFRecord($readRef->{candidates}->[0]);
            # TODO: write multi-mapped as mapq=0
            printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, '', $readRef->{candidates}->[0], 'SCMFML', $numCandidates);
        }
        return;
    }

    # question: can be inversion and TD at the same time?
    $statisticsRef->{'Single-candidate multi-fragments single-locus'}++;
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
	
    # write sam file
    # write the candidate
    printProfile (*OUTP, $readRef, $readRef->{cols}, $readRef->{candidates}->[0], 'SCMFSL', $marker, $labelIndicator);
    if (0!=$outputSam) {
        # Need to load for cigar string
        LoadMAFRecord($readRef->{candidates}->[0]);
        # TODO: write multi-mapped as mapq=0
        printSAM ($fhSam, $readRef->{read}, $readRef->{readseq}, $readRef->{readqual}, $readRef->{readsize}, $marker, $readRef->{candidates}->[0], 'SCMFSL', $numCandidates);
    }
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
