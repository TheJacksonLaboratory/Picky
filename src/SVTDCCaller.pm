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

package SVTDCCaller;

use strict;
use warnings;

our $VERSION = '0.1';

use base 'Exporter';
our @ISA = qw(SVCaller);
use utilities;


#####

my @ListPrefixCols = ('SVType', 'SVChrom', 'SVStart', 'SVEnd', 'SVSpan', 'qDiff', 'sDiff', 'ReadId', 'ReadLen', 'ftotal', 'fid1', 'fid2');
my $G_MIN_TDC_SPAN_BP = 100;
my $G_MIN_TDC_SQDELTA_BP = $G_MIN_TDC_SPAN_BP;

#####

sub new {
	my $class = shift;
	my $parametersRef = shift;
	my $file = shift;
	
	# initialize the various parameters
	$parametersRef = {} if (!defined $parametersRef);
	@{$parametersRef->{'_listPrefixCols'}} = @ListPrefixCols;
	$parametersRef->{'_tdc_min_span'} = $G_MIN_TDC_SPAN_BP if (!exists $parametersRef->{'_tdc_min_span'});
	$parametersRef->{'_tdc_min_sqdelta'} = $G_MIN_TDC_SQDELTA_BP if (!exists $parametersRef->{'_tdc_min_sqdelta'});
	
	my $self =  $class->SUPER::new($parametersRef, $file, @_);
	$self->{'_SVType'} = 'TDC';
	delete $parametersRef->{'_listPrefixCols'};
	return $self;
}

sub DESTROY {
	my $self = shift;
	
	foreach my $class (@ISA) {
		my $destroy = "${class}::DESTROY";
		$self->$destroy if $self->can($destroy);
	}
}

sub containsSV {
	my $self = shift;
	
	#sub isTDRead
	my $readId = shift;
	my $readSeq = shift;
	my $readQual = shift;
	my $readSize = shift;
	my $alignmentsRef = shift;
	
	$self->{'_svs'}=[];
	$self->{'_numsvs'}=0;
	
	my $numAligns = scalar(@{$alignmentsRef});
	return if ($numAligns<2);
	
	my @midPoints = ();
	foreach my $alignRef (@{$alignmentsRef}) {
		my $midPoint = int(($alignRef->{refEnd} + $alignRef->{refStart})/2);
		push @midPoints, $midPoint;
	}
	
	# WCH: 2017-03-24
	#      if we perform a nested loop of (i,j), we can overreport the number of TDSR
	#      it is determine that we really only to detect immediate neighbor for split-read..
	#      so, it is swapped back to for(i){j=i+1}
	my @tdDetails = ();
	for(my $i=0; $i<($numAligns-1); ++$i) {
		my $j=$i+1;
		next if ($alignmentsRef->[$i]->{refId} ne $alignmentsRef->[$j]->{refId});
		next if ($alignmentsRef->[$i]->{refStrand} ne $alignmentsRef->[$j]->{refStrand});
		
		# same chromosome and same strand
		# check the order
		if ('+' eq $alignmentsRef->[$i]->{refStrand}) {
			if ($midPoints[$i]<$midPoints[$j]) {
				# let's check the overlap
				my $iAlignRef = $alignmentsRef->[$i];
				my $jAlignRef = $alignmentsRef->[$j];
				my ($refOverlap, $refILen, $refJLen) = utilities::computeIntervalCoverage ($iAlignRef->{refStart}, $iAlignRef->{refEnd}, $jAlignRef->{refStart}, $jAlignRef->{refEnd});
				
				if ($refOverlap>0) {
					my ($qOverlap, $qILen, $qJLen) = utilities::computeIntervalCoverage ($iAlignRef->{qStart}, $iAlignRef->{qEnd}, $jAlignRef->{qStart}, $jAlignRef->{qEnd});
					# compute the span
					# we have 2 cases here.. a) both have unique part ==> span okie, b) only one part is with unique part ==> span is uncertainty
					# the 3rd case is already handled by split read
					my $iUnique = $refILen - $refOverlap;
					my $jUnique = $refJLen - $refOverlap;
					my $span = 0;
					my $chrom = '';
					my $start = 0;
					my $end = 0;
					if ($iUnique>0) {
						if ($jUnique>0) {
							# both have unique span ==> we know the span
							# ---->--- ith
							#     --->---- jth
							$span = $alignmentsRef->[$i]->{refEnd} - $alignmentsRef->[$j]->{refStart};
							$chrom = $alignmentsRef->[$i]->{refId}; $start = $alignmentsRef->[$j]->{refStart}; $end = $alignmentsRef->[$i]->{refEnd};
						} else {
							# only ith has extra
							# ------- ith,  ------- ith
							#   ---   jth,    ----- jth
							$span = $alignmentsRef->[$i]->{refEnd} - $alignmentsRef->[$j]->{refStart};
							$chrom = $alignmentsRef->[$i]->{refId}; $start = $alignmentsRef->[$j]->{refStart}; $end = $alignmentsRef->[$i]->{refEnd};
						}
					} else {
						if ($jUnique>0) {
							# only jth has extra,
							#   ---   ith,  ----- ith
							# ------- jth,  ------- jth
							$span = $alignmentsRef->[$i]->{refEnd} - $alignmentsRef->[$j]->{refStart};
							$chrom = $alignmentsRef->[$i]->{refId}; $start = $alignmentsRef->[$j]->{refStart}; $end = $alignmentsRef->[$i]->{refEnd};
						} else {
							# only way to explain this is that both ith and jth fragments overlap exactly
							# there are 3 ways that this can be computed
							# ------- ith
							# ------- jth
							$span = $alignmentsRef->[$j]->{refEnd} - $alignmentsRef->[$i]->{refStart};
							$chrom = $alignmentsRef->[$i]->{refId}; $start = $alignmentsRef->[$i]->{refStart}; $end = $alignmentsRef->[$j]->{refEnd};
						}
					}
					# c = certain, u=uncertain/unsure
					if ($refOverlap-$qOverlap>$self->{'_tdc_min_sqdelta'} && $span>=$self->{'_tdc_min_span'}) {
						my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$chrom, start=>$start, end=>$end, fragmentIds=>[$i, $j], qDiff=>-1*$qOverlap, sDiff=>-1*$refOverlap, svLen=>$span);
						# WCH : qDiff is added to get a sense of read coverage, and sDiff for reference coverage
						my $sDiff = $alignmentsRef->[$j]->{refStart} - $alignmentsRef->[$i]->{refEnd};
						my $qDiff = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
						$svItem{sDiff} = $sDiff;
						$svItem{qDiff} = $qDiff;
						# WCH: qDiff, sDiff - END
						$svItem{alignsRef} = $alignmentsRef;
						push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
						
						$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
						$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
					}
				}
			}
		} else {
			if ($midPoints[$j]<$midPoints[$i]) {
				# let's check the overlap
				my $iAlignRef = $alignmentsRef->[$i];
				my $jAlignRef = $alignmentsRef->[$j];
				my ($refOverlap, $refILen, $refJLen) = utilities::computeIntervalCoverage ($iAlignRef->{refStart}, $iAlignRef->{refEnd}, $jAlignRef->{refStart}, $jAlignRef->{refEnd});
				
				if ($refOverlap>0) {
					my ($qOverlap, $qILen, $qJLen) = utilities::computeIntervalCoverage ($iAlignRef->{qStart}, $iAlignRef->{qEnd}, $jAlignRef->{qStart}, $jAlignRef->{qEnd});
					# see comments above for '+'ve strand cases
					# they have to be mirrored for '-'ve strand scenarios
					my $iUnique = $refILen - $refOverlap;
					my $jUnique = $refJLen - $refOverlap;
					my $span = 0;
					my $chrom = '';
					my $start = 0;
					my $end = 0;
					if ($iUnique>0) {
						if ($jUnique>0) {
							# both have unique span ==> we know the span
							# ----<--- jth
							#     ---<---- ith
							$span = $alignmentsRef->[$j]->{refEnd} - $alignmentsRef->[$i]->{refStart};
							$chrom = $alignmentsRef->[$i]->{refId}; $start = $alignmentsRef->[$i]->{refStart}; $end = $alignmentsRef->[$j]->{refEnd};
						} else {
							# only ith has extra
							# ------- ith,  ------- ith
							#   ---   jth,  ----- jth
							$span = $alignmentsRef->[$j]->{refEnd} - $alignmentsRef->[$i]->{refStart};
							$chrom = $alignmentsRef->[$i]->{refId}; $start = $alignmentsRef->[$i]->{refStart}; $end = $alignmentsRef->[$j]->{refEnd};
						}
					} else {
						if ($jUnique>0) {
							# only jth has extra,
							#   ---   ith,    ----- ith
							# ------- jth,  ------- jth
							$span = $alignmentsRef->[$j]->{refEnd} - $alignmentsRef->[$i]->{refStart};
							$chrom = $alignmentsRef->[$i]->{refId}; $start = $alignmentsRef->[$i]->{refStart}; $end = $alignmentsRef->[$j]->{refEnd};
						} else {
							# only way to explain this is that both ith and jth fragments overlap exactly
							# there are 3 ways that this can be computed
							# ------- ith
							# ------- jth
							$span = $alignmentsRef->[$j]->{refEnd} - $alignmentsRef->[$i]->{refStart};
							$chrom = $alignmentsRef->[$i]->{refId}; $start = $alignmentsRef->[$i]->{refStart}; $end = $alignmentsRef->[$j]->{refEnd};
						}
					}
					if ($refOverlap-$qOverlap>$self->{'_tdc_min_sqdelta'} && $span>=$self->{'_tdc_min_span'}) {
						my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$chrom, start=>$start, end=>$end, fragmentIds=>[$i, $j], qDiff=>-1*$qOverlap, sDiff=>-1*$refOverlap, svLen=>$span);
						# WCH : qDiff is added to get a sense of read coverage, and sDiff for reference coverage
						my $sDiff = $alignmentsRef->[$i]->{refStart} - $alignmentsRef->[$j]->{refEnd};
						my $qDiff = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
						$svItem{sDiff} = $sDiff;
						$svItem{qDiff} = $qDiff;
						# WCH: qDiff, sDiff - END
						$svItem{alignsRef} = $alignmentsRef;
						push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
						
						$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
						$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
					}
				}
			}
		}
	}
	
	return $self->{'_numsvs'};
}

sub writeSVList {
	my $self = shift;
	#sub printTDCList
	my $readRef = shift;
	#my $svsRef = shift;
	my $svsRef = $self->getSVs();
	#my $SVFHsRef
	
	my $fh = $self->{'_fh'};
	for(my $i=0; $i<scalar(@{$svsRef}); ++$i) {
		my $svRef = $svsRef->[$i];
		my @cols = ($svRef->{SVType});
		push @cols, $svRef->{chrom}, $svRef->{start}+1, $svRef->{end}, $svRef->{svLen}; # 0-based to 1-based
		push @cols, $svRef->{qDiff}, $svRef->{sDiff};
		
		push @cols, $readRef->{read}, $readRef->{readsize};
		push @cols, scalar(@{$svRef->{alignsRef}});
		grep {push @cols, $_+1} @{$svRef->{fragmentIds}};
		
		foreach my $mateId (@{$svRef->{fragmentIds}}) {
			my $alignRef = $svRef->{alignsRef}->[$mateId];
			grep { push @cols, $alignRef->{$_}; } @SVCaller::ReportPrefixCols;
			push @cols, $alignRef->{refStart}+1;
			grep { push @cols, $alignRef->{$_}; } @SVCaller::ReportSuffixCols;
		}
		print $fh join("\t", @cols), "\n";
	}
}

sub toCO {
	my $self = shift;
	#sub TDCtoSAMCO
	#my $svsRef = shift;
	#my $svsRef = $self->getSVs();
	my $svsRef = $self->{'_svs'};
	
	my $retVal = $self->{'_SVType'}.'=';
	for(my $i=0; $i<$self->{'_numsvs'}; ++$i) {
		my $svRef = $svsRef->[$i];
		$retVal .= sprintf("%sf%d_%d(q=%d,s=%d,span=%d,loc=%s:%d-%d)",
		($i>0) ? ',' : '',
		$svRef->{fragmentIds}->[0]+1, $svRef->{fragmentIds}->[1]+1,
		-1*$svRef->{qDiff}, -1*$svRef->{sDiff},
		$svRef->{svLen},
		$svRef->{chrom}, $svRef->{start}+1, $svRef->{end}); # 0-based to 1-based
	}
	
	return $retVal;
}

1;

__END__
