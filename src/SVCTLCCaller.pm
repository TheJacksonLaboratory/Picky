#####
#
# Picky - Structural Variants Pipeline for (ONT) long read
#
# Copyright (c) 2016-2017  Chee-Hong WONG  The Jackson Laboratory
#
#####

package SVCTLCCaller;

use strict;
use warnings;

our $VERSION = '0.1';

use base 'Exporter';
our @ISA = qw(SVCaller);


#####

my @ListPrefixCols = ('SVType', 'SVChrom1', 'SVPos1', 'SVStrand1', 'SVChrom2', 'SVPos2', 'SVStrand2', 'qDiff', 'ReadId', 'ReadLen', 'ftotal', 'fid1', 'fid2');
##my $G_INTRA_TRANSLOCATION_MIN_SDIFF = 100000; # sdiff >100000, or >=100001
my $G_INTRA_TRANSLOCATION_MIN_SDIFF = 1000000000; # set to 1000 Mbp so that it will not be triggered

#####

sub new {
	my $class = shift;
	my $parametersRef = shift;
	my $file = shift;
	
	# initialize the various parameters
	$parametersRef = {} if (!defined $parametersRef);
	@{$parametersRef->{'_listPrefixCols'}} = @ListPrefixCols;
	$parametersRef->{'_ctlc_min_sdiff'} = $G_INTRA_TRANSLOCATION_MIN_SDIFF if (!exists $parametersRef->{'_ctlc_min_sdiff'});
	
	my $self =  $class->SUPER::new($parametersRef, $file, @_);
	$self->{'_SVType'} = 'CTCL';
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
	
	#sub containsIntraChromosomalTranslocation
	my $readId = shift;
	my $readSeq = shift;
	my $readQual = shift;
	my $readSize = shift;
	my $alignmentsRef = shift;
	
	$self->{'_svs'}=[];
	$self->{'_numsvs'}=0;
	
	my $numAligns = scalar(@{$alignmentsRef});
	return if ($numAligns<2);
	
	# check that we are on the same chromosome, same strand, and correct order!!!
	for(my $i=0; $i<($numAligns-1); ++$i) {
		my $j = $i+1;
		next if ($alignmentsRef->[$i]->{refId} ne $alignmentsRef->[$j]->{refId});
		if ($alignmentsRef->[$i]->{refStrand} eq $alignmentsRef->[$j]->{refStrand}) {
			# handle large deletion
			# handle large indel
			# same chromosome and same strand
			# check that there is a gap between these two alignments
			if ('+' eq $alignmentsRef->[$i]->{refStrand}) {
				if ($alignmentsRef->[$i]->{refEnd}<$alignmentsRef->[$j]->{refStart}) {
					my $span = $alignmentsRef->[$j]->{refStart} - $alignmentsRef->[$i]->{refEnd}; #0-based; remove + 1
					if ($span>$self->{'_ctlc_min_sdiff'}) {
						my $readGap = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd}; #0-based; remove + 1
						my $iRef = $alignmentsRef->[$i];
						my $jRef = $alignmentsRef->[$j];
						my %svItem = (SVType=>$self->{'_SVType'}, fragmentIds=>[$i, $j], alignsRef=>$alignmentsRef);
						$svItem{chrom1} = $iRef->{refId}; $svItem{strand1} = $iRef->{refStrand}; $svItem{pos1} = ('+' eq $iRef->{refStrand}) ? $iRef->{'refEnd'} : $iRef->{'refStart'}+1;
						$svItem{chrom2} = $jRef->{refId}; $svItem{strand2} = $jRef->{refStrand}; $svItem{pos2} = ('+' eq $jRef->{refStrand}) ? $jRef->{'refStart'}+1 : $jRef->{'refEnd'};
						
						# WCH : qDiff is added to get a sense of read coverage
						$svItem{sDiff} = $span;
						$svItem{qDiff} = $readGap;
						
						$svItem{alignsRef} = $alignmentsRef;
						push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
						
						$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
						$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
					}
				}
			} else {
				if ($alignmentsRef->[$j]->{refEnd}<$alignmentsRef->[$i]->{refStart}) {
					my $span = $alignmentsRef->[$i]->{refStart} - $alignmentsRef->[$j]->{refEnd}; #0-based; remove + 1
					if ($span>$self->{'_ctlc_min_sdiff'}) {
						my $readGap = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd}; #0-based; remove + 1
						my $iRef = $alignmentsRef->[$i];
						my $jRef = $alignmentsRef->[$j];
						my %svItem = (SVType=>$self->{'_SVType'}, fragmentIds=>[$i, $j], alignsRef=>$alignmentsRef);
						$svItem{chrom1} = $iRef->{refId}; $svItem{strand1} = $iRef->{refStrand}; $svItem{pos1} = ('+' eq $iRef->{refStrand}) ? $iRef->{'refEnd'} : $iRef->{'refStart'}+1;
						$svItem{chrom2} = $jRef->{refId}; $svItem{strand2} = $jRef->{refStrand}; $svItem{pos2} = ('+' eq $jRef->{refStrand}) ? $jRef->{'refStart'}+1 : $jRef->{'refEnd'};
						
						# WCH : qDiff is added to get a sense of read coverage
						$svItem{sDiff} = $span;
						$svItem{qDiff} = $readGap;
						
						$svItem{alignsRef} = $alignmentsRef;
						push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
						
						$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
						$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
					}
				}
			}
		} else {
			# different strand... possible inversion
			# using the mutually exclusive span range defined in containsInversion(...)
			
			my $iRef = $alignmentsRef->[$i];
			my $jRef = $alignmentsRef->[$j];
			my ($refOverlap, $refILen, $refJLen) = utilities::computeIntervalCoverage ($iRef->{refStart}, $iRef->{refEnd}, $jRef->{refStart}, $jRef->{refEnd});
			if ($refOverlap<=0) {
				my $iMidPoint = int(($iRef->{refEnd} + $iRef->{refStart})/2);
				my $jMidPoint = int(($jRef->{refEnd} + $jRef->{refStart})/2);
				
				if ('+' eq $alignmentsRef->[$i]->{refStrand}) {
					# condition '+'/'-'
					if ($iMidPoint < $jMidPoint) {
						# --i>-- --<j-- : fragment(j) is the inverted sequence ... but considered as ctlc
						my %svItem = (SVType=>$self->{'_SVType'}, fragmentIds=>[$i, $j], alignsRef=>$alignmentsRef);
						$svItem{chrom1} = $iRef->{refId}; $svItem{strand1} = $iRef->{refStrand}; $svItem{pos1} = ('+' eq $iRef->{refStrand}) ? $iRef->{'refEnd'} : $iRef->{'refStart'}+1;
						$svItem{chrom2} = $jRef->{refId}; $svItem{strand2} = $jRef->{refStrand}; $svItem{pos2} = ('+' eq $jRef->{refStrand}) ? $jRef->{'refStart'}+1 : $jRef->{'refEnd'};
						$svItem{sDiff} = $jRef->{refStart} - $iRef->{refEnd};
						# WCH : qDiff is added to get a sense of read coverage
						$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
						# WCH: qDiff - END
						$svItem{svLen} = abs($svItem{pos1} - $svItem{pos2});
						#$svItem{invFrag} = $j;
						$svItem{alignsRef} = $alignmentsRef;
						if ($svItem{svLen}>$self->{'_ctlc_min_sdiff'}) {
							push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
							
							$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
							$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
						}
					} else {
						# --<j-- --i>-- : fragment(i) is the inverted sequence ... but considered as ctlc
						my %svItem = (SVType=>$self->{'_SVType'}, fragmentIds=>[$i, $j], alignsRef=>$alignmentsRef);
						$svItem{chrom1} = $iRef->{refId}; $svItem{strand1} = $iRef->{refStrand}; $svItem{pos1} = ('+' eq $iRef->{refStrand}) ? $iRef->{'refEnd'} : $iRef->{'refStart'}+1;
						$svItem{chrom2} = $jRef->{refId}; $svItem{strand2} = $jRef->{refStrand}; $svItem{pos2} = ('+' eq $jRef->{refStrand}) ? $jRef->{'refStart'}+1 : $jRef->{'refEnd'};
						$svItem{sDiff} = $iRef->{refStart} - $jRef->{refEnd};
						# WCH : qDiff is added to get a sense of read coverage
						$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
						# WCH: qDiff - END
						$svItem{svLen} = abs($svItem{pos1} - $svItem{pos2});
						#$svItem{invFrag} = $i;
						$svItem{alignsRef} = $alignmentsRef;
						if ($svItem{svLen}>$self->{'_ctlc_min_sdiff'}) {
							push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
							
							$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
							$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
						}
					}
				} else {
					# condition '-'/'+'
					if ($iMidPoint < $jMidPoint) {
						# --i<-- -->j-- : fragment(i) is the inverted sequence ... but considered as ctlc
						my %svItem = (SVType=>$self->{'_SVType'}, fragmentIds=>[$i, $j], alignsRef=>$alignmentsRef);
						$svItem{chrom1} = $iRef->{refId}; $svItem{strand1} = $iRef->{refStrand}; $svItem{pos1} = ('+' eq $iRef->{refStrand}) ? $iRef->{'refEnd'} : $iRef->{'refStart'}+1;
						$svItem{chrom2} = $jRef->{refId}; $svItem{strand2} = $jRef->{refStrand}; $svItem{pos2} = ('+' eq $jRef->{refStrand}) ? $jRef->{'refStart'}+1 : $jRef->{'refEnd'};
						$svItem{sDiff} = $jRef->{refStart} - $iRef->{refEnd};
						# WCH : qDiff is added to get a sense of read coverage
						$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
						# WCH: qDiff - END
						$svItem{svLen} = abs($svItem{pos1} - $svItem{pos2});
						#$svItem{invFrag} = $i;
						$svItem{alignsRef} = $alignmentsRef;
						if ($svItem{svLen}>$self->{'_ctlc_min_sdiff'}) {
							push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
							
							$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
							$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
						}
					} else {
						# -->j-- --i<-- : fragment(j) is the inverted sequence ... but considered as ctlc
						my %svItem = (SVType=>$self->{'_SVType'}, fragmentIds=>[$i, $j], alignsRef=>$alignmentsRef);
						$svItem{chrom1} = $iRef->{refId}; $svItem{strand1} = $iRef->{refStrand}; $svItem{pos1} = ('+' eq $iRef->{refStrand}) ? $iRef->{'refEnd'} : $iRef->{'refStart'}+1;
						$svItem{chrom2} = $jRef->{refId}; $svItem{strand2} = $jRef->{refStrand}; $svItem{pos2} = ('+' eq $jRef->{refStrand}) ? $jRef->{'refStart'}+1 : $jRef->{'refEnd'};
						$svItem{sDiff} = $iRef->{refStart} - $jRef->{refEnd};
						# WCH : qDiff is added to get a sense of read coverage
						$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
						# WCH: qDiff - END
						$svItem{svLen} = abs($svItem{pos1} - $svItem{pos2});
						#$svItem{invFrag} = $j;
						$svItem{alignsRef} = $alignmentsRef;
						if ($svItem{svLen}>$self->{'_ctlc_min_sdiff'}) {
							push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
							
							$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
							$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
						}
					}
				}
			}
		}
	}
	
	return $self->{'_numsvs'};
}

sub writeSVList {
	my $self = shift;
	#sub printCTLCList
	my $readRef = shift;
	#my $svsRef = shift;
	my $svsRef = $self->getSVs();
	#my $SVFHsRef
	
	my $fh = $self->{'_fh'};
	for(my $i=0; $i<scalar(@{$svsRef}); ++$i) {
		my $svRef = $svsRef->[$i];
		my @cols = ($svRef->{SVType});
		push @cols, $svRef->{chrom1}, $svRef->{pos1}, $svRef->{strand1};
		push @cols, $svRef->{chrom2}, $svRef->{pos2}, $svRef->{strand2};
		push @cols, $svRef->{qDiff};
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
	#sub CTLCtoSAMCO
	#my $svsRef = shift;
	#my $svsRef = $self->getSVs();
	my $svsRef = $self->{'_svs'};
	
	my $retVal = $self->{'_SVType'}.'=';
	for(my $i=0; $i<$self->{'_numsvs'}; ++$i) {
		my $svRef = $svsRef->[$i];
		$retVal .= sprintf("%sf%d_%d(loc=%s:%d:%s..%s:%d:%s)",
		($i>0) ? ',' : '',
		$svRef->{fragmentIds}->[0]+1, $svRef->{fragmentIds}->[1]+1,
		$svRef->{chrom1}, $svRef->{pos1}, $svRef->{strand1},
		$svRef->{chrom2}, $svRef->{pos2}, $svRef->{strand2});
	}
	
	return $retVal;
}

1;

__END__
