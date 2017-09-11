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

package SVINVCaller;

use strict;
use warnings;

our $VERSION = '0.1';

use base 'Exporter';
our @ISA = qw(SVCaller);
use utilities;

#####

my @ListPrefixCols = ('SVType', 'SVChrom', 'SVStart', 'SVEnd', 'SVSpan', 'InvFrag', 'SVChrom1', 'SVPos1', 'SVStrand1', 'SVChrom2', 'SVPos2', 'SVStrand2', 'qDiff', 'ReadId', 'ReadLen', 'ftotal', 'fid1', 'fid2');
my $G_INVERSION_MAX_SDIFF = 100000;

#####

sub new {
	my $class = shift;
	my $parametersRef = shift;
	my $file = shift;
	
	# initialize the various parameters
	$parametersRef = {} if (!defined $parametersRef);
	@{$parametersRef->{'_listPrefixCols'}} = @ListPrefixCols;
	$parametersRef->{'_inv_max_sdiff'} = $G_INVERSION_MAX_SDIFF if (!exists $parametersRef->{'_ctlc_min_sdiff'});
	
	my $self =  $class->SUPER::new($parametersRef, $file, @_);
	$self->{'_SVType'} = 'INV';
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
	
	#sub containsInversion
	my $readId = shift;
	my $readSeq = shift;
	my $readQual = shift;
	my $readSize = shift;
	my $alignmentsRef = shift;
	
	$self->{'_svs'}=[];
	$self->{'_numsvs'}=0;
	
	my $numAligns = scalar(@{$alignmentsRef});
	return if ($numAligns<2);
	
	my @invDetails = ();
	# WCH : we change from all-to-all comparison to sequential pair-wise
	for(my $i=0; $i<($numAligns-1); ++$i) {
		my $j=$i+1;
		next if ($alignmentsRef->[$i]->{refId} ne $alignmentsRef->[$j]->{refId});
		# same chromosome and same strand
		if ($alignmentsRef->[$i]->{refStrand} ne $alignmentsRef->[$j]->{refStrand}) {
			# TODO : provide alternative SV type
			my $iRef = $alignmentsRef->[$i];
			my $jRef = $alignmentsRef->[$j];
			my ($refOverlap, $refILen, $refJLen) = utilities::computeIntervalCoverage ($iRef->{refStart}, $iRef->{refEnd}, $jRef->{refStart}, $jRef->{refEnd});
			if ($refOverlap>0) {
				# this can no longer be an INVersion strictly
				# TODO: this should be reported as INS / breakpoint? however
				
				# for now, we are using the TTLC format
				my %svItem = (SVType=>'INVB', fragmentIds=>[$i, $j], alignsRef=>$alignmentsRef);
				$svItem{chrom1} = $iRef->{refId}; $svItem{strand1} = $iRef->{refStrand}; $svItem{pos1} = ('+' eq $iRef->{refStrand}) ? $iRef->{'refEnd'} : $iRef->{'refStart'}+1;
				$svItem{chrom2} = $jRef->{refId}; $svItem{strand2} = $jRef->{refStrand}; $svItem{pos2} = ('+' eq $jRef->{refStrand}) ? $jRef->{'refStart'}+1 : $jRef->{'refEnd'};
				# WCH : qDiff is added to get a sense of read coverage
				$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
				# WCH: qDiff - END
				#next if (abs($svItem{pos1}-$svItem{pos2})>$G_INTRA_TRANSLOCATION_MIN_SDIFF);
				next if (abs($svItem{pos1}-$svItem{pos2})>$self->{'_inv_max_sdiff'});
				push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
				
			} else {
				my $iMidPoint = int(($iRef->{refEnd} + $iRef->{refStart})/2);
				my $jMidPoint = int(($jRef->{refEnd} + $jRef->{refStart})/2);
				
				if ('+' eq $alignmentsRef->[$i]->{refStrand}) {
					# condition '+'/'-'
					if ($iMidPoint < $jMidPoint) {
						# --i>-- --<j-- : fragment(j) is the inverted sequence
						my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$jRef->{refId}, start=>$iRef->{refEnd}, end=>$jRef->{refEnd}, fragmentIds=>[$i, $j]);
						$svItem{sDiff} = $jRef->{refStart} - $iRef->{refEnd};
						# WCH : qDiff is added to get a sense of read coverage
						$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
						# WCH: qDiff - END
						$svItem{svLen} = $svItem{end} - $svItem{start};
						$svItem{invFrag} = $j;
						$svItem{alignsRef} = $alignmentsRef;
						# NOTE: CTLC take precedence, not deletion
						#next if ($svItem{svLen}>$G_INTRA_TRANSLOCATION_MIN_SDIFF);
						next if ($svItem{svLen}>$self->{'_inv_max_sdiff'});
						push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
					} else {
						# --<j-- --i>-- : fragment(i) is the inverted sequence
						my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$iRef->{refId}, start=>$jRef->{refEnd}, end=>$iRef->{refEnd}, fragmentIds=>[$i, $j]);
						$svItem{sDiff} = $iRef->{refStart} - $jRef->{refEnd};
						# WCH : qDiff is added to get a sense of read coverage
						$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
						# WCH: qDiff - END
						$svItem{svLen} = $svItem{end} - $svItem{start};
						$svItem{invFrag} = $i;
						$svItem{alignsRef} = $alignmentsRef;
						# NOTE: CTLC take precedence, not deletion
						#next if ($svItem{svLen}>$G_INTRA_TRANSLOCATION_MIN_SDIFF);
						next if ($svItem{svLen}>$self->{'_inv_max_sdiff'});
						push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
					}
				} else {
					# condition '-'/'+'
					if ($iMidPoint < $jMidPoint) {
						# --i<-- -->j-- : fragment(i) is the inverted sequence
						my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$iRef->{refId}, start=>$iRef->{refStart}, end=>$jRef->{refStart}, fragmentIds=>[$i, $j]);
						$svItem{sDiff} = $jRef->{refStart} - $iRef->{refEnd};
						# WCH : qDiff is added to get a sense of read coverage
						$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
						# WCH: qDiff - END
						$svItem{svLen} = $svItem{end} - $svItem{start};
						$svItem{invFrag} = $i;
						$svItem{alignsRef} = $alignmentsRef;
						# NOTE: CTLC take precedence, not deletion
						#next if ($svItem{svLen}>$G_INTRA_TRANSLOCATION_MIN_SDIFF);
						next if ($svItem{svLen}>$self->{'_inv_max_sdiff'});
						push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
					} else {
						# -->j-- --i<-- : fragment(j) is the inverted sequence
						my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$jRef->{refId}, start=>$jRef->{refStart}, end=>$iRef->{refStart}, fragmentIds=>[$i, $j]);
						$svItem{sDiff} = $iRef->{refStart} - $jRef->{refEnd};
						# WCH : qDiff is added to get a sense of read coverage
						$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
						# WCH: qDiff - END
						$svItem{svLen} = $svItem{end} - $svItem{start};
						$svItem{invFrag} = $j;
						$svItem{alignsRef} = $alignmentsRef;
						# NOTE: CTLC take precedence, not deletion
						#next if ($svItem{svLen}>$G_INTRA_TRANSLOCATION_MIN_SDIFF);
						next if ($svItem{svLen}>$self->{'_inv_max_sdiff'});
						push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
					}
				}
			}
			
			# TODO: do we need both i and j?
			$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
			$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
			$alignmentsRef->[$j]->{SVs} = {} if (!exists $alignmentsRef->[$j]->{SVs});
			$alignmentsRef->[$j]->{SVs}->{$self->{'_SVType'}} = 1;
		}
	}
	
	return $self->{'_numsvs'};
}

sub writeSVList {
	my $self = shift;
	#sub printINVList
	my $readRef = shift;
	#my $svsRef = shift;
	my $svsRef = $self->getSVs();
	#my $SVFHsRef
	
	my $fh = $self->{'_fh'};
	for(my $i=0; $i<scalar(@{$svsRef}); ++$i) {
		my $svRef = $svsRef->[$i];
		my @cols = ($svRef->{SVType});
		
		if ('INV' eq $svRef->{SVType}) {
			push @cols, $svRef->{chrom}, $svRef->{start}+1, $svRef->{end}, $svRef->{svLen}, $svRef->{invFrag}+1; # 0-based to 1-based
			push @cols, ('.') x 2;
			push @cols, $svRef->{alignsRef}->[$svRef->{fragmentIds}->[0]]->{refStrand};
			push @cols, ('.') x 2;
			push @cols, $svRef->{alignsRef}->[$svRef->{fragmentIds}->[1]]->{refStrand};
			push @cols, $svRef->{qDiff};
		} elsif ('INVB' eq $svRef->{SVType}) {
			push @cols, ('.') x 5;
			push @cols, $svRef->{chrom1}, $svRef->{pos1}, $svRef->{strand1};
			push @cols, $svRef->{chrom2}, $svRef->{pos2}, $svRef->{strand2};
			push @cols, $svRef->{qDiff};
		} else {
			push @cols, ('.') x 12;
		}
		
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
	#sub INVtoSAMCO
	#my $svsRef = shift;
	#my $svsRef = $self->getSVs();
	my $svsRef = $self->{'_svs'};
	
	# there are two sub-classes of INV
	my $retVal = $self->{'_SVType'}.'=';
	for(my $i=0; $i<$self->{'_numsvs'}; ++$i) {
		my $svRef = $svsRef->[$i];
		if ('INVB' eq $svRef->{SVType}) {
			$retVal .= sprintf("%sf%d_%d(loc=%s:%d:%s..%s:%d:%s)",
			($i>0) ? ',' : '',
			$svRef->{fragmentIds}->[0]+1, $svRef->{fragmentIds}->[1]+1,
			$svRef->{chrom1}, $svRef->{pos1}, $svRef->{strand1},
			$svRef->{chrom2}, $svRef->{pos2}, $svRef->{strand2});
		} else {
			$retVal .= sprintf("%sf%d_%d(frag=%d,rgap=%d,span=%d,loc=%s:%d-%d)",
			($i>0) ? ',' : '',
			$svRef->{fragmentIds}->[0]+1, $svRef->{fragmentIds}->[1]+1,
			$svRef->{invFrag}+1,
			$svRef->{sDiff},
			$svRef->{svLen},
			$svRef->{chrom}, $svRef->{start}, $svRef->{end});
		}
	}
	
	return $retVal;
}

1;

__END__
