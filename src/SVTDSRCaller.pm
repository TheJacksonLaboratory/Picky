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

package SVTDSRCaller;

use strict;
use warnings;

our $VERSION = '0.1';

use base 'Exporter';
our @ISA = qw(SVCaller);


#####

my @ListPrefixCols = ('SVType', 'SVChrom', 'SVStart', 'SVEnd', 'SVSpan', 'qDiff', 'ReadId', 'ReadLen', 'ftotal', 'fid1', 'fid2');

#####

sub new {
	my $class = shift;
	my $parametersRef = shift;
	my $file = shift;
	
	# initialize the various parameters
	$parametersRef = {} if (!defined $parametersRef);
	@{$parametersRef->{'_listPrefixCols'}} = @ListPrefixCols;
	
	my $self =  $class->SUPER::new($parametersRef, $file, @_);
	$self->{'_SVType'} = 'TDSR';
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

	#sub containsTDSplitRead
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
	for(my $i=0; $i<($numAligns-1); ++$i) {
		my $j=$i+1;
		next if ($alignmentsRef->[$i]->{refId} ne $alignmentsRef->[$j]->{refId});
		next if ($alignmentsRef->[$i]->{refStrand} ne $alignmentsRef->[$j]->{refStrand});
		
		# same chromosome and same strand
		# check the order
		if ('+' eq $alignmentsRef->[$i]->{refStrand}) {
			if ($midPoints[$i]>$midPoints[$j]) {
				my $iAlignRef = $alignmentsRef->[$i];
				my $jAlignRef = $alignmentsRef->[$j];
				my ($refOverlap, $refILen, $refJLen) = utilities::computeIntervalCoverage ($iAlignRef->{refStart}, $iAlignRef->{refEnd}, $jAlignRef->{refStart}, $jAlignRef->{refEnd});
				my ($qOverlap, $qILen, $qJLen) = utilities::computeIntervalCoverage ($iAlignRef->{qStart}, $iAlignRef->{qEnd}, $jAlignRef->{qStart}, $jAlignRef->{qEnd});
				
				my $span = $alignmentsRef->[$i]->{refEnd} - $alignmentsRef->[$j]->{refStart}; #0-based; remove + 1
				
				my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$alignmentsRef->[$i]->{refId}, start=>$alignmentsRef->[$j]->{refStart}, end=>$alignmentsRef->[$i]->{refEnd},
				fragmentIds=>[$i, $j], qDiff=>-1*$qOverlap, sDiff=>-1*$refOverlap, svLen=>$span);
				# WCH : qDiff is added to get a sense of read coverage, and sDiff for reference coverage
				$svItem{sDiff} = $alignmentsRef->[$j]->{refStart} - $alignmentsRef->[$i]->{refEnd};
				$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
				# WCH: qDiff, sDiff - END
				$svItem{alignsRef} = $alignmentsRef;
				push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
				
				$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
				$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
			}
		} else {
			if ($midPoints[$j]>$midPoints[$i]) {
				my $iAlignRef = $alignmentsRef->[$i];
				my $jAlignRef = $alignmentsRef->[$j];
				my ($refOverlap, $refILen, $refJLen) = utilities::computeIntervalCoverage ($iAlignRef->{refStart}, $iAlignRef->{refEnd}, $jAlignRef->{refStart}, $jAlignRef->{refEnd});
				my ($qOverlap, $qILen, $qJLen) = utilities::computeIntervalCoverage ($iAlignRef->{qStart}, $iAlignRef->{qEnd}, $jAlignRef->{qStart}, $jAlignRef->{qEnd});
				
				my $span = $alignmentsRef->[$j]->{refEnd} - $alignmentsRef->[$i]->{refStart}; #0-based; remove + 1
				my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$alignmentsRef->[$i]->{refId}, start=>$alignmentsRef->[$i]->{refStart}, end=>$alignmentsRef->[$j]->{refEnd},
				fragmentIds=>[$i, $j], qDiff=>-1*$qOverlap, sDiff=>-1*$refOverlap, svLen=>$span);
				# WCH : qDiff is added to get a sense of read coverage, and sDiff for reference coverage
				$svItem{sDiff} = $alignmentsRef->[$i]->{refStart} - $alignmentsRef->[$j]->{refEnd};
				$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
				# WCH: qDiff, sDiff - END
				$svItem{alignsRef} = $alignmentsRef;
				push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
				
				$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
				$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
			}
		}
	}
	
	return $self->{'_numsvs'};
}

sub writeSVList {
	my $self = shift;
	#sub printTDSRList
	my $readRef = shift;
	#my $svsRef = shift;
	my $svsRef = $self->getSVs();
	#my $SVFHsRef
	
	my $fh = $self->{'_fh'};
	for(my $i=0; $i<scalar(@{$svsRef}); ++$i) {
		my $svRef = $svsRef->[$i];
		my @cols = ($svRef->{SVType});
		push @cols, $svRef->{chrom}, $svRef->{start}+1, $svRef->{end}, $svRef->{svLen}; # 0-based to 1-based
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
	#sub TDSRtoSAMCO
	#my $svsRef = shift;
	#my $svsRef = $self->getSVs();
	my $svsRef = $self->{'_svs'};
	
	my $retVal = $self->{'_SVType'}.'=';
	for(my $i=0; $i<$self->{'_numsvs'}; ++$i) {
		my $svRef = $svsRef->[$i];
		my $sOverQ = (0==$svRef->{qDiff}) ? -1*$svRef->{sDiff} : $svRef->{sDiff} / $svRef->{qDiff};
		$retVal .= sprintf("%sf%d_%d(q=%d,s=%d,s/q=%.2f,span=%d,loc=%s:%d-%d)",
		($i>0) ? ',' : '',
		$svRef->{fragmentIds}->[0]+1, $svRef->{fragmentIds}->[1]+1,
		-1*$svRef->{qDiff}, -1*$svRef->{sDiff},
		$sOverQ,
		$svRef->{svLen},
		$svRef->{chrom}, $svRef->{start}, $svRef->{end});
	}
	
	return $retVal;
}

1;

__END__
