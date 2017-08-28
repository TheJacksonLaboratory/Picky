#####
#
# Picky - Structural Variants Pipeline for (ONT) long read
#
# Copyright (c) 2016-2017  Chee-Hong WONG  The Jackson Laboratory
#
#####

package SVTTLCCaller;

use strict;
use warnings;

our $VERSION = '0.1';

use base 'Exporter';
our @ISA = qw(SVCaller);


#####

my @ListPrefixCols = ('SVType', 'SVChrom1', 'SVPos1', 'SVStrand1', 'SVChrom2', 'SVPos2', 'SVStrand2', 'qDiff', 'ReadId', 'ReadLen', 'ftotal', 'fid1', 'fid2');

#####

sub new {
	my $class = shift;
	my $parametersRef = shift;
	my $file = shift;
	
	# initialize the various parameters
	$parametersRef = {} if (!defined $parametersRef);
	@{$parametersRef->{'_listPrefixCols'}} = @ListPrefixCols;
	
	my $self =  $class->SUPER::new($parametersRef, $file, @_);
	$self->{'_SVType'} = 'TTLC';
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
	
	#sub containsInterChromosomalTranslocation
	my $readId = shift;
	my $readSeq = shift;
	my $readQual = shift;
	my $readSize = shift;
	my $alignmentsRef = shift;
	
	$self->{'_svs'}=[];
	$self->{'_numsvs'}=0;
	
	my $numAligns = scalar(@{$alignmentsRef});
	return if ($numAligns<2);
	
	for(my $i=0; $i<($numAligns-1); ++$i) {
		my $j = $i + 1;
		my $iRef = $alignmentsRef->[$i];
		my $jRef = $alignmentsRef->[$j];
		
		if ($iRef->{refId} ne $jRef->{refId}) {
			my %svItem = (SVType=>$self->{'_SVType'}, fragmentIds=>[$i, $j], alignsRef=>$alignmentsRef);
			$svItem{chrom1} = $iRef->{refId}; $svItem{strand1} = $iRef->{refStrand}; $svItem{pos1} = ('+' eq $iRef->{refStrand}) ? $iRef->{'refEnd'} : $iRef->{'refStart'}+1;
			$svItem{chrom2} = $jRef->{refId}; $svItem{strand2} = $jRef->{refStrand}; $svItem{pos2} = ('+' eq $jRef->{refStrand}) ? $jRef->{'refStart'}+1 : $jRef->{'refEnd'};
			# WCH : qDiff is added to get a sense of read coverage
			$svItem{qDiff} = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd};
			# WCH: qDiff - END
			push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
			
			# TODO: why isn't this included for TTLC?
			# TODO: verify that TTLC is reported even when there is not {SVs} set!!!
			$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
			$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
			
		} else {
			# same chromosome!!!
		}
	}
	
	return $self->{'_numsvs'};
}

sub writeSVList {
	my $self = shift;
	#sub printTTLCList
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
	#sub TTLCtoSAMCO
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
