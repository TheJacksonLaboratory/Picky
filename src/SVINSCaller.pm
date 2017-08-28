#####
#
# Picky - Structural Variants Pipeline for (ONT) long read
#
# Copyright (c) 2016-2017  Chee-Hong WONG  The Jackson Laboratory
#
#####

package SVINSCaller;

use strict;
use warnings;

our $VERSION = '0.1';

use base 'Exporter';
our @ISA = qw(SVCaller);
use utilities;


#####

my @ListPrefixCols = ('SVType', 'SVChrom', 'SVStart', 'SVEnd', 'SVSpan', 'qDiff', 'sDiff', 'ReadId', 'ReadLen', 'ftotal', 'fid1', 'fid2');
my $G_INS_MIN_QDIFF = 20; # qDiff >20, or >=21
my $G_INS_MAX_SDIFF = 20; # sDiff <20, or <=19

#####

sub new {
	my $class = shift;
	my $parametersRef = shift;
	my $file = shift;
	
	# initialize the various parameters
	$parametersRef = {} if (!defined $parametersRef);
	@{$parametersRef->{'_listPrefixCols'}} = @ListPrefixCols;
	$parametersRef->{'_ins_min_qdiff'} = $G_INS_MIN_QDIFF if (!exists $parametersRef->{'_ins_min_qdiff'});
	$parametersRef->{'_ins_max_sdiff'} = $G_INS_MAX_SDIFF if (!exists $parametersRef->{'_ins_max_sdiff'});
	
	my $self =  $class->SUPER::new($parametersRef, $file, @_);
	$self->{'_SVType'} = 'INS';
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
	
	#sub containsInsertionRead
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
		next if ($alignmentsRef->[$i]->{refStrand} ne $alignmentsRef->[$j]->{refStrand});
		
		# same chromosome and same strand
		# check that there is a gap between these two alignments
		if ('+' eq $alignmentsRef->[$i]->{refStrand}) {
			if ($alignmentsRef->[$i]->{qEnd}<$alignmentsRef->[$j]->{qStart}) {
				my $span = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd}; # 0-based; + 1
				if ($span>$self->{'_ins_min_qdiff'}) {
					my $refGap = $alignmentsRef->[$j]->{refStart} - $alignmentsRef->[$i]->{refEnd}; # 0-based; + 1
					if ($refGap<$self->{'_ins_max_sdiff'}) {
						# NOTE: If this is a TDC/TDSR, there is no need to indicate that it is an INS anymore
						if (exists $alignmentsRef->[$i]->{SVs}) {
							next if (exists $alignmentsRef->[$i]->{SVs}->{TDC});
							next if (exists $alignmentsRef->[$i]->{SVs}->{TDSR});
						}
						my $start = undef; my $end = undef;
						if ($refGap<=0) {
							# insertion happens at the end of fragment(i)
							$start = $alignmentsRef->[$i]->{refEnd}; $end = $alignmentsRef->[$i]->{refEnd}+1;
						} else {
							# there is a segment that is replaced by the insertion
							$start = $alignmentsRef->[$i]->{refEnd}; $end = $alignmentsRef->[$j]->{refStart};
						}
						
						my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$alignmentsRef->[$i]->{refId}, start=>$start, end=>$end,
						fragmentIds=>[$i, $j], qDiff=>$span, sDiff=>$refGap, svLen=>$span);
						$svItem{alignsRef} = $alignmentsRef;
						push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
						
						$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
						$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
					}
				}
			}
		} else {
			if ($alignmentsRef->[$i]->{qEnd}<$alignmentsRef->[$j]->{qStart}) {
				my $span = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd}; # 0-based; + 1
				if ($span>$self->{'_ins_min_qdiff'}) {
					my $refGap = $alignmentsRef->[$i]->{refStart} - $alignmentsRef->[$j]->{refEnd}; # 0-based; + 1
					if ($refGap<$self->{'_ins_max_sdiff'}) {
						# NOTE: If this is a TDC/TDSR, there is no need to indicate that it is an INS anymore
						if (exists $alignmentsRef->[$i]->{SVs}) {
							next if (exists $alignmentsRef->[$i]->{SVs}->{TDC});
							next if (exists $alignmentsRef->[$i]->{SVs}->{TDSR});
						}
						my $start = undef; my $end = undef;
						if ($refGap<=0) {
							# insertion happens at the end of fragment(j)
							# NOTE: should use fragment(i) but fragment(j) will be the same as '+' set up
							$start = $alignmentsRef->[$j]->{refEnd}; $end = $alignmentsRef->[$j]->{refEnd}+1;
						} else {
							# there is a segment that is replaced by the insertion
							$start = $alignmentsRef->[$j]->{refEnd}; $end = $alignmentsRef->[$i]->{refStart};
						}
						
						my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$alignmentsRef->[$i]->{refId}, start=>$start, end=>$end,
						fragmentIds=>[$i, $j], qDiff=>$span, sDiff=>$refGap, svLen=>$span);
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
	#sub printINSList
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
	#sub INStoSAMCO
	#my $svsRef = shift;
	#my $svsRef = $self->getSVs();
	my $svsRef = $self->{'_svs'};
	
	my $retVal = $self->{'_SVType'}.'=';
	for(my $i=0; $i<$self->{'_numsvs'}; ++$i) {
		my $svRef = $svsRef->[$i];
		my $sOverQ = (0==$svRef->{qDiff}) ? -1*$svRef->{sDiff} : $svRef->{sDiff} / $svRef->{qDiff};
		$retVal .= sprintf("%sf%d_%d(rgap=%d,span=%d,loc=%s:%d-%d)",
		($i>0) ? ',' : '',
		$svRef->{fragmentIds}->[0]+1, $svRef->{fragmentIds}->[1]+1,
		$svRef->{sDiff},
		$svRef->{svLen},
		$svRef->{chrom}, $svRef->{start}, $svRef->{end});
	}
	
	return $retVal;
}

1;

__END__
