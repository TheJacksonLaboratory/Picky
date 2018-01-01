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

package SVINDELCaller;

use strict;
use warnings;

our $VERSION = '0.1';

use base 'Exporter';
our @ISA = qw(SVCaller);
use utilities;


#####

my @ListPrefixCols = ('SVType', 'SVChrom', 'SVStart', 'SVEnd', 'SVSpan', 'qDiff', 'sDiff', 'ReadId', 'ReadLen', 'ftotal', 'fid1', 'fid2');
my $G_INDEL_MIN_QDIFF = 20; # qDiff >20, or >=21   ## SVDELCaller.pm: my $G_DEL_MIN_SDIFF = 20; # sDiff >20, or >=21
my $G_INDEL_MIN_SDIFF = 20; # sDiff >20, or >=21   ## SVDELCaller.pm: my $G_DEL_MAX_QDIFF = 20; # qDiff <20, or <=19
my $G_REMOVE_HOMOPOLYMER_ARTIFACT = 0;

#####

sub new {
	my $class = shift;
	my $parametersRef = shift;
	my $file = shift;
	
	# initialize the various parameters
	$parametersRef = {} if (!defined $parametersRef);
	@{$parametersRef->{'_listPrefixCols'}} = @ListPrefixCols;
	$parametersRef->{'_indel_min_qdiff'} = $G_INDEL_MIN_QDIFF if (!exists $parametersRef->{'_indel_min_qdiff'});
	$parametersRef->{'_indel_min_sdiff'} = $G_INDEL_MIN_SDIFF if (!exists $parametersRef->{'_indel_min_sdiff'});
	$parametersRef->{'_remove_homopolymer_artifact'} = $G_REMOVE_HOMOPOLYMER_ARTIFACT if (!exists $parametersRef->{'_remove_homopolymer_artifact'});
	$parametersRef->{'_genome_sequence_hash'} = undef if (!exists $parametersRef->{'_genome_sequence_hash'});
	
	my $self =  $class->SUPER::new($parametersRef, $file, @_);
	$self->{'_SVType'} = 'INDEL';
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
	
	#sub containsInsertionAndDeletionRead
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
			if ($alignmentsRef->[$i]->{refEnd}<$alignmentsRef->[$j]->{refStart}) {
				my $span = $alignmentsRef->[$j]->{refStart} - $alignmentsRef->[$i]->{refEnd}; #0-based; remove + 1
				if ($span>=$self->{'_indel_min_sdiff'}) {
					my $readGap = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd}; #0-based; remove + 1
					if ($readGap>=$self->{'_indel_min_qdiff'}) {
						if ($readGap<=$self->{'_ctlc_min_sdiff'}) {
							# NOTE: CTLC take precedence, not deletion
							next if (exists $alignmentsRef->[$i]->{SVs} && exists $alignmentsRef->[$i]->{SVs}->{CTLC});
							
							# check for homopolymer issue
							if (0!=$self->{'_remove_homopolymer_artifact'}) {
								my $affected = utilities::isAffectedByHomopolymer($readId, $readSeq, $readQual, $readSize, $alignmentsRef->[$i], $alignmentsRef->[$j], $self->{'_genome_sequence_hash'}, $self);
								next if (0!=$affected);
							}
							
							my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$alignmentsRef->[$i]->{refId}, start=>$alignmentsRef->[$i]->{refEnd}, end=>$alignmentsRef->[$j]->{refStart},
							fragmentIds=>[$i, $j], qDiff=>$readGap, sDiff=>$span, svLen=>$span);
							$svItem{alignsRef} = $alignmentsRef;
							push @{$self->{'_svs'}}, \%svItem; $self->{'_numsvs'}++;
							
							$alignmentsRef->[$i]->{SVs} = {} if (!exists $alignmentsRef->[$i]->{SVs});
							$alignmentsRef->[$i]->{SVs}->{$self->{'_SVType'}} = 1;
						}
					}
				}
			}
		} else {
			if ($alignmentsRef->[$j]->{refEnd}<$alignmentsRef->[$i]->{refStart}) {
				my $span = $alignmentsRef->[$i]->{refStart} - $alignmentsRef->[$j]->{refEnd}; #0-based; remove + 1
				if ($span>=$self->{'_indel_min_sdiff'}) {
					my $readGap = $alignmentsRef->[$j]->{qStart} - $alignmentsRef->[$i]->{qEnd}; #0-based; remove + 1
					if ($readGap>=$self->{'_indel_min_qdiff'}) {
						if ($readGap<=$self->{'_ctlc_min_sdiff'}) {
							# NOTE: CTLC take precedence, not deletion
							next if (exists $alignmentsRef->[$i]->{SVs} && exists $alignmentsRef->[$i]->{SVs}->{CTLC});
							
							# check for homopolymer issue
							if (0!=$self->{'_remove_homopolymer_artifact'}) {
								my $affected = utilities::isAffectedByHomopolymer($readId, $readSeq, $readQual, $readSize, $alignmentsRef->[$i], $alignmentsRef->[$j], $self->{'_genome_sequence_hash'}, $self);
								next if (0!=$affected);
							}
							
							my %svItem = (SVType=>$self->{'_SVType'}, chrom=>$alignmentsRef->[$i]->{refId}, start=>$alignmentsRef->[$j]->{refEnd}, end=>$alignmentsRef->[$i]->{refStart},
							fragmentIds=>[$i, $j], qDiff=>$readGap, sDiff=>$span, svLen=>$span);
							$svItem{alignsRef} = $alignmentsRef;
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
	#sub printINDELList
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
	#sub INDELtoSAMCO
	#my $svsRef = shift;
	#my $svsRef = $self->getSVs();
	my $svsRef = $self->{'_svs'};
	
	my $retVal = $self->{'_SVType'}.'=';
	for(my $i=0; $i<$self->{'_numsvs'}; ++$i) {
		my $svRef = $svsRef->[$i];
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
