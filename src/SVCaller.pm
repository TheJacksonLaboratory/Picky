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

package SVCaller;

use strict;
use warnings;

our $VERSION = '0.1';

use base 'Exporter';

#####

our @AlignmentCols = ('score', 'EG2', 'E', 'Identity', 'IdentityP', 'Mismatch', 'MismatchP', 'Deletion', 'DeletionP', 'Insertion', 'InsertionP', 'qStart', 'qEnd', 'qStrand', 'qALen', 'qP', 'refId', 'refStart', 'refEnd', 'refStrand', 'refALen');
our @ReportPrefixCols = ('score', 'EG2', 'E', '=', '%=', 'X', '%X', 'D', '%D', 'I', '%I', 'qStart', 'qEnd', 'qStrand', 'qALen', 'q%', 'refId');
our @ReportSuffixCols = ('refEnd', 'refStrand', 'refALen');

my @ListPrefixCols = ('SVType', 'sameChrom', 'sameStrand', 'strandCode', 'qdiff', 'sdiff', 'ReadId', 'ReadLen', 'ftotal', 'fid1', 'fid2');

#####


#####
# base class - all specific-SV Caller should dervied from
# - set up the list file ; can have different headers
# - containsSV()
# - writeSVList()
# - toCO()

# new(\%parameters, $filename)
sub new {
	my $class = shift;
	my $parametersRef = shift;
	my $file = shift;
	
	my $self = {'_SVType'=>'NONE', '_file'=>undef, '_fh'=>undef, '_svs'=>[], '_numsvs'=>0};
	
	# initialize the various parameters
	while (my ($key, $valueRef) = each %{$parametersRef}) {
		$self->{$key} = $valueRef;
	}
	if ((!defined $parametersRef) || (!exists $parametersRef->{'_listPrefixCols'})) {
		@{$self->{'_listPrefixCols'}} = @ListPrefixCols;
	} else {
		@{$self->{'_listPrefixCols'}} = @{$parametersRef->{'_listPrefixCols'}};
	}
	
	# set up the file if necessary
	if (defined $file) {
		$self->{'_file'} = $file;
		open my $fh, ">$file" || die "Fail to open $file\n$!\n";
		$self->{'_fh'} = $fh;
	
		# WCH : qDiff is added to get a sense of read coverage, and reference coverage (sDiff)
		my @cols = @{$self->{'_listPrefixCols'}};
		grep { push @cols, $_.'_1' } @AlignmentCols;
		grep { push @cols, $_.'_2' } @AlignmentCols;
		print $fh '#', join("\t", @cols), "\n";
	}
	
	bless $self, $class;
	return $self;
}

sub DESTROY {
	my $self = shift;
	$self->{'_fh'}->close() if ($self->{'_fh'});
}

# base class does not detect any SV and thus classified as 'NONE'
sub containsSV {
	my $self = shift;
	my $readId = shift;
	my $readSeq = shift;
	my $readQual = shift;
	my $readSize = shift;
	my $alignmentsRef = shift;
	
	$self->{'_svs'}=[];
	$self->{'_numsvs'}=0;
	return 0;
}

sub writeSVList {
	my $self = shift;
	#sub printNONEList
	#my ($readRef, $alignmentsRef, $SVFHsRef) = @_;
	my $readRef = shift;
	my $alignmentsRef = shift;
	
	my $numAligns = scalar(@{$alignmentsRef});
	return if ($numAligns<2);
	
	for(my $i=0; $i<($numAligns-1); ++$i) {
		my $j = $i + 1;
		my $iRef = $alignmentsRef->[$i];
		my $jRef = $alignmentsRef->[$j];
		
		my $toCompute = 1;
		my @cols = ($self->{'_SVType'});
		if ($iRef->{refId} ne $jRef->{refId}) {
			# diff chrom
			$toCompute = 0;
			push @cols, 'N';
		} else {
			push @cols, 'Y';
		}
		if ($iRef->{refStrand} ne $jRef->{refStrand}) {
			# diff strand
			push @cols, 'N';
		} else {
			push @cols, 'Y';
		}
		push @cols, $iRef->{refStrand}.$jRef->{refStrand};
		
		my $qdiff = 'n.a.';
		my $sdiff = 'n.a.';
		if (0!=$toCompute) {
			# query diff
			$qdiff = $jRef->{qStart} - $iRef->{qEnd};
			# subject diff
			if ('+' eq $iRef->{refStrand}) {
				if ('+' eq $jRef->{refStrand}) {
					# +/+
					$sdiff = $jRef->{refStart} - $iRef->{refEnd};
				} else {
					# +/-
					if ($iRef->{refStart}<=$jRef->{refStart}) {
						$sdiff = $jRef->{refStart} - $iRef->{refStart};
					} else {
						$sdiff = $jRef->{refEnd} - $iRef->{refEnd};
					}
				}
			} else {
				if ('+' eq $jRef->{refStrand}) {
					# -/+
					if ($iRef->{refStart}<=$jRef->{refStart}) {
						$sdiff = $jRef->{refStart} - $iRef->{refStart};
					} else {
						$sdiff = $jRef->{refEnd} - $iRef->{refEnd};
					}
				} else {
					# -/-
					$sdiff = $iRef->{refStart} - $jRef->{refEnd};
				}
			}
		}
		push @cols, $qdiff, $sdiff;
		
		push @cols, $readRef->{read}, $readRef->{readsize};
		push @cols, $numAligns, $i+1, $j+1;
		
		foreach my $alignRef ($iRef, $jRef) {
			grep { push @cols, $alignRef->{$_}; } @ReportPrefixCols;
			push @cols, $alignRef->{refStart}+1;
			grep { push @cols, $alignRef->{$_}; } @ReportSuffixCols;
		}
		my $fh = $self->{'_fh'};
		print $fh join("\t", @cols), "\n";
	}
}

sub getSVType {
	my $self = shift;
	return $self->{'_SVType'};
}

sub getSVs {
	my $self = shift;
	return $self->{'_svs'};
}

sub getNumberOfSVs {
	my $self = shift;
	return $self->{'_numsvs'};
}

sub toCO {
	my $self = shift;
	return '';
}

1;

__END__
