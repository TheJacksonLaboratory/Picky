#####
#
# Picky - Structural Variants Pipeline for (ONT) long read
#
# Copyright (c) 2016-2017  Chee-Hong WONG  The Jackson Laboratory
#
#####

package hashFastq;

BEGIN {
	$VERSION = '0.1';
}

#####
# hashFastq.pm
#
# To assign a human-friendly read id in place of the unique 128-bits number
# (aka read_uuid) assigned by the MinKNOW software for each read.
#
# The human-friendly read id composed of 7 elements as follow:
#   Element 1: experiment id; e.g. WTD09
#   Element 2: unique running serial number from 1 assigned to each read; e.g. 1..n
#   Element 3: "1D" or "2D" status of read; e.g. 2D
#   Element 4: "P"ass or "F"ail status of read;	P
#   Element 5: Read code where
#              1 = Read template sequence is present
#              2 = Read complement sequence is present
#              3 = Both read template and complement sequence are present
#   Element 6: Read sequence code where
#              T=Template
#              C=Complement
#              M=Merged consensus 2D sequence
#   Element 7: Read sequence length prefixed with 'L'; e.g. L23284
#
#
# Input:
#   pass read fastq file
#   fail read fastq file
#
# Output:
#   <oprefix>.2D.fastq
#   <oprefix>.2Dsupport.fastq
#   <oprefix>.1D.fastq
#   <oprefix>.strange.fastq
#   Hash table (<oprefix>-hash.xls)
#     Column 1: read_uuid
#     Column 2: fast5
#     Column 3: PassFail
#     Column 4: 1D_2D
#     Column 5: template_hash_id
#     Column 6: complement_hash_id
#     Column 7: 2d_hash_id
#
#####

use strict;
use Getopt::Long;
use base 'Exporter';
our @EXPORT = qw(runHashFastq);

my $G_USAGE = "
$0 hashFq --pfile <passFQFile> --ffile <failFQFile> --oprefix <outputPrefix>

  --pfile STR    pass .fastq file
  --ffile STR    fail .fastq file
  --oprefix STR  prefix to output filename
";

sub runHashFastq {
	my $pfile = undef;
	my $ffile = undef;
	my $outputPrefix = undef;
	
	GetOptions (
	"pfile=s" => \$pfile,
	"ffile=s"   => \$ffile,
	"oprefix=s"   => \$outputPrefix,
	) or die("Error in command line arguments\n$G_USAGE");
	
	die "Please specify the output prefix\n$G_USAGE" if (!defined $outputPrefix);
	die "Please specify either the pass or fail fastq file\n$G_USAGE" if (!defined $pfile && !defined $ffile);
	
	# 2-pass approach
	# 1st pass: gather attributes
	# 2nd pass: to write out the fastq files
	
	my $logFile = sprintf("%s.run.txt", $outputPrefix);
	open LFILE, ">$logFile" || die "Fail to open $logFile\n$!\n";
	
	# 1st pass
	my %reads = ();
	parseONTPoreToolsFastq($pfile, 'P', \%reads) if (defined $pfile);
	parseONTPoreToolsFastq($ffile, 'F', \%reads) if (defined $ffile);
	
	# output a master list
	my $hashFile = sprintf("%s-hash.xls", $outputPrefix);
	open HFILE, ">$hashFile" || die "Fail to open $hashFile\n$!\n";
	my @hashCols = ('read_uuid', 'fast5', 'PassFail', '1D_2D', 'template_hash_id', 'complement_hash_id', '2d_hash_id');
	print HFILE '#', join("\t", @hashCols), "\n";
	
	# 2nd pass
	my $usable2DFile = sprintf("%s.2D.fastq", $outputPrefix);
	my $supportFile = sprintf("%s.2Dsupport.fastq", $outputPrefix);
	my $usable1DFile = sprintf("%s.1D.fastq", $outputPrefix);
	my $strangeFile = sprintf("%s.strange.fastq", $outputPrefix);
	
	open U2DFILE, ">$usable2DFile" || die "Fail to open $usable2DFile\n$!\n";
	open SFILE, ">$supportFile" || die "Fail to open $supportFile\n$!\n";
	open U1DFILE, ">$usable1DFile" || die "Fail to open $usable1DFile\n$!\n";
	open EFILE, ">$strangeFile" || die "Fail to open $strangeFile\n$!\n";
	
	# favor 'pass' read
	secondPassONTPoreToolsFastq($pfile, 'P', $outputPrefix, \%reads) if (defined $pfile);
	# include 'fail' read
	secondPassONTPoreToolsFastq($ffile, 'F', $outputPrefix, \%reads) if (defined $ffile);
	
	close LFILE;
	close HFILE;
	close U2DFILE;
	close SFILE;
	close U1DFILE;
	close EFILE;
}

# skip those which we have skipped in the first pass
sub secondPassONTPoreToolsFastq {
    my ($file, $passfail, $libId, $readsRef) = @_;

	my @hashColIDs = ('read_uuid', 'fast5', 'PassFail', '1D_2D', 'template_hash_id', 'complement_hash_id', '2d_hash_id');
	my %hashRec = ();
	my $currHashId = undef;
    open INFILE, $file || die "Fail to open $file";
    while (<INFILE>) {
        my $id = $_; chomp($id);
        my $seq = <INFILE>; chomp($seq);
        my $ignore = <INFILE>;
        my $qual = <INFILE>; chomp($qual);
        
        $id =~ s/^@//;
        my @ids = split(/\s+/, $id);
        my @readids = split(/\_/, $ids[0]);
        my $readType = lc($readids[-1]);
        my ($channel, $readid) = $ids[1] =~ /_ch(\d+)_read(\d+)_strand\d*/;
        die "Fail to locate channel in '$ids[1]'" if (!defined $channel);
        die "Fail to locate read id in '$ids[1]'" if (!defined $readid);
		
        my $readuid = $readids[0];
        die "'$readuid' does not exist?!" if (!exists $readsRef->{$readuid});
        
		my @paths = split(/\//, $ids[2]);
		my $fast5 = $paths[-1];
		die "Fast5 file cannot be located in 3rd component $id!\n" if ($fast5 !~/\.fast5$/);
		$fast5 =~ s/\.fast5$//;
		
        # check the list to get the various code to output
        # A) 2D-pass, 2D-fail
        # B) 2D-pass-template, 2D-pass-complement, 2D-fail-template, 2D-fail-complement
        # C) 1D-pass-template, 1D-pass-complement
        # D) 1D-fail-template, 1D-fail-complement
        # E) 1D-pass-template & 1D-pass-complement: QUESTION: why not 2D-pass???
        # F) 1D-fail-template & 1D-fail-complement: QUESTION: why not 2D-fail???
        #
        # A), C), D) as a group ==> 2Dn1D
        # B) is for checking result of A) ==> 2Dsupport
        # E), F) consider as erroneous cases - to exclude ==> strange
        #
        # A), B) : R01:1:2D:P:3:T, R01:1:2D:P:3:C, R01:1:2D:P:3:M
        # A), B) : R01:2:2D:F:3:T, R01:2:2D:F:3:C, R01:2:2D:F:3:M
        # C) R01:3:1D:P:1:T, R01:4:1D:P:2:C
        # D) R01:5:1D:F:1:T, R01:6:1D:F:2:C
        # E) R01:7:1D:P:3:T, R01:7:1D:P:3:C
        # F) R01:8:1D:F:3:T, R01:8:1D:F:3:C
        #
        # format: <run_id>:<read_id>:<1D/2D>:<P_ass/F_ail>:<Template/Complement_bit>:<T_emplate,C_omplement,M_erge>
        #

        my $readRecRef = $readsRef->{$readuid};
        die "'$readuid''s $readType does not exist?!" if (!exists $readRecRef->{$readType});
        
        my $seqLen = length($seq);
		
		# col 1: read_uuid
		# col 2: fast5
		# col 3: PassFail
		# col 4: 1D_2D
		# col 5: template_hash_id
		# col 6: complement_hash_id
		# col 7: 2d_hash_id

		if ($currHashId ne $readuid) {
			if (defined $currHashId) {
				# write the hash record
				my @cols = (); grep { push @cols, $hashRec{$_}; } @hashColIDs;
				print HFILE join("\t", @cols), "\n";
			}
			grep { $hashRec{$_} = '.'; } @hashColIDs;
			$hashRec{'read_uuid'} = $readuid;
			$hashRec{'fast5'} = $fast5;
			$hashRec{'PassFail'} = $passfail;
			$hashRec{'1D_2D'} = (exists $readRecRef->{'2d'}) ? '2D' : '1D';
			$currHashId = $readuid
		}
		
        my $hashid = undef;
        if ('2d' eq $readType) {
            # A)
            if (!exists $readRecRef->{$readType.'_p'}) {
                $hashid = sprintf("%s:%d:2D:%s:3:M:L%d", $libId, $readRecRef->{readoffset}, $passfail, $seqLen);
                print U2DFILE '@',$hashid,"\n";
                print U2DFILE $seq, "\n";
                print U2DFILE "+\n";
                print U2DFILE $qual, "\n";
				
				$hashRec{'2d_hash_id'} = $hashid;
            } else {
                print LFILE "Pass#2: Read type '$readType' already processed! Skipping.\n";
            }
        } elsif ('template' eq $readType) {
            if (exists $readRecRef->{'2d'}) {
                # B)
                if (!exists $readRecRef->{$readType.'_p'}) {
                    $hashid = sprintf("%s:%d:2D:%s:3:T:L%d", $libId, $readRecRef->{readoffset}, $passfail, $seqLen);
                    print SFILE '@',$hashid,"\n";
                    print SFILE $seq, "\n";
                    print SFILE "+\n";
                    print SFILE $qual, "\n";
                    
					$hashRec{'template_hash_id'} = $hashid;
                } else {
                    print LFILE "Pass#2: Read type '$readType' already processed! Skipping.\n";
                }
            } elsif ('complement' eq $readType) {
                # E), F)
                if (!exists $readRecRef->{$readType.'_p'}) {
                    $hashid = sprintf("%s:%d:1D:%s:3:T:L%d", $libId, $readRecRef->{readoffset}, $passfail, $seqLen);
                    print EFILE '@',$hashid,"\n";
                    print EFILE $seq, "\n";
                    print EFILE "+\n";
                    print EFILE $qual, "\n";
                    
					$hashRec{'complement_hash_id'} = $hashid;
                } else {
                    print LFILE "Pass#2: Read type '$readType' already processed! Skipping.\n";
                }
            } else {
                # C), D)
                if (!exists $readRecRef->{$readType.'_p'}) {
                    $hashid = sprintf("%s:%d:1D:%s:1:T:L%d", $libId, $readRecRef->{readoffset}, $passfail, $seqLen);
                    print U1DFILE '@',$hashid,"\n";
                    print U1DFILE $seq, "\n";
                    print U1DFILE "+\n";
                    print U1DFILE $qual, "\n";
                    
					$hashRec{'template_hash_id'} = $hashid;
                } else {
                    print LFILE "Pass#2: Read type '$readType' already processed! Skipping.\n";
                }
            }
        } elsif ('complement' eq $readType) {
            if (exists $readRecRef->{'2d'}) {
                # B)
                if (!exists $readRecRef->{$readType.'_p'}) {
                    $hashid = sprintf("%s:%d:2D:%s:3:C:L%d", $libId, $readRecRef->{readoffset}, $passfail, $seqLen);
                    print SFILE '@',$hashid,"\n";
                    print SFILE $seq, "\n";
                    print SFILE "+\n";
                    print SFILE $qual, "\n";
                    
					$hashRec{'complement_hash_id'} = $hashid;
                } else {
                    print LFILE "Pass#2: Read type '$readType' already processed! Skipping.\n";
                }
            } elsif ('template' eq $readType) {
                # E), F)
                if (!exists $readRecRef->{$readType.'_p'}) {
                    $hashid = sprintf("%s:%d:1D:%s:3:C:L%d", $libId, $readRecRef->{readoffset}, $passfail, $seqLen);
                    print EFILE '@',$hashid,"\n";
                    print EFILE $seq, "\n";
                    print EFILE "+\n";
                    print EFILE $qual, "\n";
                    
					$hashRec{'template_hash_id'} = $hashid;
                } else {
                    print LFILE "Pass#2: Read type '$readType' already processed! Skipping.\n";
                }
            } else {
                # C), D)
                if (!exists $readRecRef->{$readType.'_p'}) {
                    $hashid = sprintf("%s:%d:1D:%s:2:C:L%d", $libId, $readRecRef->{readoffset}, $passfail, $seqLen);
                    print U1DFILE '@',$hashid,"\n";
                    print U1DFILE $seq, "\n";
                    print U1DFILE "+\n";
                    print U1DFILE $qual, "\n";
                    
					$hashRec{'complement_hash_id'} = $hashid;
                } else {
                    print LFILE "Pass#2: Read type '$readType' already processed! Skipping.\n";
                }
            }
        } else {
            die "Unrecognized read type '$readType' for '$readuid'\n";
        }
    }
    close INFILE;
	if (defined $currHashId) {
		# write the hash record
		my @cols = (); grep { push @cols, $hashRec{$_}; } @hashColIDs;
		print HFILE join("\t", @cols), "\n";
	}
}

sub parseONTPoreToolsFastq {
    my ($file, $passfail, $readsRef) = @_;
    
    $readsRef->{'G_Count'} = 0 if (! exists $readsRef->{'G_Count'});
    my $counter = $readsRef->{'G_Count'};
    
    open INFILE, $file || die "Fail to open $file";
    while (<INFILE>) {
        my $id = $_; chomp($id);
        my $seq = <INFILE>; chomp($seq);
        my $ignore = <INFILE>;
        my $qual = <INFILE>; chomp($qual);

        $id =~ s/^@//;
        my @ids = split(/\s+/, $id);
        my @readids = split(/\_/, $ids[0]);
        my $readType = lc($readids[-1]);
        my ($channel, $readid) = $ids[1] =~ /_ch(\d+)_read(\d+)_strand\d*/;
        die "Fail to locate channel in '$ids[1]'" if (!defined $channel);
        die "Fail to locate read id in '$ids[1]'" if (!defined $readid);
        
        my $readuid = $readids[0];
        if (!exists $readsRef->{$readuid}) {
            $counter++;
            $readsRef->{$readuid} = {readuid=>$readuid, channel=>$channel, readid=>$readid, passfail=>$passfail, readoffset=>$counter};
        } else {
            # let's check that we have the same values!!!
            my $readRecRef = $readsRef->{$readuid};
            if (exists $readRecRef->{'2d'}) {
                print LFILE "Pass#1: Read uid '$readuid' already has a 2D record!\n";
            }
            if ($readRecRef->{channel} ne $channel) {
                print LFILE "Pass#1: Different channel '$ids[1]' vs '$readRecRef->{channel}' (exp) in line '$_'";
                next;
            }
            if ($readRecRef->{readid} ne $readid) {
                print LFILE "Pass#1: Different read id '$ids[1]' vs '$readRecRef->{readid}' (exp) in line '$_'";
                next;
            }
            if ($readRecRef->{passfail} ne $passfail) {
                print LFILE "Pass#1: Different read type '$passfail' vs '$readRecRef->{passfail}' (exp), keeping read type '$passfail'\n\tline '$_'";
                next;
            }
        }
        
        my $readRecRef = $readsRef->{$readuid};
        if (exists $readRecRef->{$readType}) {
            print LFILE "Pass#1: Read type '$readType' already exists! Skipping.\n";
        } else {
            $readRecRef->{$readType} = {};
        }
    }
    close INFILE;
    
    $readsRef->{'G_Count'} = $counter;
}


1;

__END__
