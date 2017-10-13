################################################################################
##  Author:
##      Guoyan Zhao <gzhao@pathology.wustl.edu>
#******************************************************************************
#* Copyright (C) Washington University in St. Louis Developed by Guoyan Zhao,
#* Washington University in St. Louis, Department of Pathology and Immunology.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php.
#*
###############################################################################


#!/usr/bin/perl
use strict;

my $usage = '
This script will check the Tanta masked file.
Some sequences have only/lots of Ns because masked by Tantan.
1) Sequences that do not have greater than 50 nt of consecutive
sequence without N will be put into file .cdhit_out.masked.badSeq 
2) Sequences with >= 40% of total length of being masked will be put 
into file .cdhit_out.masked.RepeatLowComplexSeq
3) output original sequnce into .cdhit_out.masked.goodSeq file.

perl script <dir> <sample name>
<dir> = full path to the file without the last "/"
<sample name> = sample name

';
die $usage unless scalar @ARGV == 2;
my ( $dir, $sample_name ) = @ARGV;

####################################################################################################
# get name of the sample and path to the data
if ($dir =~/(.+)\/$/) {
	$dir = $1;
}

my $original_fa_file = $dir."/".$sample_name.".QCed.cdhit.tantan.fa";
my %original_fa_seq = ();
&read_FASTA_data($original_fa_file, \%original_fa_seq);

my $total_seq = 0;
my $good_seq = 0;
my $bad_seq = 0;
my $OutFile1 = $dir."/".$sample_name.".QCed.cdhit.tantan.goodSeq.fa";
my $OutFile2 = $dir."/".$sample_name.".QCed.cdhit.tantan.badSeq.fa";

open (OUT1, ">$OutFile1") or die "can not open $OutFile1\n";
open (OUT2, ">$OutFile2") or die "can not open $OutFile2\n";

my $maskedFile = $dir."/".$sample_name.".QCed.cdhit.tantan.fa";;
my %seq = ();
&read_FASTA_data($maskedFile, \%seq);

# check for contiguous bases >= 50 bp (non-Ns) 
foreach my $read_id (keys %seq) {
	$total_seq++;
	my $seq_temp = $seq{$read_id};
	my $goodQuality=$seq_temp=~/[ACTG]{50,}/; 
	if($goodQuality) {
		print OUT1 ">$read_id\n";
#		print OUT1 $seq{$read_id}, "\n";
		print OUT1 $original_fa_seq{$read_id}, "\n";
		$good_seq++;
	}
	else { # output bad sequence
		print OUT2 ">$read_id\n";
		print  OUT2 "$seq{$read_id}\n";
		$bad_seq++;
	}
}

print OUT2 "total seq = $total_seq\n";
print OUT2 "good seq = $good_seq\n";
print OUT2 "bad seq = $bad_seq\n";

close(OUT1);
close(OUT2);

exit;

############################################################################
sub read_FASTA_data () {
    my ($fastaFile, $hash_ref) = @_;

    #keep old read seperator and set new read seperator to ">"
    my $oldseperator = $/;
    $/ = ">";
	 
    open (FastaFile, $fastaFile) or die "Can't Open FASTA file: $fastaFile";
    while (my $line = <FastaFile>){
		# Discard blank lines
        if ($line =~ /^\s*$/) {
			next;
		}
		# discard the first line which only has ">", keep the rest
		elsif ($line ne ">") {
			chomp $line;
		    my @rows = ();
		    @rows = split (/\n/m, $line);	
		    my $seqName = shift @rows;
			my @temp = split (/\s/, $seqName);
			$seqName = shift @temp;
		    my $Seq = join("", @rows);
		    $Seq =~ s/\s//g; #remove white space
		    $hash_ref->{$seqName} = $Seq;
#			print "name = $seqName\n";
#			print "seq = \\$Seq\\\n";
		}
    }

	close FastaFile;
    #reset the read seperator
    $/ = $oldseperator;
}
