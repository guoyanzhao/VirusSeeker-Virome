################################################################################
##  Author:
##      Guoyan Zhao <gzhao@pathology.wustl.edu>
#*
#* Copyright (C) Washington University in St. Louis Developed by Guoyan Zhao,
#* Washington University in St. Louis, Department of Pathology and Immunology.
#*
###############################################################################

#!/usr/bin/perl
use strict;

my $usage = '
This script will read sequence in RepeatMasker.goodSeq.RefGenome.unmapped.unmasked.fasta
( fasta file with unmasked sequences), obtain the 
corresponding masked sequence in the file .QCed.cdhit.fa 
and output to a file RepeatMasker.goodSeq.RefGenome.unmapped.masked.fasta

perl script <sample dir> <sample name>
<sample dir> = full path of the folder holding files for this sample
               without last "/"
<sample name> = sample name

';
die $usage unless scalar @ARGV == 2;
my ( $dir, $SampleName ) = @ARGV;
my $percent_masked_cutoff = 0.4;


my $maskedFile = $dir."/".$SampleName.".RepeatMasker.goodSeq.masked.fa";
my %maskedSeq = ();
&read_FASTA_data($maskedFile, \%maskedSeq);

my $Unmasked_read_file = $dir."/".$SampleName."_VIRUSDB_HIT_Bacteria_unmapped.fasta";
my @Unmasked_read = ();
# get read ID
open (FastaFile, $Unmasked_read_file) or die "Can't Open FASTA file: $Unmasked_read_file";
while (my $line = <FastaFile>){
	if ($line =~ />/) {
		chomp $line;
		my @temp = split (/\s/, $line);
		my $seqName = shift @temp;
		$seqName =~ s/>//;
#		print "seq name is /$seqName/\n";
		push @Unmasked_read, $seqName;
	}
}
#&read_FASTA_data($Unmasked_read_file, \%Unmasked_read);

my $OutFile1 = $dir."/".$SampleName."_VIRUSDB_HIT_Bacteria_unmapped_masked.fasta";
open (OUT1, ">$OutFile1") or die "can not open $OutFile1\n";


foreach my $read_id ( @Unmasked_read) {
	print OUT1 ">$read_id\n";
	print OUT1 $maskedSeq{$read_id}, "\n";
}

close(OUT1);
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
		# discard comment lines
		elsif ($line =~ /#/) {
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
