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
This script will read the assignment report files in the given 
directory and generate a summary report for a given library. It will report 
in each library, for each category, how many total sequence were 
assigned to this category, how many were assigned by BLASTN, how many
were assigned by BLASTX.

It will also filter the virus lineage, leave out virus that are phage.
It will rank the virus lineage by range of percent ID from low to high. 

It will generate a .InterestingReads report about the details of each lineage
and the sequence of reads belong to this lineage.

perl script <sample folder> <sample name>
<sample folder> = full path to the folder holding files for a given sample 
<sample name> = sample name

';

die $usage unless scalar @ARGV == 2;
my ( $dir, $lib_name ) = @ARGV;



# generate viral reads assignment summary
my $out = $dir."/".$lib_name.".AssignmentSummary";
open (OUT, ">$out") or die "can not open file $out!\n";

my $seq_file = $dir."/".$lib_name.".RepeatMasker.goodSeq.unmasked.fa";
my %sequences = (); # read_ID => sequence
my %sequences_desc = (); # read_ID => description
&read_FASTA_data($seq_file, \%sequences, \%sequences_desc);
 
my %ID_low = ();  # lineage => lowest percent identity to hits
my %ID_high = (); # lineage => highest percent identity to hits

my $oldSeperator = $/;
$/ = "###########\n";
my $AssignmentReport_file = $dir."/".$lib_name.".AssignmentReport";
open (IN, $AssignmentReport_file) or die "can not open file $AssignmentReport_file!\n";
my $line = <IN>;
print OUT $line, "\n\n";

while (<IN>) { # read in every lineage information
	if ($_ =~ /^\s*$/) { # skip blank line
		next;
	}
	elsif ($_ =~ /Finished Assignment Report/) { next; }

	my @lines = split("\n", $_);
	my $lineage = shift @lines;
	$lineage = shift @lines;
	my $high = 0;
	my $low = 100;
	my %readID_Identity = (); # readID => percent ID
	my %readID_desc = (); # readID => description of the read
	for (my $i = 0; $i <= $#lines; $i++) {
#		print "line is \\$lines[$i]\\\n";
		if ($lines[$i] =~ /^\s*$/) { next; } # skip empty line
		elsif ($lines[$i] =~ /QueryName/) { next; } # skip description line
		elsif ($lines[$i] =~ /reads from/) { next; } # skip "reads from ..." line
		elsif ($lines[$i] =~ /^#+/) { next; } # skip # line
		elsif ($lines[$i] =~ /^>/) { 
			$i++;
			next; 
		}  # skip sequence line
		my ($read_ID, $Qlength, $hitName, $hitLen, $hitDesc, $alnLen, $ID, $hitS, $hitE, $e) = split("\t", $lines[$i]);
#		print "read_ID  = $read_ID, Qlength=$Qlength, hitName= $hitName, hitLen=$hitLen, hitDesc = \\$hitDesc\\, alnLen= $alnLen\\, ID = \\$ID\\, hitS= \\$hitS\\, hitE=\\$hitE\\, e = \\$e\\\n";
		$ID =~ s/%//;
		if($ID > $high) { $high = $ID;}
		if($ID < $low) { $low = $ID;}
		if (defined ($readID_Identity{$read_ID})) {
			if ($ID > $readID_Identity{$read_ID}) { 
				$readID_Identity{$read_ID} = $ID;
				$readID_desc{$read_ID} = $lines[$i];	
			}
		}
		else {
			$readID_Identity{$read_ID} = $ID;	
			$readID_desc{$read_ID} = $lines[$i];
		}
	}
	$ID_low{$lineage} = $low;
	$ID_high{$lineage} = $high;
}
close IN;

foreach my $key (sort {$ID_low{$a} <=> $ID_low{$b}} keys %ID_low) {
	printf OUT  ("%s\t[%4.1f, %4.1f]%\n", $key, $ID_low{$key}, $ID_high{$key});
}
print OUT "# Finished Assignment Summary\n";
close OUT;

########################################################################
# generate ambiguous reads assignment summary
my $out2 = $dir."/".$lib_name.".AmbiguousReadsAssignmentSummary";
open (OUT2, ">$out2") or die "can not open file $out2!\n";

my %ID_low_AmbiguousReads = ();  # lineage => lowest percent identity to hits
my %ID_high_AmbiguousReads = (); # lineage => highest percent identity to hits

my $AssignmentReport_file_AmbiguousReads = $dir."/".$lib_name.".AmbiguousReadsAssignmentReport";
open (IN2, $AssignmentReport_file_AmbiguousReads) or die "can not open file $AssignmentReport_file_AmbiguousReads!\n";
<IN2>; # read in the first line
while (<IN2>) {
	if ($_ =~ /^\s*$/) { # skip blank line
		next;
	}

	my @lines = split("\n", $_);
	my $lineage = shift @lines;
	$lineage = shift @lines;
	my $high = 0;
	my $low = 100;
	my %readID_Identity_AmbiguousReads = (); # readID => percent ID
	my %readID_desc_AmbiguousReads = (); # readID => description of the read
	for (my $i = 0; $i <= $#lines; $i++) {
#		print "line is \\$lines[$i]\\\n";
		if ($lines[$i] =~ /^\s*$/) { next; } # skip empty line
		elsif ($lines[$i] =~ /QueryName/) { next; } # skip description line
		elsif ($lines[$i] =~ /reads from/) { next; } # skip "reads from ..." line
		elsif ($lines[$i] =~ /^#+/) { next; } # skip # line
		elsif ($lines[$i] =~ /^>/) { 
			$i++;
			next; 
		}  # skip sequence line
		
		my ($read_ID, $Qlength, $hitName, $hitLen, $hitDesc, $alnLen, $ID, $hitS, $hitE, $e) = split("\t", $lines[$i]);
#		print "read_ID  = $read_ID, Qlength=$Qlength, hitName= $hitName, hitLen=$hitLen, hitDesc = \\$hitDesc\\, alnLen= $alnLen\\, ID = \\$ID\\, hitS= \\$hitS\\, hitE=\\$hitE\\, e = \\$e\\\n";
		$ID =~ s/%//;
		if($ID > $high) { $high = $ID;}
		if($ID < $low) { $low = $ID;}
		if (defined ($readID_Identity_AmbiguousReads{$read_ID})) {
			if ($ID > $readID_Identity_AmbiguousReads{$read_ID}) { 
				$readID_Identity_AmbiguousReads{$read_ID} = $ID;
				$readID_desc_AmbiguousReads{$read_ID} = $lines[$i];	
			}
		}
		else {
			$readID_Identity_AmbiguousReads{$read_ID} = $ID;	
			$readID_desc_AmbiguousReads{$read_ID} = $lines[$i];
		}
	}
#	print "lineage is \\$lineage\\\n";
#	print "low = $low, high = $high\n";

	$ID_low_AmbiguousReads{$lineage} = $low;
	$ID_high_AmbiguousReads{$lineage} = $high;
}
close IN;

foreach my $key (sort {$ID_low_AmbiguousReads{$a} <=> $ID_low_AmbiguousReads{$b}} keys %ID_low_AmbiguousReads) {
	printf OUT2  ("%s\t[%4.1f, %4.1f]%\n", $key, $ID_low_AmbiguousReads{$key}, $ID_high_AmbiguousReads{$key});
}
print OUT2 "# Finished Assignment Summary\n";
close OUT2;
$/ = $oldSeperator;

exit;

#####################################################################
sub read_FASTA_data () {
    my ($fastaFile, $seq_hash_ref, $seq_desc_hash_ref) = @_;

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
			my $temp = shift @rows;
			my @temp_arr = split(/\s/, $temp);
			my $contigName = shift @temp_arr;
			my $contigSeq = join("", @rows);
			$contigSeq =~ s/\s//g; #remove white space
			$seq_hash_ref->{$contigName} = $contigSeq;
			$seq_desc_hash_ref->{$contigName} = join(" ", @temp_arr);
#			print " name = \\$contigName\\, seq  = \\$contigSeq\\\n\n";
		}
    }

    # check
#	 foreach my $key (keys %fastaSeq){
#	      print "Here is the key for fasta seq: $key \t $fastaSeq{$key}\n";
#	 }

    #reset the read seperator
    $/ = $oldseperator;
}

