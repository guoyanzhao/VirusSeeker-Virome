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
This script will read the phage assignment report files in the given 
directory and generate a phage summary report for a given library. 
It will rank the phage lineage by range of percent ID from low to high. 

perl script <sample folder> <sample name>
<sample folder> = full path to the folder holding files for a given sample 
<sample name> = sample name

';

die $usage unless scalar @ARGV == 2;
my ( $dir, $lib_name ) = @ARGV;


# generate viral reads assignment summary
my $out = $dir."/".$lib_name.".PhageAssignmentSummary";
open (OUT, ">$out") or die "can not open file $out!\n";

=head1
my $seq_file = $dir."/".$lib_name.".QCed.cdhit.fa";
my %sequences = (); # read_ID => sequence
my %sequences_desc = (); # read_ID => description
&read_FASTA_data($seq_file, \%sequences, \%sequences_desc);
=cut

my %ID_low = ();  # lineage => lowest percent identity to hits
my %ID_high = (); # lineage => highest percent identity to hits

my $oldSeperator = $/;
$/ = "###########\n";
my $AssignmentReport_file = $dir."/".$lib_name.".PhageAssignmentReport";
open (IN, $AssignmentReport_file) or die "can not open file $AssignmentReport_file!\n";
my $line = <IN>;
print OUT $line, "\n";

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

$/ = $oldSeperator;

exit;

=head1
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
=cut
