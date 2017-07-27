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

my $usage = "
This script will read corresponding files in the given director and 
generate a report which contains SampleDescription, SequenceReport,
AssignmentSummary, InterestingReads.

perl $0 <run folder> <program version>
<run folder> = full path of the folder holding files for this sequence run
<program version> = e.g. 0.063

";
die $usage unless scalar @ARGV == 2;
my ( $dir, $version ) = @ARGV;

my @temp = split("/", $dir);
my $run_name = pop @temp;
my $outFile = $dir."/Analysis_Report_".$run_name.".txt";
open (OUT, ">$outFile") or die "can not open file $outFile!\n";

my ($wkday,$month,$day,$time,$year) = split(/\s+/, localtime);
print OUT "Illumina V${version}; Processing date: $day-$month-$year\n";
print OUT "cd-hit parameter: -c 0.95 -n 8 -G 0 -aS 0.95 -g 1 -r 1 -M 9000 -d 0\n";
print OUT "MegaBLAST against ref genome database parameter: -evalue 1e-9 -word_size 16  -show_gis  \n";
print OUT "BLASTn against NT VIRUSDB database parameter: -evalue 1e-4 -show_gis -task blastn \n";
print OUT "BLASTx against NR VIRUSDB parameter: -evalue 1e-2 -show_gis  \n";
print OUT "MegaBLAST against full NT database parameter:\n";
print OUT " -evalue 1e-9 -word_size 16  -show_gis -num_descriptions 50 -num_alignments 50  \n";
print OUT "BLASTn against nt database parameter:\n";
print OUT " -evalue 1e-9 -show_gis -task blastn -num_descriptions 50 -num_alignments 50   \n";
print OUT "BLASTx against nr database parameter: \n";
print OUT  "-evalue 1e-2 -show_gis -num_descriptions 50 -num_alignments 50   \n";
print OUT "database version: Human genomic: human_g1k_v37.fasta\n";
print OUT "NR, NT, taxonomy: 20160802\n";
print OUT "E-value cutoff for significant hit:\n";
print OUT "BLASTn VIRUSDB 1e-5, MegaBLAST ref genome, MegaBLAST NT: e-10\n";
print OUT "BLASTn NT: e-10, BLASTx VIRUSDB and BLASTx NR: e-3\n";

my $c = "**************************************************************************\n";
my $c2 = "#########################################################################\n\n";
print OUT $c;
###################################################
print OUT "How to read this file:\n\n";

print OUT "For the Summary section:\n";
print OUT "column 1: sample name\n";
print OUT "column 2: total number of reads obtained for this sample (pairs)\n\n";

print OUT "If there is any viral sequence identified in this sample, it will show up under \n";
print OUT "the information of this sample. There are 3 columns to describe the viral sequences identified in this sample:\n";
print OUT "column 1: number of viral sequences (could be read or contig)\n";
print OUT "column 2: range of percentage identity to blast hits. Some times one sequence can hit multiple \n";
print OUT "sequences in the database with the same e value. Only the first best hit percent identity is taken.\n";
print OUT "column 3: name of the virus (the last rank in the taxonomy lineage of the virus)\n\n";
 
###################################################
print OUT "For the Sequence Report section: It was divided into three subsections.\n\n";

print OUT "First Subsection:\n";
print OUT "column 1: #totalPair: total number of reads in pairs obtained for this sample.\n";
print OUT "column 2: #unjoined: number of reads not joined by fastq-join\n";
print OUT "column 3: #joined: number of reads joined by fastq-join \n";
print OUT "column 4: %total: percentage of joined reads over total pair of reads\n";
print OUT "column 5: #TotalAfterStitch: number of reads (unjoined+joined) after fastq-join \n";
print OUT "column 6: #QC: number of reads left after quality control (removing low quality nt and read, low complexity read, exact duplicates etc)\n";
print OUT "column 7: %Stitched: percentage of reads left after quality control (reads left after quality control divided by total number of reads after join read 1 and 2 (column 5)\n";
print OUT "column 8: #unique: number of unique reads after removing redundancy using CD-HIT\n";
print OUT "column 9: %Stitched: percentage of deduplicated reads (reads left after deduplication divided by total number of reads after join read 1 and 2 (column 5)\n";
print OUT "column 10: sample: name of the sample\n\n";

print OUT "Second Subsection:\n";
print OUT "column 1: #TANTAN_good: number of good quality reads after TANTAN filtering low complexity sequence.\n";
print OUT "column 2: %unique: percentage of #TANTAN_good reads over total number of deduplicated reads (column 8 in the first subsection\n";
print OUT "column 3: #RM_good: number of good quality reads after RepeatMasker filtering low complexity sequence.\n";
print OUT "column 4: %unique: percentage of #RM_good reads over total number of deduplicated reads (column 8 in the first subsection\n";
print OUT "column 5: #RefGenomeFiltered: number of reads after reference genome filtering.\n";
print OUT "column 6: %unique: percentage of #RM_good reads over total number of deduplicated reads (column 8 in the first subsection\n";
print OUT "column 7: #BN_VIRUSDB_Hit: number of reads have significant hit in BLASTn Virus-only database\n";
print OUT "column 8: #BX_VIRUSDB_Hit: number of reads have significant hit in BLASTx Virus-only database\n";
print OUT "column 9: #Total_VIRUSDB_Hit: number of reads have significant hit in both BLASTn and BLASTx Virus-only database (candidate eukaryotic viral sequences)\n";
print OUT "column 10: %unique: percentage of #RM_good reads over total number of deduplicated reads (column 8 in the first subsection\n";
print OUT "column 11: VIRUSDB_Hit_Bac: number of candidate eukaryotic viral reads mapped to bacteria\n";
print OUT "column 12: %VIRUSDB_Hit: percentage of VIRUSDB_Hit_Bac over total number of candidate eukaryotic viral sequences (column 9\n";
print OUT "column 13: VIRUSDB_Hit_NON_Bac: number of candidate eukaryotic viral reads not mapped to bacteria\n";
print OUT "column 14: %VIRUSDB_Hit: percentage of VIRUSDB_Hit_NON_Bac over total number of candidate eukaryotic viral sequences (column 9\n";
print OUT "column 15: sample: name of the sample\n\n";


print OUT "Third Subsection: information about BLAST results.\n";
print OUT "column 1: #BNinput: number of reads entered BLASTn NT step\n";
print OUT "column 2: %VIRUSDB_Hit: percentage of #BNinput over total number of candidate \n eukaryotic viral sequences (column 3)\n";
print OUT "column 3: BXinput: number of reads entered BLASTx NR step\n";
print OUT "column 4: %VIRUSDB_Hit: percentage of BXinput over total number of candidate \n eukaryotic viral sequences (column 3)\n";
print OUT "column 5: Unassigned: number of reads that do not have sequence similarity with \n anything in the database.\n";
print OUT "column 6: %VIRUSDB_Hit: percentage of Unassigned over total number of candidate \n eukaryotic viral sequences (column 3)\n";
print OUT "column 7: sample name \n\n";


###################################################
print OUT "For the Taxonomy Assignment section:\n"; 
print OUT "Describes the number of sequences assigned to each category.\n\n";
print OUT "If any viral sequence is detected in this sample, one line of description for each viral lineage will\n";
print OUT "be under the assignment information. The description includes the full taxonomy lineage of the virus, \n";
print OUT "the total number of reads that are assigned to this taxonomy lineage and the percentage of identify to the\n";
print OUT "most related viruses.\n\n";

print OUT "For the Interesting Reads section:\n";
print OUT "Sample name\n";
print OUT "Viral lineage, total number of reads assigned to this lineage, range of percentage of identity to the most related viruses\n\n";
print OUT "Description of the reads and the top hit BLAST information:\n";
print OUT "The fields are: Query Name, Query Length, Hit Name, Hit Length, Hit Description, Alignment Length, Percent Identity, Hit Start, Hit End, e value\n\n";
print OUT "Sequences in FASTA format\n\n";
print OUT $c;


print OUT "Viral sequence Summary:\n\n";
&generate_SampleDescription( $dir );
print OUT "End of Summary\n\n";
print OUT $c ;

print OUT "Ambiguous sequence Summary:\n\n";
&generate_SampleDescription_ambiguous( $dir );
print OUT "End of Summary\n\n";
print OUT $c ;

print OUT "\n\nSequence Report\n\n";
&generate_SequenceReport( $dir );
print OUT "End of Sequence Report\n\n";
print OUT $c ;

print OUT "\n\nTaxonomy Assignment:\n\n";
&generate_AssignmentSummary( $dir );
print OUT "End of Assignment\n\n";
print OUT $c ;

print OUT "\n\nViral Reads\n\n";
&generate_InterestingReads( $dir );
print OUT "End of Viral Reads\n\n";
print OUT $c ;

print OUT "\n\nAmbiguous Reads\n\n";
&generate_InterestingReads_ambiguous( $dir );
print OUT "End of Ambiguous Reads\n\n";
print "\n";

print OUT "# Finished\n";

exit;

############################################################################
sub generate_SampleDescription {
	my ($dir) = @_;

	print OUT $dir,"\n";
	printf OUT  "%15s\t%15s\t%40s\n", "ViralRead", "PercentIDrange", "IdentifiedVirus";

	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	foreach my $name (sort {$a cmp $b} @files) {
		if (!($name eq ".") && !($name eq "..")) {
		# name is either file name or sample name (directory)
			my $full_path = $dir."/".$name;
			if (-d $full_path) { # is a directory, sample directory
				# get total number of sequences in the sample
				my $SE1 = $full_path."/".$name."_SE1.RemoveAdapter.fastq";
				my $total_seq  = 0;
				if (-e $SE1) {
					my $line_number_SE1 = `wc -l $SE1`;
					$total_seq  = $line_number_SE1/4;
				}

				# print out report for this sample
				printf OUT "%30s\t%8d\n", $name,  $total_seq;
				my $Summary_file = $full_path."/".$name.".AssignmentSummary";
				if (-e $Summary_file) {
					open (IN, $Summary_file) or die "can not open file $Summary_file!\n";
					foreach (1..21) {
						<IN>;
					}
					while (<IN>) {
						if ($_ =~ /^\s*$/) { # empty line
							next;
						}
						elsif ($_ =~ /# Finished Assignment Summary/) { 
							next;
						}
						else {
							chomp $_;
							my $number_reads = 0;
							my $range = "";
							my @temp = split(/\t/, $_);
							my $range = pop @temp;
							my $info = pop @temp;
							my $virus_info = pop @temp;
							my $virus = "";
							if ($info =~ /total number of reads: (\d+)/) {
								$number_reads = $1;
							}	
							if ($virus_info =~ /hit does not have taxonomy entry/) {
								my @temp2 = split (",", $virus_info);
								$virus = shift @temp2;
							}
							else {
								my @temp2 = split(";", $virus_info);
								$virus = pop @temp2;
							}
							printf OUT "%15d\t%15s\t%40s\n", $number_reads, $range, $virus;
						}
					}
					close IN;
				}
				else {
					print OUT "$Summary_file does not exist!\n";
				}
			}
		}
	}
}

############################################################################
sub generate_SampleDescription_ambiguous {
	my ($dir) = @_;

	printf OUT  "%15s\t%15s\t%40s\n", "AmbiguousRead", "PercentIDrange", "VirusLineage";

	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	foreach my $name (sort {$a cmp $b} @files) {
		if (!($name eq ".") && !($name eq "..")) {
		# name is either file name or sample name (directory)
			my $full_path = $dir."/".$name;
			if (-d $full_path) { # is a directory, sample directory
				# get total number of sequences in the sample
				my $SE1 = $full_path."/".$name."_SE1.RemoveAdapter.fastq";
				my $total_seq  = 0;
				if (-e $SE1) {
					my $line_number_SE1 = `wc -l $SE1`;
					$total_seq  = $line_number_SE1/4;
				}

				# print out report for this sample
				printf OUT "%30s\t%8d\n", $name,  $total_seq;
				my $Summary_file = $full_path."/".$name.".AmbiguousReadsAssignmentSummary";
				if (-e $Summary_file) {
					open (IN, $Summary_file) or die "can not open file $Summary_file!\n";
					while (<IN>) {
						if ($_ =~ /^\s*$/) { # empty line
							next;
						}
						elsif ($_ =~ /# Finished Assignment Summary/) { 
							next;
						}
						else {
							chomp $_;
							my $number_reads = 0;
							my $range = "";
							my @temp = split(/\t/, $_);
							my $range = pop @temp;
							my $info = pop @temp;
							my $virus_info = pop @temp;
							my $virus = "";
							if ($info =~ /total number of reads: (\d+)/) {
								$number_reads = $1;
							}	
							if ($virus_info =~ /hit does not have taxonomy entry/) {
								my @temp2 = split (",", $virus_info);
								$virus = shift @temp2;
							}
							else {
								my @temp2 = split(";", $virus_info);
								$virus = pop @temp2;
							}
							printf OUT "%15d\t%15s\t%40s\n", $number_reads, $range, $virus;
						}
					}
					close IN;
				}
				else {
					print OUT "$Summary_file does not exist!\n";
				}
			}
		}
	}
}
##########################################################################
sub generate_SequenceReport {
	my ( $dir ) = @_;

	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	# use the same format as in the script generating assignment report
	printf OUT "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#totalPair", "#unjoined", "#joined", "%total", "#TotalAfterStitch", "#QC", "%Stitched", "#unique", "%Stitched", "sample" ;

	foreach my $name (sort {$a cmp $b} @files) {
		# name is either file name or sample name (directory)
		my $Summary_file = $dir."/".$name."/".$name.".AssignmentSummary";
		if (-e $Summary_file) {
			open (IN, $Summary_file) or die "can not open file $Summary_file!\n";
			while (my $line = <IN>) {
				if ($line =~ /#total/) {
					$line = <IN>;
					chomp $line;
					print OUT $line, "\t", $name, "\n";
					next;
				}		
			}
		}
	}

	print OUT $c2 ;
	printf OUT "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#TANTAN__good", "%unique", "#RM_good", "%unique","#RefGenomeFiltered", "%unique",  "#BN_VIRUSDB_Hit", "#BX_VIRUSDB_Hit", "Total_VIRUSDB_Hit", "%unique","VIRUSDB_Hit_Bac", "%VIRUSDB_Hit", "VIRUSDB_Hit_NON_Bac", "%VIRUSDB_Hit",  "sample";

	foreach my $name (sort {$a cmp $b} @files) {
		# name is either file name or sample name (directory)
		my $Summary_file = $dir."/".$name."/".$name.".AssignmentSummary";
		if (-e $Summary_file) {
			open (IN, $Summary_file) or die "can not open file $Summary_file!\n";
			while (my $line = <IN>) {
				if ($line =~ /#RM_good/) {
					$line = <IN>;
					chomp $line;
					print OUT $line, "\t", $name, "\n";
					next;
				}
			}
		}
	}
	print OUT $c2 ;

	printf OUT "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#BNinput",  "%VIRUSDB_Hit" , "#BXinput", "%VIRUSDB_Hit" , "#unassigned", "%VIRUSDB_Hit", "sample" ;

	foreach my $name (sort {$a cmp $b} @files) {
		# name is either file name or sample name (directory)
		my $Summary_file = $dir."/".$name."/".$name.".AssignmentSummary";
		if (-e $Summary_file) {
			open (IN, $Summary_file) or die "can not open file $Summary_file!\n";
			while (my $line = <IN>) {
				if ($line =~ /#BNinput/) {
					$line = <IN>;
					chomp $line;
					print OUT $line, "\t", $name, "\n";
					next;
				}
			}
		}
	}
	print OUT $c2 ;

}


#####################################################################
# Assignment Summary
sub generate_AssignmentSummary {
	my ( $dir ) = @_;

	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	foreach my $name (sort {$a cmp $b} @files) {
		# name is either file name or sample name (directory)
		my $full_path = $dir."/".$name;

		if (!($name eq ".") && !($name eq "..")) {
#		if (!($name =~ /\./)) {
			if (-d $full_path) { # is a directory
				my $Summary_file = $full_path."/".$name.".AssignmentSummary";
				if (-e $Summary_file) {
					open (IN, $Summary_file) or die "can not open file $Summary_file!\n";
					while (<IN>) {
						if ($_ =~ /# Finished Assignment Summary/) { 
							next;
						}

						print OUT $_;
					}
				}		
				print OUT $c2 ;
			}
		}
	}
}


############################################################################
sub count_num_of_seq () {
	my ($fastaFile) = @_;
	my $count = 0;
	
    open (FastaFile, $fastaFile) or die "Can't Open FASTA file: $fastaFile";
    while (my $line = <FastaFile>){
		if ($line =~ ">") {
			$count++;
		}
	}
	close FastaFile;

	return $count;
}

####################################################################################
# viral sequence
sub generate_InterestingReads {
	my ( $dir ) = @_;

	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	foreach my $name (sort {$a cmp $b} @files) {
		# name is either file name or sample name (directory)
		my $full_path = $dir."/".$name;

		if (!($name eq ".") && !($name eq "..")) {
#		if (!($name =~ /\./)) {
			if (-d $full_path) { # is a directory
				print OUT $name, "\n";
				my $tempF = $full_path."/".$name.".AssignmentReport";
				if ( -e $tempF ) {
					open (IN, $tempF) or die "can not open file $tempF!\n";
					foreach (1..20) {
						<IN>;
					}
					while (<IN>) {
						if ($_ =~ /# Finished Assignment Report/) { 
							next;
						}

						print OUT $_;
					}
					close IN;
				}
				else {
					print OUT "$name does not have .AssignmentReport file!\n";
				}		
				print OUT $c2;
			}
		}
	}
}

####################################################################################
# ambiguous sequence
sub generate_InterestingReads_ambiguous {
	my ( $dir ) = @_;

	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	foreach my $name (sort {$a cmp $b} @files) {
		# name is either file name or sample name (directory)
		my $full_path = $dir."/".$name;

		if (!($name eq ".") && !($name eq "..")) {
#		if (!($name =~ /\./)) {
			if (-d $full_path) { # is a directory
				print OUT $name, "\n";
				my $tempF = $full_path."/".$name.".AmbiguousReadsAssignmentReport";
				if ( -e $tempF ) {
					open (IN, $tempF) or die "can not open file $tempF!\n";
					while (<IN>) {
						if ($_ =~ /# Finished Assignment Report/) { 
							next;
						}

						print OUT $_;
					}
					close IN;
				}
				else {
					print OUT "$name does not have .AssignmentReport file!\n";
				}		
				print OUT $c2;
			}
		}
	}
}
