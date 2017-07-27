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
This script will read corresponding files in the given directory and 
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
my $outFile = $dir."/Analysis_Report_Phage_".$run_name.".txt";
open (OUT, ">$outFile") or die "can not open file $outFile!\n";

my ($wkday,$month,$day,$time,$year) = split(/\s+/, localtime);
print OUT "VirusSeeker-VirusOnly V${version}; Processing date: $day-$month-$year\n";
my $c = "**************************************************************************\n";
my $c2 = "#########################################################################\n\n";
###################################################
print OUT "How to read this file:\n\n";

print OUT "For the Summary section:\n";
print OUT "column 1: sample name\n";
print OUT "column 2: total number of reads obtained for this sample (pairs)\n\n";

print OUT "If there is any viral sequence identified in this sample, it will show up under the information \n";
print OUT "of this sample. There are 3 columns to describe the viral sequences identified in this sample:\n";
print OUT "column 1: number of viral sequences (could be read or contig)\n";
print OUT "column 2: range of percentage identity to blast hits. Some times one sequence can hit multiple \n";
print OUT "sequences in the database with the same e value. Only the first best hit percent identity is taken.\n";
print OUT "column 3: name of the virus (the last rank in the taxonomy lineage of the virus)\n\n";
 
print OUT "Summary:\n\n";
&generate_SampleDescription( $dir );
print OUT "End of Summary\n\n";
print OUT $c ;

print OUT "\n\nPhage Reads\n\n";
&generate_InterestingReads_phage( $dir );
print OUT "End of Phage Reads\n\n";
print "\n";

print OUT "# Finished\n";

exit;

############################################################################
sub generate_SampleDescription {
	my ($dir) = @_;

	print OUT $dir,"\n";
	printf OUT  "%15s\t%15s\t%40s\n", "PhageRead", "PercentIDrange", "IdentifiedVirus";

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
				printf OUT "%s\t%8d\n", $name,  $total_seq;
				my $Summary_file = $full_path."/".$name.".PhageAssignmentSummary";
				if (-e $Summary_file) {
					open (IN, $Summary_file) or die "can not open file $Summary_file!\n";
					<IN>; # the directory path line
					while (<IN>) {
						if ($_ =~ /^\s*$/) { # empty line
							next;
						}
						elsif ($_ =~ /# Finished Assignment Summary/) { 
							next;
						}
						elsif ($_ =~ /#/) { 
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
							printf OUT "%15d\t%20s\t%40s\n", $number_reads, $range, $virus;
						}
					}
				}
				else {
					print OUT "$Summary_file does not exist!\n";
				}
			}
		}
	}
}

####################################################################################
# phage sequence
sub generate_InterestingReads_phage {
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
				my $tempF = $full_path."/".$name.".PhageAssignmentReport";
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
					print OUT "$name does not have .PhageAssignmentReport file!\n";
				}		
				print OUT $c2;
			}
		}
	}
}

