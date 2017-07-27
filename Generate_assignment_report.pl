################################################################################
##  Author:
##      Guoyan Zhao <gzhao@pathology.wustl.edu>
#******************************************************************************
#* Copyright (C) Washington University in St. Louis Developed by Guoyan Zhao,
#* Washington University in St. Louis, Department of Pathology and Immunology.
#*
###############################################################################

#!/usr/bin/perl
use strict;
use Switch;
#use feature qw(switch);
use Bio::SearchIO;

my $usage = '
This script will read corresponding files in the given director and 
generate a report. 

perl script <sample dir> 
<sample dir> = full path to the directory holding files for the given 
               library without the last "/"
               e.g. ~/tools/Illumina_VirusDiscoveryPipeline/data/MiSeq_run_1/I10_12310_Project404_CSF_Glaser_Encephalitis_V11T00919_TruSeq-AD007_CAGATC

';
die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;

my @temp = split("\/", $dir);
my $sample_name = pop @temp;
# print "lib is $sample_name\n";

####################################
# read in original sequences
my $fasta_file = $dir."/".$sample_name.".RemoveAdapter.stitched.prinseq.fasta"; #

my %seq = (); # read_ID => sequence
my %seq_desc = (); # read_ID => description
&read_FASTA_data($fasta_file, \%seq, \%seq_desc);

####################################
# output files 
my $out1 = $dir."/".$sample_name.".AssignmentReport";
open (OUT1, ">$out1") or die "can not open file $out1!\n";
my $out2 = $dir."/".$sample_name.".ViralReads_all.fa";
open (OUT2, ">$out2") or die "can not open file $out2!\n";

my $out3 = $dir."/".$sample_name.".AmbiguousReadsAssignmentReport";
open (OUT3, ">$out3") or die "can not open file $out3!\n";
my $out4 = $dir."/".$sample_name.".AmbiguousReads_all.fa";
open (OUT4, ">$out4") or die "can not open file $out4!\n";

#my $out5 = $dir."/".$sample_name.".PhageReads_all.fa";
#open (OUT5, ">$out5") or die "can not open file $out5!\n";


##############################################################################
# variable to keep all the blast output files
my @blast_files_MegaBlast_NT = (); # all MegaBlast output files
my @blast_files_blastn = (); # all blastn.out files
my @blast_files_blastx = (); # all blastx.out files
my @blast_files_blastx_VIRUSDB = (); # all blastx.out files

# all the viral read ID
my %viral_reads_MegaBlast_NT = ();
my %viral_reads_blastn = ();
my %viral_reads_blastx = ();
my %viral_reads_blastx_VIRUSDB = ();

# viral_lineage => number of reads assigned to this lineage in the sample
my %viral_lineage_to_num_reads = ();
my %viral_reads_blast_info =();    # readID => information about this read
my %viral_reads_in_lineage_MegaBlast_NT = ();    # lineage => [read ID]
my %viral_reads_in_lineage_blastn = ();    # lineage => [read ID]
my %viral_reads_in_lineage_blastx = ();    # lineage => [read ID]

# all the ambiguous viral read ID
my %ambiguous_reads_MegaBlast_NT = ();
my %ambiguous_reads_blastn = ();
my %ambiguous_reads_blastx = ();

# ambiguous reads viral lineage => number of reads assigned to this lineage in the sample
my %ambiguous_reads_lineage_to_num_reads = ();
my %ambiguous_reads_blast_info =();    # readID => information about this read
my %ambiguous_reads_in_lineage_MegaBlast_NT = ();    # lineage => [read ID]
my %ambiguous_reads_in_lineage_blastn = ();    # lineage => [read ID]
my %ambiguous_reads_in_lineage_blastx = ();    # lineage => [read ID]

my @unassigned_reads = ();
my @phage_reads = ();

##################################
# generate unassigned reads file
my $BLASTX_NR_dir = $dir."/".$sample_name."_BLASTX_NR";
my $unassigned_reads_file = $dir."/".$sample_name.".unassigned_reads.fa";
if (-f $unassigned_reads_file) {
	unlink $unassigned_reads_file;
}
my $com = "cat $BLASTX_NR_dir/*BXfiltered.fa >> $unassigned_reads_file";
system ($com);

# category => num of sequence assigned to this category by blastn
my %Assignment_blastn = ( 
	"Bacteria" => 0,
	"Fungi" => 0,
	"Homo" => 0,
	"Mus" => 0,
	"Phage" => 0,
	"Viruses" => 0,
	"other" => 0,
	"unassigned" => 0,
	"Ambiguous" => 0,
);

# category => num of sequence assigned to this category by MegaBLAST 
my %Assignment_MegaBlast_NT = ();
foreach my $key (keys %Assignment_blastn) {
	$Assignment_MegaBlast_NT{$key} = 0;
}
# category => num of sequence assigned to this category by blastx 
my %Assignment_blastx = ();
foreach my $key (keys %Assignment_blastn) {
	$Assignment_blastx{$key} = 0;
}

######################################################################################
# get statistics of the sample
# to get total raw number of reads
my $SE1 = $dir."/".$sample_name."_SE1.RemoveAdapter.fastq";
my $SE2 = $dir."/".$sample_name."_SE2.RemoveAdapter.fastq";
my $line_number_SE1 = `wc -l $SE1`;
my $SE1_total_NumReads_raw = $line_number_SE1/4;
my $line_number_SE2 = `wc -l $SE2`;
my $SE2_total_NumReads_raw = $line_number_SE2/4;

# to get number of reads after stiching
my $unstitched_SE1 = $dir."/".$sample_name.".RemoveAdapter.un1.fq";
my $unstitched_SE2 = $dir."/".$sample_name.".RemoveAdapter.un2.fq";
my $joined = $dir."/".$sample_name.".RemoveAdapter.join.fq";
$line_number_SE1 = `wc -l $unstitched_SE1`;
my $un1_NumReads = $line_number_SE1/4;
$line_number_SE2 = `wc -l $unstitched_SE2`;
my $un2_NumReads = $line_number_SE2/4;
my $line_number_join = `wc -l $joined`;
my $joined_NumReads = $line_number_join/4;
my $total_stiched =  $un1_NumReads + $un2_NumReads + $joined_NumReads;

# to get number of reads after quality control
my $file = $dir."/".$sample_name.".RemoveAdapter.stitched.prinseq.fastq";
my $line_number = `wc -l $file`;
my $NumReads_QC = $line_number/4;

# to get number of unique reads after CD-HIT
$file = $dir."/".$sample_name.".QCed.cdhit.fa";
my $unique = `grep ">" $file |wc -l`;

# to get number of good reads after TANTAN
$SE1 = $dir."/".$sample_name.".QCed.cdhit.tantan.goodSeq.fa";
my $TANTAN_goodSeq = `grep ">" $SE1 |wc -l`;

# to get number of filtered and low complexity reads after RepeatMasker
$SE1 = $dir."/".$sample_name.".RepeatMasker.badSeq.fa";
my $RepeatMasker_input_total = 0;
my $RepeatMasker_good = 0;
my $RepeatMasker_filtered = 0;
my $RepeatMasker_LowComplexity = 0;

my $temp = `grep "total seq = " $SE1`;
if ($temp =~ /total seq = (\d+)/) {
	$RepeatMasker_input_total = $1;
}
$temp = `grep "good seq = " $SE1`;
if ($temp =~ /good seq = (\d+)/) {
	$RepeatMasker_good = $1;
}
$temp = `grep "bad seq = " $SE1`;
if ($temp =~ /bad seq = (\d+)/) {
	$RepeatMasker_filtered = $1;
}
$temp = `grep "Repeat and Low complexicity seq = " $SE1`;
if ($temp =~ /Repeat and Low complexicity seq = (\d+)/) {
	$RepeatMasker_LowComplexity = $1;
}

# to get number of reads not mapped to reference genome
$SE1 = $dir."/".$sample_name.".RepeatMasker.goodSeq.RefGenome.unmapped.masked.fasta";
my $unmapped = `grep ">" $SE1 |wc -l`;

# to get number of reads mapped to reference genome 
my $mapped = $RepeatMasker_good - $unmapped;

# to get number of reads after MegaBlast against reference genome
$SE1 = $dir."/".$sample_name.".RefGenomeFiltered.fa";
my $RefGenomeFiltered = `grep ">" $SE1 |wc -l`;

# to get number of viral hit in BLASTN VIRUSDB step
$file = $dir."/".$sample_name.".BLASTN_VIRUSDB_HIT.fa";
my $BLASTN_VIRUSDB_Hit_number = `grep ">" $file |wc -l`;

# to get number of viral hit in BLASTX VIRUSDB step
$file = $dir."/".$sample_name.".BLASTX_VIRUSDB_HIT.fa";
my $BLASTX_VIRUSDB_Hit_number = `grep ">" $file |wc -l`;

# to get total number of viral hit in BLASTN VIRUSDB and BLASTX VIRUSDB step
$file = $dir."/".$sample_name."_VIRUSDB_HIT.fa";
my $VIRUSDB_Hit_number = `grep ">" $file |wc -l`;

# to get number of seq NOT mapped to Bacteria reference genome
$file = $dir."/".$sample_name."_VIRUSDB_HIT_Bacteria_unmapped.fasta";
my $VIRUSDB_Hit_NONE_Bacteria = `grep ">" $file |wc -l`;
my $VIRUSDB_Hit_Bacteria = $VIRUSDB_Hit_number - $VIRUSDB_Hit_NONE_Bacteria;

# to get number of reads after MegaBLAST against full NT
$file = $dir."/".$sample_name.".VIRUSDB_HIT_MegaBLAST_NT_Filtered.fa";
my $MegaBLAST_Filtered = `grep ">" $file |wc -l`;

# get number of sequence after BLASTN NT filter and classified at this step
$file = $dir."/".$sample_name."_BLASTN_NT_Filtered.fa";
my $num_BLASTN_Filtered = `grep ">" $file |wc -l`;
my $num_BLASTN_Classified = $MegaBLAST_Filtered - $num_BLASTN_Filtered;

# get number of sequence classified at BLASTX NR
$file = $dir."/".$sample_name.".unassigned_reads.fa";
my $num_unassigned = `grep ">" $file |wc -l`;
my $num_BLASTX_NR_Classified = $num_BLASTN_Filtered - $num_unassigned;

####################################################################################
# to obtain viral read information
# MegaBLAST NT step
my $MegaBLAST_dir = $dir."/".$sample_name."_MegaBLAST_NT";
# enter directory where blast results resides
opendir (DIR, $MegaBLAST_dir) or die "can not open dir $MegaBLAST_dir!\n";
foreach my $blast_file (readdir DIR) {
	if ($blast_file =~ /megablast\.parsed$/) {
		# print "blastn parsed file $blast_file\n";
		my $parsed = $MegaBLAST_dir."/".$blast_file;
		my $blast_out = $blast_file;
		$blast_out =~ s/\.megablast\.parsed/\.megablast\.out/;
		$blast_out = $MegaBLAST_dir."/".$blast_out;
		push @blast_files_MegaBlast_NT, $blast_out;

		&collect_information($parsed, \%Assignment_MegaBlast_NT, \%viral_reads_MegaBlast_NT, \%viral_reads_in_lineage_MegaBlast_NT, \%viral_lineage_to_num_reads, \%ambiguous_reads_MegaBlast_NT, \%ambiguous_reads_in_lineage_MegaBlast_NT, \%ambiguous_reads_lineage_to_num_reads, \@unassigned_reads, \%viral_reads_blast_info, \%ambiguous_reads_blast_info, \@phage_reads );
	}
}
closedir DIR;

# BLASTN step
my $BLASTN_dir = $dir."/".$sample_name."_BLASTN_NT";
# enter directory where blastn results resides
opendir (BNDIR, $BLASTN_dir) or die "can not open dir $BLASTN_dir!\n";
foreach my $blast_file (readdir BNDIR) {
	if ($blast_file =~ /blastn\.parsed$/) {
		# print "blastn parsed file $blast_file\n";
		my $blast_out = $blast_file;
		$blast_out =~ s/\.blastn\.parsed/\.blastn\.out/;
		$blast_out = $BLASTN_dir."/".$blast_out;
		push @blast_files_blastn, $blast_out;
		my $parsed = $BLASTN_dir."/".$blast_file;
		&collect_information($parsed, \%Assignment_blastn, \%viral_reads_blastn, \%viral_reads_in_lineage_blastn, \%viral_lineage_to_num_reads, \%ambiguous_reads_blastn, \%ambiguous_reads_in_lineage_blastn, \%ambiguous_reads_lineage_to_num_reads, \@unassigned_reads, \%viral_reads_blast_info, \%ambiguous_reads_blast_info,\@phage_reads );
	}
}
closedir BNDIR;

# BLASTX NR step
my $BLASTX_dir = $dir."/".$sample_name."_BLASTX_NR";
# enter directory where blastx results resides
opendir (BXDIR, $BLASTX_dir) or die "can not open dir $BLASTX_dir!\n";
foreach my $blast_file (readdir BXDIR) {
	if ($blast_file =~ /blastx\.parsed$/) {
		# print "blast parsed file $blast_file\n";
		my $blast_out = $blast_file;
		$blast_out =~ s/\.blastx\.parsed/\.blastx\.out/;
		$blast_out = $BLASTX_dir."/".$blast_out;
		push @blast_files_blastx, $blast_out;
		my $parsed = $BLASTX_dir."/".$blast_file;
		&collect_information($parsed, \%Assignment_blastx, \%viral_reads_blastx, \%viral_reads_in_lineage_blastx, \%viral_lineage_to_num_reads, \%ambiguous_reads_blastx, \%ambiguous_reads_in_lineage_blastx, \%ambiguous_reads_lineage_to_num_reads, \@unassigned_reads, \%viral_reads_blast_info, \%ambiguous_reads_blast_info,\@phage_reads );
	}
}
closedir BXDIR;


#####################################################################################
# print out report for this library
print OUT1 $dir, "\n";

printf OUT1 "%s\t%s\t%s\t%s\t%s\t", "#totalPair", "#unjoined", "#joined", "%total", "#TotalAfterStitch"; 
printf OUT1  "%s\t%s\t%s\t%s\n", "#QC", "%Stitched" , "#unique", "%Stitched";

printf OUT1 "%d\t%d\t%d\t%7.2f\t%d\t", $SE1_total_NumReads_raw, $un1_NumReads, $joined_NumReads, $joined_NumReads*100/$SE1_total_NumReads_raw, $total_stiched;
printf OUT1 "%7d\t%7.2f\t%7d\t%7.2f\n\n", $NumReads_QC, $NumReads_QC*100/$total_stiched, $unique, $unique*100/$total_stiched;   

printf OUT1 "%s\t%s\t%s\t%s\t%s\t%s\t", "#TANTAN_good", "%unique", "#RM_good", "%unique", "#RefGenomeFiltered", "%unique"; 
printf OUT1 "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "#BN_VIRUSDB_Hit", "#BX_VIRUSDB_Hit", "Total_VIRUSDB_Hit", "%unique", "VIRUSDB_Hit_Bac", "%VIRUSDB_Hit", "VIRUSDB_Hit_NON_Bac", "%VIRUSDB_Hit" ; 
printf OUT1 "%7d\t%7.2f\t%7d\t%7.2f\t%7d\t%7.2f\t", $TANTAN_goodSeq, $TANTAN_goodSeq*100/$unique, $RepeatMasker_good, $RepeatMasker_good*100/$unique, $RefGenomeFiltered, $RefGenomeFiltered*100/$unique; 

printf OUT1 "%d\t%7d\t%7d\t%7.2f\t%7d\t%7.2f\t%7d\t%7.2f\n\n", $BLASTN_VIRUSDB_Hit_number,  $BLASTX_VIRUSDB_Hit_number, $VIRUSDB_Hit_number, $VIRUSDB_Hit_number*100/$unique, $VIRUSDB_Hit_Bacteria, $VIRUSDB_Hit_Bacteria*100/$VIRUSDB_Hit_number,  $VIRUSDB_Hit_NONE_Bacteria, $VIRUSDB_Hit_NONE_Bacteria*100/$VIRUSDB_Hit_number;

printf OUT1 "%s\t%s\t%s\t%s\t%s\t%s\n", "#BNinput", "%VIRUSDB_Hit", "BXinput", "%VIRUSDB_Hit", "Unassigned", "%VIRUSDB_Hit";

printf OUT1 "%d\t%7.2f\t%7d\t%7.2f\t%7d\t%7.2f\n\n",$MegaBLAST_Filtered, $MegaBLAST_Filtered*100/$VIRUSDB_Hit_number, $num_BLASTN_Filtered, $num_BLASTN_Filtered*100/$VIRUSDB_Hit_number, $num_unassigned, $num_unassigned*100/$VIRUSDB_Hit_number;

########################################################################################
# assignment in different categories
printf OUT1 "%12s\t%7s\t%9s\t%7s\t%7s\t%7s\n", "category", "total", "MegaBLAST", "BN_NT", "BX_NR" ;
foreach my $key (sort {$a cmp $b } keys %Assignment_blastx) {
	printf OUT1 "%12s\t%7d\t%9d\t%7d\t%7d\n", $key, $Assignment_blastn{$key}+$Assignment_blastx{$key}+$Assignment_MegaBlast_NT{$key},  $Assignment_MegaBlast_NT{$key}, $Assignment_blastn{$key} , $Assignment_blastx{$key};
}

##############################################################################
# viral lineage, number of viral reads in this lineage, BLAST alginment statistics.

my $c1 = "############################################################\n\n";
print OUT1 $c1;

foreach my $lineage (sort {$viral_lineage_to_num_reads{$a} <=> $viral_lineage_to_num_reads{$b}} keys %viral_lineage_to_num_reads) {
	print OUT1 $lineage, "\ttotal number of reads: ", $viral_lineage_to_num_reads{$lineage}, "\n\n";
	print OUT1 "QueryName\tQuerylength\t         HitName       \tHitLen\t                             HitDesc                       \tAlnLen\t%ID\tHitStart\tHitEnd\te\n";

	if (defined $viral_reads_in_lineage_MegaBlast_NT{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_MegaBlast_NT{$lineage}}) {
			print OUT1 "reads from MegaBlast NT:\n";
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_MegaBlast_NT{$lineage}}) {
				print OUT1 $viral_reads_blast_info{$read}, "\n";
			}
		}
		print OUT1 "\n";
	}

	if (defined $viral_reads_in_lineage_blastn{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_blastn{$lineage}}) {
			print OUT1 "reads from blastn:\n";
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_blastn{$lineage}}) {
				print OUT1 $viral_reads_blast_info{$read}, "\n";
			}
		}
		print OUT1 "\n";
	}

	if (defined $viral_reads_in_lineage_blastx{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_blastx{$lineage}}) {
			print OUT1 "reads from blastx:\n";
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_blastx{$lineage}}) {
				print OUT1 $viral_reads_blast_info{$read}, "\n";
			}
		}
		print OUT1 "\n";
	}

	# print out viral sequences in this taxonomy lineage
	if (defined $viral_reads_in_lineage_MegaBlast_NT{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_MegaBlast_NT{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_MegaBlast_NT{$lineage}}) {
				print OUT1 ">$read $seq_desc{$read}\n";
				print OUT1 $seq{$read}, "\n";
			}
		}
	}

	if (defined $viral_reads_in_lineage_blastn{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_blastn{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_blastn{$lineage}}) {
				print OUT1 ">$read $seq_desc{$read}\n";
				print OUT1 $seq{$read}, "\n";
			}
		}
	}

	if (defined $viral_reads_in_lineage_blastx{$lineage}) {
		if (scalar @{$viral_reads_in_lineage_blastx{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$viral_reads_in_lineage_blastx{$lineage}}) {
				print OUT1 ">$read $seq_desc{$read}\n";
				print OUT1 $seq{$read}, "\n";
			}
		}
	}

	print OUT1 $c1;
}

##########################################################################
# generate ambiguous reads report
print OUT3 $c1;
foreach my $lineage (sort {$ambiguous_reads_lineage_to_num_reads{$a} <=> $ambiguous_reads_lineage_to_num_reads{$b}} keys %ambiguous_reads_lineage_to_num_reads) {
	print OUT3 $lineage, "\ttotal number of reads: ", $ambiguous_reads_lineage_to_num_reads{$lineage}, "\n\n";
	print OUT3 "QueryName\tQuerylength\t         HitName       \tHitLen\t                             HitDesc                       \tAlnLen\t%ID\tHitStart\tHitEnd\te\n";

	if (defined $ambiguous_reads_in_lineage_MegaBlast_NT{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_MegaBlast_NT{$lineage}}) {
			print OUT3 "reads from MegaBlast NT:\n";
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_MegaBlast_NT{$lineage}}) {
				print OUT3 $ambiguous_reads_blast_info{$read}, "\n";
			}
		}
	}

	if (defined $ambiguous_reads_in_lineage_blastn{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_blastn{$lineage}}) {
			print OUT3 "reads from blastn:\n";
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_blastn{$lineage}}) {
				print OUT3 $ambiguous_reads_blast_info{$read}, "\n";
			}
		}
		print OUT3 "\n";
	}

	if (defined $ambiguous_reads_in_lineage_blastx{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_blastx{$lineage}}) {
			print OUT3 "reads from blastx:\n";
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_blastx{$lineage}}) {
				print OUT3 $ambiguous_reads_blast_info{$read}, "\n";
			}
		}
		print OUT3 "\n";
	}

	##########################################################################
	# print out sequences
	if (defined $ambiguous_reads_in_lineage_MegaBlast_NT{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_MegaBlast_NT{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_MegaBlast_NT{$lineage}}) {
				print OUT3 ">$read $seq_desc{$read}\n";
				print OUT3 $seq{$read}, "\n";
			}
		}
	}

	if (defined $ambiguous_reads_in_lineage_blastn{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_blastn{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_blastn{$lineage}}) {
				print OUT3 ">$read $seq_desc{$read}\n";
				print OUT3 $seq{$read}, "\n";
			}
		}
	}

	if (defined $ambiguous_reads_in_lineage_blastx{$lineage}) {
		if (scalar @{$ambiguous_reads_in_lineage_blastx{$lineage}}) {
			foreach my $read (sort {$a cmp $b} @{$ambiguous_reads_in_lineage_blastx{$lineage}}) {
				print OUT3 ">$read $seq_desc{$read}\n";
				print OUT3 $seq{$read}, "\n";
			}
		}
	}

	print OUT3 $c1;
}


###############################################################################
# get all the viral reads and put the sequence into output file:
foreach my $lineage (keys %viral_lineage_to_num_reads) {
	foreach my $read (@{$viral_reads_in_lineage_MegaBlast_NT{$lineage}}) {
		print OUT2 ">$read $seq_desc{$read}\n";
		print OUT2 $seq{$read}, "\n";
	}
	foreach my $read (@{$viral_reads_in_lineage_blastn{$lineage}}) {
		print OUT2 ">$read $seq_desc{$read}\n";
		print OUT2 $seq{$read}, "\n";
	}
	foreach my $read (@{$viral_reads_in_lineage_blastx{$lineage}}) {
		print OUT2 ">$read $seq_desc{$read}\n";
		print OUT2 $seq{$read}, "\n";
	}
}	
	
###############################################################################
# get all the ambiguous reads and put the sequence into output file:
foreach my $lineage (keys %ambiguous_reads_lineage_to_num_reads) {
	foreach my $read (@{$ambiguous_reads_in_lineage_MegaBlast_NT{$lineage}}) {
		print OUT4 ">$read $seq_desc{$read}\n";
		print OUT4 $seq{$read}, "\n";
	}
	foreach my $read (@{$ambiguous_reads_in_lineage_blastn{$lineage}}) {
		print OUT4 ">$read $seq_desc{$read}\n";
		print OUT4 $seq{$read}, "\n";
	}
	foreach my $read (@{$ambiguous_reads_in_lineage_blastx{$lineage}}) {
		print OUT4 ">$read $seq_desc{$read}\n";
		print OUT4 $seq{$read}, "\n";
	}
}

#####################
# output phage sequences
#foreach my $read (@phage_reads) {
#	print OUT5 ">$read $seq_desc{$read}\n";
#	print OUT5 $seq{$read}, "\n";
#}	
print OUT1 "# Finished Assignment Report\n";

exit;

#####################################################################################
# collecte information from given BLAST parsed output file
sub collect_information {
	################################## modified 04/2014
	my ($blast_parsed_file, $category_hash_ref, $viral_reads_hash_ref, $viral_reads_in_lineage_hash_ref, $viral_lineage_to_num_reads_hash_ref,  $ambiguous_reads_hash_ref, $ambiguous_reads_in_lineage_hash_ref, $ambiguous_reads_lineage_to_num_reads_hash_ref, $unassigned_reads_arr_ref, $viral_reads_blast_info_ref, $ambiguous_reads_blast_info_ref, $phage_reads_arr_ref ) = @_;

	##################################
	open (IN, $blast_parsed_file) or die "can not open file $blast_parsed_file!\n";
	while (<IN>) {
		if ($_ =~ /^#/) { # skip comment line
			next;
		}
		chomp;
		my @info = split("\t", $_);
		my $read_ID = shift @info;
		my $length = shift @info;
		my $category = shift @info;
		my $lineage = shift @info;
		my $desc = join ("\t", @info);

#		my ($read_ID, $length, $category, $lineage, $hit_name, $hit_length, $hit_desc, $hsp_len, $e_value) = split("\t", $_);
#		print "readID = \\$read_ID\\, length = \\$length\\, category = \\$category\\, lineage = \\$lineage\\, desc = \\$desc\\\n";
		switch ($category ) {
			case ("Bacteria") { $category_hash_ref->{"Bacteria"}++;	}
			case ("Fungi") { $category_hash_ref->{"Fungi"}++ ;}
			case ("Homo") { $category_hash_ref->{"Homo"}++ ;}
			case ("Mus") { $category_hash_ref->{"Mus"}++ ;}
			case ("Phage") {$category_hash_ref->{"Phage"}++; }
			case ("Viruses") { $category_hash_ref->{"Viruses"}++; }
			case ("other") {$category_hash_ref->{"other"}++; }
			case ("unassigned") {$category_hash_ref->{"unassigned"}++; }
			case ("Ambiguous") {$category_hash_ref->{"Ambiguous"}++; } 
		}

		$desc = $read_ID."\t".$length."\t".$desc; 
		if ($category eq "Viruses") {
			$viral_reads_hash_ref->{$read_ID} = 1;
			$viral_reads_blast_info_ref->{$read_ID} = $desc;
			if (!(defined $viral_reads_in_lineage_hash_ref->{$lineage})) {
				$viral_reads_in_lineage_hash_ref->{$lineage} = [$read_ID];
			}
			else {
				push @{$viral_reads_in_lineage_hash_ref->{$lineage}}, $read_ID;
			}
		
			if (defined $viral_lineage_to_num_reads_hash_ref->{$lineage}) {
				$viral_lineage_to_num_reads_hash_ref->{$lineage}++;
			}
			else {
				$viral_lineage_to_num_reads_hash_ref->{$lineage} = 1;
			}
		################################## modified 2014.04
		}elsif ($category eq "Ambiguous"){
			$ambiguous_reads_hash_ref->{$read_ID} = 1;
			$ambiguous_reads_blast_info_ref->{$read_ID} = $desc;
			if (!(defined $ambiguous_reads_in_lineage_hash_ref->{$lineage})) {
				$ambiguous_reads_in_lineage_hash_ref->{$lineage} = [$read_ID];
			}
			else {
				push @{$ambiguous_reads_in_lineage_hash_ref->{$lineage}}, $read_ID;
			}
		
			if (defined $ambiguous_reads_lineage_to_num_reads_hash_ref->{$lineage}) {
				$ambiguous_reads_lineage_to_num_reads_hash_ref->{$lineage}++;
			}
			else {
				$ambiguous_reads_lineage_to_num_reads_hash_ref->{$lineage} = 1;
			}
		#############################################################
		}elsif ($category eq "Phage"){		
			push @{$phage_reads_arr_ref}, $read_ID;
		##################################
		}elsif ($category eq "unassigned") {
			push @{$unassigned_reads_arr_ref}, $read_ID;
		}
	}
	close IN;
}


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

