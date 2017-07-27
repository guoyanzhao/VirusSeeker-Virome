##---------------------------------------------------------------------------##
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
#color code
my $red = "\e[31m";
my $green = "\e[32m";
my $yellow = "\e[33m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $gray = "\e[37m";
my $normal = "\e[0m"; 


my $usage = "
This script will drive the execution of the pipeline 
for every sample in the given directory.

perl $0 <script to be run> <run folder> <ref genome><use checkpointing><step>
<script to be run> = name of the script to be run
<run folder> = full path of the folder holding files for this sequence run
               without the last \"/\"
<ref genome> = 1. Human genome
               2. Mouse
               3. Worm (C. elegans, C. briggsae)
               4. Worm (C. brenneri)
               5. Mouse lemur (Microcebus_murinus)
<use checkpointing> = 1. yes, 0. no

<step_number> = [1..36] run this pipeline step by step. (running the whole pipeline if step number is 0)
$red	[1] Remove Adapter
$green	[2] Stitch SE1 and SE2
$yellow	[3] Sequence quality control 
$purple	[4] Convert fastq format to fasta format	
$cyan	[5] Run CD-HIT to cluster similar sequences 
$red	[6] Tantan mask low complexity sequence
	[7] Sequence QC after Tantan
$cyan	[8]   Prepare for RepeatMasker
	[9]   Submit RepeatMasker job array
	[10]   Sequence quality control after RepeatMasker

$yellow	[11] Align to reference genome ( default Human, can be changed)
	[12] Extract unmapped reads
$gray	[13] Split for MegaBLAST against reference genome (default Human)
	[14]   Run MegaBLAST against reference genome (default Human)
	[15]   Parse MegaBLAST against reference genome output file (default Human)

$red	[16] Extract Ref genome filtered 
$gray	[17] Split files for BLASTn against VIRUSDB
	[18] Submit BLASTn VIRUSDB job array
	[19] Parse BLASTn VIRUSDB result

$red	[20] Split files for BLASTx VIRUSDB
	[21] Submit BLASTx VIRUSDB job array
	[22] Parse BLASTx VIRUSDB result

$yellow	[23] Pool sequence for mapping to Bacteria reference genome
	[24] Align to Bacteria reference genome 
	[25] Extract unmapped reads, convert to fasta format

$purple	[26] Split for MegaBLAST against NT 
	[27] Run MegaBLAST 
	[28] Parse MegaBLAST output file 

$yellow	[29]  Split for BLASTN NT
	[30]  Submit BLASTN job array
	[31]  Parse BlastN result

$purple	[32]  Pool and split files for BlastX
	[33]  Submit BlastX job array
	[34]  Parse BlastX result

$cyan	[35]  Generate assignment report for the sample
	[36]  Generate assignment summary for the sample

$red	[37]  Generate phage report
	[38]  Generate phage summary

$normal

";

die $usage unless scalar @ARGV == 5;
my ( $script, $dir, $ref_genome, $use_checkpoint, $step ) = @ARGV;

# change parameters based on which machine this script is running on
my $HOME = $ENV{HOME};

my $run_script_path = `dirname $0`;
chomp $run_script_path;
$run_script_path = "/usr/bin/perl ".$run_script_path."/";

#print "run scrip path is $run_script_path\n\n";

opendir(DH, $dir) or die "Can not open dir $dir!\n";
foreach my $name (readdir DH) {
	if (!($name eq ".") && !($name eq "..")) {
		my $full_path = $dir."/".$name;
		if (-d $full_path) { # is a directory
			my $com = $run_script_path.$script." $full_path $ref_genome $use_checkpoint  $step \n";
			print $com, "\n";
			system( $com );
		}
	}
}

exit;

