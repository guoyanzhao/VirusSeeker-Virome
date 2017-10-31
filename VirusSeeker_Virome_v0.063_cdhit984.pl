#!/usr/bin/perl
##################################################################################
# Contact: Guoyan Zhao
# Email: gzhao@pathology.wustl.edu
# Changes over previous version:
# 1) At step [22] Pool sequence for mapping to Bacteria reference genome
# Pool only viral hits instead of viral hits and all the BLAST VIRUSDB filtered
# reads.
# 2) remove MegaBLAST against Bacteria + VirusDB, add MegaBLAST against NT as in
# VirusSeeker (VirusOnly) pipeline
# 3) modified number of files to split into for Repeatmasker, BLASTN and BLASTX 
# VirusOnly-database, MegaBLAST NT, BLAST NT, NR  steps.
# changed number of sequence files at each step to decrease number of job slots
# requested. 
# added checkpointing at all necessary steps.
##################################################################################

use strict;
use warnings;

my $version = 0.063;
my $cdhit_cutoff = 0.984;

#color code
my $red = "\e[31m";
my $green = "\e[32m";
my $yellow = "\e[33m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $gray = "\e[37m";
my $normal = "\e[0m"; 

#usage information
(my $usage = <<OUT) =~ s/\t+//g;
This script will run the Metagenomic pipeline use Slurm Workload Manager.

Pipeline version: $version
$yellow		Usage: perl $0 <sample_folder><ref genome><use checkpointing> <step_number> $normal

<sample_folder> = full path of the folder holding files for the sample
<ref genome> = 1. Human genome
               2. Mouse
               3. Worm (C. elegans, C. briggsae)
               4. Worm (C. brenneri)
               5. Mouse lemur (Microcebus_murinus)
<use checkpointing> = 1. yes, 0. no

<step_number> = [1..38] run this pipeline step by step. (running the whole pipeline if step number is 0)
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

$purple [26] Split for MegaBLAST against NT 
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
OUT

die $usage unless scalar @ARGV == 4;
my ($sample_dir, $ref_genome_choice, $use_checkpoint,  $step_number) = @ARGV;
die $usage unless (($step_number >= 0)&&($step_number <= 38)) ;


#####################################################################################
# get name of the sample and path to the data
if ($sample_dir =~/(.+)\/$/) {
	$sample_dir = $1;
}
my $sample_name = (split(/\//,$sample_dir))[-1];
print $sample_name,"\n\n";

#####################################################################################
# path and name of databases
my $db_BN = "/scratch/ref/dwlab/nt_20160802/nt";
my $db_BX = "/scratch/ref/dwlab/nr_20160802/nr";
# Virus database created by extracting all viral sequences from NR and then
# CD-HIT clustered with 98% ID
my $NR_VIRUS = "/scratch/ref/dwlab/VirusDBNR_20160802/VirusDBNR_20160802_ID98.fa";
my $NT_VIRUS = "/scratch/ref/dwlab/VirusDBNT_20160802/VirusDBNT_20160802_ID98.fa";

####################################################################################
# software path
my $cdhit = "cd-hit-est";
# -c sequence identity threshold, default 0.9. This is the default cd-hit's "global sequence 
# identity" calculated as: number of identical amino acids in alignment divided by the full 
# length of the shorter sequence
# -n word_length, default 5, see user's guide for choosing it
# -G use global sequence identity, default 1. If set to 0, then use local sequence identity, 
# calculated as : number of identical amino acids in alignment divided by the length of the 
# alignment NOTE!!! don't use -G 0 unless you use alignment coverage controls see options -aL, -AL, -aS, -AS
# -aS alignment coverage for the shorter sequence, default 0.0. If set to 0.9, the alignment 
# must covers 90% of the sequence
# -g 1 or 0, default 0. By cd-hit's default algorithm, a sequence is clustered to the first 
# cluster that meet the threshold (fast mode). If set to 1, the program will cluster it into 
# the most similar cluster that meet the threshold (accurate but slow mode)
# -r 1 or 0, default 0, if set to 1, comparing both strand (++, +-)
# -M max available memory (Mbyte), default 400. Here requested for 20 GB
# -d length of description in .clstr file, default 20 if set to 0, it takes the fasta defline and stops at first space
# -T number of threads, default 1; with 0, all CPUs will be used
#my $cdhit_parameter = " -c 1 -n 8 -G 0 -aS 1 -g 1 -r 1 -M 20480 -d 0 -T 6 ";
# 99% ID over 99% overlap
# $cdhit_cutoff is the identify and overlap cutoff
my $cdhit_parameter = " -c $cdhit_cutoff -n 8 -G 0 -aS $cdhit_cutoff -g 1 -r 1 -M 20480 -d 0 -T 8 ";
my $repeat_masker = "RepeatMasker";
my $blastn = "blastn";
my $blastx = "blastx";
my $prinseq = "/opt/apps/prinseq/0.20.4/prinseq-lite.pl";
my $Tantan = "tantan ";
my $bwa = "bwa";
my $samtools = "samtools";
my $bamToFastq = "bamToFastq";

###################################################################################
# directory suffix constants
my $MegaBLAST_dir_host = $sample_dir."/".$sample_name."_MegaBLAST_HOST";
my $repeatmasker_dir = $sample_dir."/".$sample_name."_RepeatMasker";
my $BLASTN_VIRUSDB_DIR_SUFFIX = "_BLASTN_VIRUSDB"; 	# blastn against virus-only (NT extracted) database
my $BLASTX_VIRUSDB_DIR_SUFFIX = "_BLASTX_VIRUSDB"; # blastx against virus-only (NR extracted) database

my $MegaBLAST_dir_NT = $sample_dir."/".$sample_name."_MegaBLAST_NT";
my $BLASTN_NT_dir = $sample_dir."/".$sample_name."_BLASTN_NT";
my $BLASTX_NR_dir =$sample_dir."/".$sample_name."_BLASTX_NR";

my $Bacteria_genome = "/scratch/ref/dwlab/Bacteria_ref/Bacteria_ref_genome.fna";

# reference genome taxonomy classification and database location.
# Has to change $refrence_genome_taxonomy and $reference_genome based on the data 
# being analyzed!!!
my $refrence_genome_taxonomy = "";
my $reference_genome = "";
my $last_jobid ="";

if ($ref_genome_choice == 1) {
	$refrence_genome_taxonomy = "Homo"; # use Bacteria, Homo, Phage, Fungi, Mus, other

	# path to the reference genome 
	$reference_genome = "/scratch/ref/dwlab/HumanGenome_BWA/human_g1k_v37.fasta";
}
elsif ($ref_genome_choice == 2) {
	$refrence_genome_taxonomy = "Mus"; # use Bacteria, Homo, Phage, Fungi, Mus, other
	$reference_genome = "/scratch/ref/dwlab/MouseGenome/MouseGenome_UCSC_mm10_2012_12_12.fa";
}
elsif ($ref_genome_choice == 3) {
	$refrence_genome_taxonomy = "other"; # use bacteria, homo, phage, fungi, mus, other
	$reference_genome = "/scratch/ref/dwlab/CelegansCbriggsae/Celegans_WS220_Cbriggsae_WS234.fa";
}
elsif ($ref_genome_choice == 4) {
	$refrence_genome_taxonomy = "other"; # use bacteria, homo, phage, fungi, mus, other
	$reference_genome = "/scratch/ref/dwlab/Cbrenneri/C_brenneri_genome.fa";
}
elsif ($ref_genome_choice == 5) {
	$refrence_genome_taxonomy = "other"; # use Bacteria, Homo, Phage, Fungi, Mus, other
	$reference_genome = "/srv/cgs/data/gzhao/MouseLemur_genome/Microcebus_murinus.micMur1.71.dna_sm.toplevel.fa";
}

####################################################################################
# Here are the parameters that can be adjusted according to computer cluster 
# configuration.
# Use this as max number of reads per file as input for assembler. 
####################################################################################
# To run jobs faster, split large fasta files to small ones. Split to specific number 
# of files instead of specific sequences in each small file, because the number of job 
# array cannot be determined if spliting to specific number of sequences in each file. 
# Job number is required by qsub ${SLURM_ARRAY_TASK_ID} the minimum size of each file is 4kb

# Repeatmasker
# After Brian modified RepeatMasker code, the execution time dropped from 20 min to 2 min 
# for 2400 sequences. For 24000 sequences, should take about 2 min.
# minimum number of sequences per file when spliting into smaller files
my $num_seq_per_file_Repeatmasker = 36000;

# the number of job files for RepeatMasker. The most optimal number would be an 
# estimated number that makes the number of sequences in each file close to 2400. 
#my $file_number_of_RepeatMasker = 1000; 
my $file_number_of_RepeatMasker = 30; 

# MegaBLAST
# 20,000 reads use ~15 min in MegaBLAST against reference genome. If too many reads in 
# one file the memory usage is too big. Job may be killed before it finishes. 
my $num_seq_per_file_MegaBLAST_host = 40000;
# If 1 sample is sequenced in 1 MiSeq run, 400 file would be max number of files
# based on percentage of reads left for analysis at each step from MiSeq run 1 data. 
# This number should be bigger to make sure each file does not have too many reads.
#my $file_number_of_MegaBLAST = 3000;
my $file_number_of_MegaBLAST_host = 25;

# BLASTN virus only database. 10,000 sequence takes ~ 10 min
# Specify the number of sequences to be in each file 
my $num_seq_per_file_BLASTN_VIRUSDB = 10000;
# Specify the total number of files to be splited into
my $file_number_of_BLASTN_VIRUSDB = 100; 

# BLASTX virus only database.
# Specify the number of sequences to be in each file 
my $num_seq_per_file_BLASTX_VIRUSDB = 1000;
# Specify the total number of files to be splited into
my $file_number_of_BLASTX_VIRUSDB = 800; 

# MegaBLAST NT
# 2,000 reads use ~40 min in MegaBLAST against reference genome. If too many reads in 
# one file the memory usage is too big. Job may be killed before it finishes. 
my $num_seq_per_file_MegaBLAST_NT = 1000;
# If 1 sample is sequenced in 1 MiSeq run, 400 file would be max number of files
# based on percentage of reads left for analysis at each step from MiSeq run 1 data. 
# This number should be bigger to make sure each file does not have too many reads.
my $file_number_of_MegaBLAST_NT = 50;

# BLASTN: 1000 reads takes ~ 6 hours. 
# minimum number of sequences per file
my $num_seq_per_file_BLASTN = 200;
# Specify the number of small fasta files to split to from a large file for BLASTN
# The most optimal number would be an estimated number that makes the number of 
# sequences in each file close to 200.
my $file_number_of_BLASTN = 150; 
#my $file_number_of_BLASTN = 1; 

# BLASTX: 1000 reads takes ~ 15 hours.
# Specify the number of sequences to be in each file for BLASTX NR
my $num_seq_per_file_BLASTX = 80;
# Specify the number of small fasta files to split to from a large file for BLASTX
my $file_number_of_BLASTX = 200; 
#my $file_number_of_BLASTX = 200; 

####################################################################################
# Everything else below should be automated.
#my $HOME = $ENV{HOME};

# store job files here
my $job_files_dir = $sample_dir."/job_script";
if (! -d $job_files_dir) {
	`mkdir $job_files_dir`;
}

#create a folder to store SGE output and error files
#my $sge_files_dir = $sample_dir."/SGE_DIR";
#if (! -d $sge_files_dir) {
	#`mkdir $sge_files_dir`;
#}
#create a folder to store SGE output and error files
my $SLURM_files_dir = $sample_dir."/SLURM_DIR";
if (! -d $SLURM_files_dir) {
	`mkdir $SLURM_files_dir`;
}
#create a folder to store status files
my $status_log = $sample_dir."/status_log";
if (! -d $status_log) {
	`mkdir $status_log`;
}

my $run_script_path = `dirname $0`;
chomp $run_script_path;
$run_script_path = "perl ".$run_script_path."/";
#print "run script path is : $run_script_path \n";
 
# variables
my $current_job_file1 = "no"; #cannot be empty
my $current_job_file2 = "no"; #cannot be empty
my $hold_job_file1 = "";
my $hold_job_file2 = "";

####################################################################################
# run the whole pipeline
if ($step_number == 0) {
	# 1
	&remove_adapter();
	# 2 
	&stitching();
#        exit;
	# 3
	&quality_control();
	# 4
	&fastqtofasta(); 
	# 5
	&run_CDHIT(); 
	# 6
	&run_Tantan(); 
	# 7
	&seq_QC_TANTAN(); 
	# 8
	&prepare_for_RM();
	# 9
	&submit_job_array_RM();
	# 10
	&seq_RMQC();
	# 11
	&map_Host_Reference_Genome_BWA();
	# 12
	&Extract_Host_unmapped_reads();
	# 13
	&split_for_MegaBLAST_RefGenome();
	# 14
	&MegaBLAST_RefGenome();
	# 15
	&parse_MegaBLAST_RefGenome();
	# 16
	&extract_RefGenomeFiltered_reads();
	# 17
	&split_for_BLASTN_VIRUSDB();
	# 18
	&submit_job_array_BLASTN_VIRUSDB();
	# 19
	&parse_BLASTN_VIRUSDB();
	# 20
	&split_for_BLASTX_VIRUSDB();
	# 21
	&submit_job_array_BLASTX_VIRUSDB();
	# 22
	&parse_BLASTX_VIRUSDB();
	# 23
	&pool_seq_for_mapping_Bacteria();
	# 24
	&map_Bacteria_Genome();
	# 25
	&Extract_Bacteria_unmapped_reads();
	# 26
	&split_for_MegaBLAST_NT();
	# 27
	&submit_job_array_MegaBLAST_NT();
	# 28
	&parse_MegaBLAST_NT();
	# 29
	&pool_split_for_BlastN();
	# 30
	&submit_job_array_BLASTN();
	# 31
	&parse_BLASTN();
	# 32 
	&pool_split_for_BLASTX();
	# 33
	&submit_job_array_BLASTX();
	# 34
	&parse_BLASTX();
	# 35
	&generate_assignment_report();
	# 36
	&generate_assignment_summary();
	# 37
	&generate_phage_report();
	# 38
	&generate_phage_summary();

}elsif ($step_number == 1) {
	&remove_adapter();
}elsif ($step_number == 2) {
	&stitching();
}elsif ($step_number == 3) {
	&quality_control();
}elsif ($step_number == 4) {
	&fastqtofasta(); 
}elsif ($step_number == 5) {
	&run_CDHIT(); 
}elsif ($step_number == 6) {
	&run_Tantan(); 
}elsif ($step_number == 7) {
	&seq_QC_TANTAN(); 
}elsif ($step_number == 8){
	&prepare_for_RM();
}elsif ($step_number == 9) {
	&submit_job_array_RM();
}elsif ($step_number == 10) {
	&seq_RMQC();
}elsif ($step_number == 11) {
	&map_Host_Reference_Genome_BWA();
}elsif ($step_number == 12) {
	&Extract_Host_unmapped_reads();
}elsif ($step_number == 13) {
	&split_for_MegaBLAST_RefGenome();
}elsif ($step_number == 14) {
	&MegaBLAST_RefGenome();
}elsif ($step_number == 15) {
	&parse_MegaBLAST_RefGenome();
}elsif ($step_number == 16) {
	&extract_RefGenomeFiltered_reads();
}elsif ($step_number == 17) {
	&split_for_BLASTN_VIRUSDB();
}elsif ($step_number == 18){
	&submit_job_array_BLASTN_VIRUSDB();
}elsif ($step_number == 19){
	&parse_BLASTN_VIRUSDB();
}elsif ($step_number == 20){
	&split_for_BLASTX_VIRUSDB();
}elsif ($step_number == 21) {
	&submit_job_array_BLASTX_VIRUSDB();
}elsif ($step_number == 22) {
	&parse_BLASTX_VIRUSDB();
}elsif ($step_number == 23) {
	&pool_seq_for_mapping_Bacteria();
}elsif ($step_number == 24) {
	&map_Bacteria_Genome()
}elsif ($step_number == 25) {
	&Extract_Bacteria_unmapped_reads();
}elsif ($step_number == 26) {
	&split_for_MegaBLAST_NT();
}elsif ($step_number == 27) {
	&submit_job_array_MegaBLAST_NT();
}elsif ($step_number == 28) {
	&parse_MegaBLAST_NT();
}elsif ($step_number == 29) {
	&pool_split_for_BlastN();
}elsif ($step_number == 30) {
	&submit_job_array_BLASTN();
}elsif ($step_number == 31){
	&parse_BLASTN();
}elsif ($step_number == 32) {
	&pool_split_for_BLASTX();
}elsif ($step_number == 33) {
	&submit_job_array_BLASTX();
}elsif ($step_number == 34) {
	&parse_BLASTX();
}elsif ($step_number == 35){
	&generate_assignment_report();
}elsif ($step_number == 36){
	&generate_assignment_summary();
}elsif ($step_number == 37){
	&generate_phage_report();
}elsif ($step_number == 38) {
	&generate_phage_summary();
}else{
	die $usage;
}

exit;

###########################################################################
# Remove adapter sequences from reads
sub remove_adapter {
#	$hold_job_file1 = "no";
#	$hold_job_file2 = "no";
	$current_job_file1 = "j1_".$sample_name."_RemoveAdapter1.sh";
#	$current_job_file2 = "j1_".$sample_name."_RemoveAdapter2.sh";

	my $Adapter_file = $sample_dir."/".$sample_name."_adapter.txt";
	my @adaper_seqs = ();
	open(IN, "<$Adapter_file") or die "$Adapter_file does not exist!!!";
	while (<IN>) {
		if ($_ =~ /^\s$/) { # skip blank line
			next;
		}
		chomp;
		push @adaper_seqs, $_;
	}
	close IN;

	my $command_line_adapter_seq = "";
	foreach my $adaper (@adaper_seqs) {
		$command_line_adapter_seq .= " -b \"".$adaper."\" ";
	}
	 
	######################################################################
	open(RemoveAdapter, ">$job_files_dir/$current_job_file1") or die $!;
	print RemoveAdapter "#!/bin/bash\n";
	print RemoveAdapter "#SBATCH --array=1-2\n";
	print RemoveAdapter "#SBATCH --mem-per-cpu=20G\n";
	print RemoveAdapter "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out\n";
        print RemoveAdapter "module load cutadapt\n";
	print RemoveAdapter "IN=$sample_dir"."/".$sample_name."_SE\${SLURM_ARRAY_TASK_ID}.fastq.gz\n";
	print RemoveAdapter "OUTFILE=$sample_dir"."/".$sample_name."_SE\${SLURM_ARRAY_TASK_ID}.RemoveAdapter.fastq\n";
	print RemoveAdapter "REPORT=$sample_dir"."/".$sample_name."_SE\${SLURM_ARRAY_TASK_ID}.RemoveAdapter.report\n\n";

	# Check if RemoveAdapter finished successfully. If finished successfully, resubmitting 
	# this job will skip RemoveAdapter step and will not overwrite the results
	print RemoveAdapter "if [ ! -e $status_log/j1_RemoveAdapter\${SLURM_ARRAY_TASK_ID}_finished ]\n";
	print RemoveAdapter "then\n";
	print RemoveAdapter " 	cutadapt ".$command_line_adapter_seq." -n 5 -O 5  -o \${OUTFILE}  \${IN} > \${REPORT}   \n";
	print RemoveAdapter "	OUT=\$?\n"; 
	print RemoveAdapter '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print RemoveAdapter "	then\n";
	print RemoveAdapter "		echo \"Fatal Error trying to run cutadapt.\"  \n";
	print RemoveAdapter "		exit 1\n"; 
	# if prinseq finished succussfully, make a file named QC_finished
	print RemoveAdapter "	else\n";
	# RemoveAdapter finished succussfully, make a file named RemoveAdapter1_finished
	print RemoveAdapter "		touch $status_log/j1_RemoveAdapter\${SLURM_ARRAY_TASK_ID}_finished\n";
	print RemoveAdapter "	fi\n\n";
	print RemoveAdapter "fi\n\n";
	close RemoveAdapter;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

##################################################################################
# Stitch SE1 and SE2 reads together
sub stitching {
	$current_job_file1 = "j2_".$sample_name."_Stitch.sh";
	$current_job_file2 = $current_job_file1;
		
	######################################################################
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";         
	print STCH "#SBATCH --mem-per-cpu=10G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load ea-utils\n";
	print STCH "IN1=$sample_dir"."/".$sample_name."_SE1.RemoveAdapter.fastq\n";
	print STCH "IN2=$sample_dir"."/".$sample_name."_SE2.RemoveAdapter.fastq\n";
	print STCH "OUTFILE=$sample_dir"."/".$sample_name.".RemoveAdapter.%.fq\n";
	print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".Stitching_time.txt\n";
	print STCH "OVERLAPLENGTH=$sample_dir"."/".$sample_name.".Stitching_OVERLAPLENGTH.txt\n";
	my $stitching_report_file = $sample_dir."/".$sample_name.".stitching.report";
	my $full_sample_name=$sample_dir."/".$sample_name.".RemoveAdapter"; 

	# check whether stitching finished successfully. If finished successfully, 
	# resubmitting this job will skip stitching step and will not overwrite the results
	print STCH "if [ ! -e $status_log/j2_Stitching_finished ]\n";
	print STCH "then\n";
	print STCH "	date > \${TIMEFILE}\n\n";
	print STCH "	fastq-join -p 5 -m 10 -r \${OVERLAPLENGTH}  \${IN1} \${IN2} -o \${OUTFILE} > $stitching_report_file\n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error trying to run fastq-join.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n\n";

	# the fastq-join software produce 3 files: join, un1 and un2. 
	# Need to put all the reads into a single file
	print STCH "	".$run_script_path."Combine_joined_unjoin_reads.pl $full_sample_name\n";
	print STCH "	OUT2=\$?\n"; 
	print STCH '	if [ ${OUT2} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error trying to post process fastq-join files.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n\n";

	# Stitching finished succussfully, make the log file 
	print STCH "	touch $status_log/j2_Stitching_finished\n";
	print STCH "	date >> \${TIMEFILE}\n";
	print STCH "fi\n\n";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

###########################################################################
# Quality control: trim low quality nucleotide, remove low quality reads, remove low 
# complexity reads, remove poly A/T tail etc.
# starting from this step use $current_job_file and $hold_job_file

sub quality_control {
	$current_job_file1 = "j3_".$sample_name."_QC.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";         
	print STCH "#SBATCH --mem-per-cpu=10G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load prinseq\n";
	print STCH "IN=$sample_dir"."/".$sample_name.".RemoveAdapter.stitched.fastq\n";
	print STCH "OUTFILE=$sample_dir"."/".$sample_name.".RemoveAdapter.stitched.prinseq\n";
	print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".prinseq_time.txt\n";

	print STCH "if [ ! -e $status_log/j3_QC_finished ]\n";
	print STCH "then\n";
	print STCH "date > \${TIMEFILE}\n";
	# -phred64: Not required for Illumina 1.8+, Sanger, Roche/454, Ion Torrent, PacBio data
	# QC score cutoff 30
#	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -min_qual_score 20  -out_good \${OUTFILE} ";

	# variant of QC score cutoff 30, min_qual_score to 15
#	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -min_qual_score 15  -out_good \${OUTFILE} ";

	# variant of QC score cutoff 30, min_qual_score to 10
#	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -min_qual_score 10  -out_good \${OUTFILE} ";

	# variant of QC score cutoff 30, min_qual_score to 3
#	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -min_qual_score 3  -out_good \${OUTFILE} ";


	# variant of QC score cutoff 30,  min_qual_score to 0
#	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30   -out_good \${OUTFILE} ";


	# variation of QC score cutoff 25,  min_qual_score to 10
	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 25 -trim_qual_left 25 -trim_qual_right 25 -min_qual_score 10  -out_good \${OUTFILE} ";

	# variation of QC score cutoff 25,  min_qual_score to 3
#	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 25 -trim_qual_left 25 -trim_qual_right 25 -min_qual_score 3  -out_good \${OUTFILE} ";


	# variation of QC score cutoff 25,  min_qual_score to 0
#	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 25 -trim_qual_left 25 -trim_qual_right 25   -out_good \${OUTFILE} ";

	# QC score cutoff 20, min_qual_score 15
#	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 20 -trim_qual_left 20 -trim_qual_right 20 -min_qual_score 15  -out_good \${OUTFILE} ";


	# variation of QC score cutoff 20,  min_qual_score to 10
#	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 20 -trim_qual_left 20 -trim_qual_right 20 -min_qual_score 10  -out_good \${OUTFILE} ";

	# variation of QC score cutoff 20,  min_qual_score to 3
#	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 20 -trim_qual_left 20 -trim_qual_right 20 -min_qual_score 3  -out_good \${OUTFILE} ";

	# variation of QC score cutoff 20,  min_qual_score to 0
#	my $com = $prinseq." -fastq \${IN}  -no_qual_header  -min_len 50 -ns_max_p 4 -min_gc 10 -max_gc 90 -derep  14  -lc_method  dust  -lc_threshold  8 -trim_tail_left  5 -trim_tail_right  5 -trim_ns_left 1  -trim_ns_right  1 -min_qual_mean 20 -trim_qual_left 20 -trim_qual_right 20   -out_good \${OUTFILE} ";

	print STCH "	", $com, "\n"; 
	# check if prinseq finished successfully. If finished successfully, resubmitting 
	# this job will skip prinseq step and will not overwrite the results
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error trying to run prinseq.\"  \n";
	print STCH "		exit 1\n"; 
	# if prinseq finished succussfully, make a file named QC_finished
	print STCH "	else\n";
	print STCH "		date >> \${TIMEFILE}\n";
	print STCH '		echo "', $com, '"  >> ${TIMEFILE} ', "\n";
	print STCH "		touch $status_log/j3_QC_finished\n";
	print STCH "	fi\n";
	print STCH "fi\n\n";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#################################################
# convert fastq to fasta format
sub fastqtofasta{
	$current_job_file1 = "j4_".$sample_name."_FastqToFasta.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";         
	print STCH "#SBATCH --mem-per-cpu=20G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "IN=$sample_dir"."/".$sample_name.".RemoveAdapter.stitched.prinseq.fastq\n";
	print STCH "OUTFILE=$sample_dir"."/".$sample_name.".RemoveAdapter.stitched.prinseq.fasta\n";

 	print STCH "OUT=1\n";
	print STCH "if [ ! -e $status_log/j4_FastqToFasta_finished ]\n";
	print STCH "then\n";
	print STCH "	".$run_script_path."FASTQ_to_FASTA.pl \${IN} > \${OUTFILE}\n"; 
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error trying to run j4 FastqToFasta conversion.\"  \n";
	print STCH "		exit 1\n"; 
	# if finished succussfully, make a file named QC_finished
	print STCH "	else\n";
	print STCH "		touch $status_log/j4_FastqToFasta_finished\n";
	print STCH "	fi\n";
	print STCH "fi\n\n";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
 }

#######################################################################################
# This step will cluster sequences using CD-HIT to reduce redundancy.
sub run_CDHIT{
	# this is the job script
	my $job_script = $job_files_dir."/j5_".$sample_name."_CDHIT_script.sh";
	open(STCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print STCH "#!/bin/bash\n";         
	print STCH "CDHIT_IN=$sample_dir"."/".$sample_name.".RemoveAdapter.stitched.prinseq.fasta \n";
	print STCH "CDHIT_OUT=$sample_dir"."/".$sample_name.".QCed.cdhit.fa \n";
	print STCH "CDHIT_REPORT=".$sample_dir."/".$sample_name.".RemoveAdapter.stitched.prinseq.cdhitReport\n";
	my $cdhitReport = $sample_dir."/".$sample_name.".RemoveAdapter.stitched.prinseq.cdhitReport";

	print STCH "set -x\n";
	print STCH "if [ ! -e $status_log/j5_CDHIT_finished ]\n";
	print STCH "then\n";

	#if CD_HIT output file does not exist, run and check the completeness of output, 
	# -T 8 to request 8 CPU, use -pe shared 8  in qusub to request 8 CPU. -M is memory limit, 20G = 20480 M
	# the -M specified memory size should be in consistant with the "h_vmem" specified memory request for 
	# the job submission.
	print STCH '	if [ ! -e ${CDHIT_REPORT} ]',"\n"; # no report file
	print STCH "	then\n";
	print STCH "		$cdhit -i \${CDHIT_IN} -o \${CDHIT_OUT} ".$cdhit_parameter." > \${CDHIT_REPORT} ", "\n";
	print STCH "		OUT=\$?\n"; 
	print STCH '		if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to run CDHIT.\"  \n";
	print STCH "			exit 1\n"; 
	# if finished succussfully, make a file named QC_finished
	print STCH "		else\n";
	print STCH "			touch $status_log/j5_CDHIT_finished\n";
	print STCH "		fi\n";

	#if CD_HIT report file exists, check the completeness of output
	print STCH "	else\n";
	print STCH '		tail -5 ${CDHIT_REPORT} | grep "program completed"', "\n";
	print STCH '		OUT=$?',"\n";
	print STCH '		if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "		then\n";
	print STCH "			$cdhit -i \${CDHIT_IN} -o \${CDHIT_OUT} ".$cdhit_parameter." > \${CDHIT_REPORT} ", "\n";
	print STCH "			OUT=\$?\n"; 
	print STCH '			if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to run CDHIT.\"  \n";
	print STCH "				exit 1\n"; 
	# if finished succussfully, make a file named QC_finished
	print STCH "			else\n";
	print STCH "				touch $status_log/j5_CDHIT_finished\n";
	print STCH "			fi\n";
	print STCH "		else\n"; # cd-hit report file is finished
	print STCH "			touch $status_log/j5_CDHIT_finished\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi\n";
	close STCH;

	# here to submit
	$current_job_file1 = "j5_".$sample_name."_CDHIT.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";         
#	print STCH "#SBATCH --mem-per-cpu=20G\n";
	print STCH "#SBATCH --mem=20G\n"; # request total memory of 20G
	print STCH "#SBATCH --cpus-per-task=8\n"; # request 8 CPU
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".error\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load checkpoint\n";
	print STCH "module load cd-hit\n";
	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}
	close STCH;
    $last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#######################################################################################
# This step will run Tantan to mask low complexity sequence
sub run_Tantan{
	$current_job_file1 = "j6_".$sample_name."_Tanta.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";         
	print STCH "#SBATCH --mem-per-cpu=20G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load tantan\n";
	print STCH "#!/bin/bash\n\n";
	print STCH "TANTAN_IN=$sample_dir"."/".$sample_name.".QCed.cdhit.fa \n";
	print STCH "TANTAN_OUT=$sample_dir"."/".$sample_name.".QCed.cdhit.tantan.fa \n";

	print STCH "if [ ! -e $status_log/j6_TANTAN_finished ]\n";
	print STCH "then\n";
	print STCH "	$Tantan -x N  \${TANTAN_IN} >  \${TANTAN_OUT} ", "\n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error trying to run Tantan.\"  \n";
	print STCH "		exit 1\n"; 
	# if finished succussfully, make log file
	print STCH "	else\n";
	print STCH "		touch $status_log/j6_TANTAN_finished\n";
	print STCH "	fi\n";
	print STCH "fi\n\n";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
sub seq_QC_TANTAN {
	$current_job_file1 = "j7_".$sample_name."_QC_TANTAN.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";         
	print STCH "#SBATCH --mem-per-cpu=16G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "SAMPLE_DIR=".$sample_dir."\n";

	print STCH "if [ ! -e $status_log/j7_QC_TANTAN_finished ]\n";
	print STCH "then\n";
	print STCH "	".$run_script_path."Tantan_Sequence_Quality_Control.pl ".$sample_dir."\n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error trying to run Tantan.\"  \n";
	print STCH "		exit 1\n"; 
	# if finished succussfully, make log file
	print STCH "	else\n";
	print STCH "		touch $status_log/j7_QC_TANTAN_finished \n";
	print STCH "	fi\n";
	print STCH "fi\n\n";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}


#####################################################################################
# This step will split unique reads to small files for RepeatMasker.
sub prepare_for_RM {
	$current_job_file1 = "j8_".$sample_name."_PRM.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";         
	print STCH "#SBATCH --mem-per-cpu=10G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "IN=".${sample_name}.".QCed.cdhit.tantan.goodSeq.fa\n";
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "RM_DIR=".$repeatmasker_dir."\n\n";

	print STCH "if [ ! -e $status_log/j8_finished_prepare_repeatmasker_split ]\n";
	print STCH "then\n";

	# if the directory already exist, remove it.
	print STCH "	if [ -d \${RM_DIR} ]\n";
	print STCH "	then\n";
	print STCH "		rm -rf \${RM_DIR} \n";
	print STCH "	fi\n";

	# make repeatmasker directory
	print STCH "	mkdir \${RM_DIR} \n\n";

	# split into smaller files
	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d $sample_dir -i  \${IN} -o $repeatmasker_dir -f $file_number_of_RepeatMasker -n $num_seq_per_file_Repeatmasker \n";
	print STCH "	CHECK1=\$?\n"; # if finishes value is 0
	print STCH "	".$run_script_path."check_split_RepeatMasker.pl \${SAMPLE_DIR}\n";
	print STCH "	CHECK2=\$?\n"; # if finishes correctly value is 0
	print STCH "	if [ \${CHECK1} -ne 0 -o \${CHECK2} -ne 0 ] \n"; 
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error at preparing for RepeatMasker.\"  \n";
	print STCH "		exit 1\n"; 
	# if finished succussfully, make log file
	print STCH "	else\n";
	print STCH "		touch $status_log/j8_finished_prepare_repeatmasker_split \n";
	print STCH "	fi\n";
	print STCH "fi\n\n";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

########################################################################################
# submit RepeatMasker job array
sub submit_job_array_RM {
	# this is the job script
	my $job_script = $job_files_dir."/j9_".$sample_name."_RM_script.sh";
	open(JSTCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print JSTCH "RMOUT=$repeatmasker_dir/${sample_name}.QCed.cdhit.tantan.goodSeq.fa_file\${SLURM_ARRAY_TASK_ID}.fasta.masked\n";#
	print JSTCH "RMIN=$repeatmasker_dir/${sample_name}.QCed.cdhit.tantan.goodSeq.fa_file\${SLURM_ARRAY_TASK_ID}.fasta\n";#
	print JSTCH "RMOTHER=$repeatmasker_dir/${sample_name}.QCed.cdhit.tantan.goodSeq.fa_file\${SLURM_ARRAY_TASK_ID}.fasta.out\n";
	print JSTCH "RM_dir=".$repeatmasker_dir."\n\n";
	
	print JSTCH 'if [ -e $RMIN ]',"\n"; # input file exist
	print JSTCH "then\n";
	print JSTCH '	if [ ! -e $RMOTHER ]',"\n"; # don't have RepeatMasker output ".out" file, means RepeatMasker never ran or finished
	print JSTCH "	then\n";
	print JSTCH "		$repeat_masker -pa 4 \$RMIN \n"; # run RepeatMasker, -pa number of processors to use
	print JSTCH '     	CHECK=$?',"\n"; # if finished correctly the value will be 0
	print JSTCH "		if [ \${CHECK} -ne 0 ] \n"; 
	print JSTCH "		then\n";
	print JSTCH "			echo \"Fatal Error at RepeatMasker.\"  \n";
	print JSTCH "			exit 1\n"; 
	# if finished succussfully, make log file
	print JSTCH "		fi\n";
	print JSTCH "	fi\n\n";

 	print JSTCH '	if [ ! -e $RMOUT ]',"\n"; #sometimes repeatmasker does not find any repeat in input files, in these cases no .masked file will be generated.
	print JSTCH "	then\n";
	print JSTCH '		cp ${RMIN} ${RMOUT}',"\n";
	print JSTCH "	fi\n";
	print JSTCH "fi\n";
	close JSTCH;

	# here to submit
	$current_job_file1 = "j9_".$sample_name."_RM.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem=10G\n"; # request total memory of 10G
	print STCH "#SBATCH --cpus-per-task=4\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "#SBATCH --workdir=$SLURM_files_dir\n";
	print STCH "#SBATCH --array=1-$file_number_of_RepeatMasker\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load checkpoint\n";
	print STCH "module load ncbi-blast \n";
	print STCH "module load repeatmasker\n";
	print STCH "set -x\n";
	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
sub seq_RMQC {
	$current_job_file1 = "j10_".$sample_name."_RMQC.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=8G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "QC_Record=".$sample_dir."/".$sample_name.".RepeatMasker_check.txt\n\n";

	print STCH "if [ ! -e $status_log/j10_RMQC_output_finished ]\n";
	print STCH "then\n";
	print STCH "	".$run_script_path."RepeatMasker_SequenceQualityControl.pl ".$sample_dir."\n";
	print STCH '	CHECK1=$?',"\n"; # if finishes value is 0
	print STCH "	".$run_script_path."RepeatMasker_check_SequenceQualityControl.pl \${SAMPLE_DIR} > \${QC_Record} \n";
	print STCH '	CHECK2=$?',"\n"; # if finishes correctly value is 0
	print STCH '	if [ ${CHECK1} -ne 0 -o ${CHECK2} -ne 0 ]',"\n"; # either one is not finished correctly
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error at RepeatMasker SequenceQualityControl.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	else\n";
	print STCH "		touch $status_log/j10_RMQC_output_finished \n";
	print STCH "	fi\n";
	print STCH "fi\n";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

################################################################################
sub map_Host_Reference_Genome_BWA{
	# this is the job script
	my $job_script = $job_files_dir."/j11_".$sample_name."_BWA_map_host_genome_script.sh";
	open(JSTCH, ">$job_script") or die $!;
	chmod 0755, $job_script;
	print JSTCH "#!/bin/bash\n";         
	print JSTCH "IN=$sample_dir"."/".$sample_name.".RepeatMasker.goodSeq.unmasked.fa\n";
	print JSTCH "OUTFILE=$sample_dir"."/".$sample_name.".RepeatMasker.goodSeq.RefGenome.mapped.sam\n";
	print JSTCH "TIMEFILE=$sample_dir"."/".$sample_name.".RepeatMasker.goodSeq.RefGenome.map_time.txt\n";

	# Check if mapping finished successfully. If finished successfully, 
	# resubmitting this job will skip mapping step and will not overwrite 
	# the results.
	print JSTCH "if [ ! -e $status_log/j11_BWA_map_host_genome_finished ]\n";
	print JSTCH "then\n";
	print JSTCH "	date > \${TIMEFILE}\n";
	# use -L 100,100 to reduce spuriously mapped sequences
	print JSTCH "	$bwa mem -L 100,100 -k 15 -t \$SLURM_CPUS_PER_TASK  $reference_genome  \${IN} >  \${OUTFILE}  \n";
	print JSTCH "	OUT=\$?\n"; 
	print JSTCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print JSTCH "	then\n";
	print JSTCH "		echo \"Fatal Error trying to run bwa.\"  \n";
	print JSTCH "		exit 1\n"; 
	print JSTCH "	else\n";
	print JSTCH "		touch $status_log/j11_BWA_map_host_genome_finished\n";
	print JSTCH "	fi\n";
	#if mapping finished succussfully, make a file named finished
	print JSTCH "	date >> \${TIMEFILE}\n";
	print JSTCH "fi\n\n";
	close JSTCH;

	# here to submit
	$current_job_file1 = "j11_".$sample_name."_BWA_map_host_genome.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
#	print STCH "#SBATCH --mem-per-cpu=10G\n";
	print STCH "#SBATCH --mem=20G\n"; # request total memory of 20G
	print STCH "#SBATCH --cpus-per-task=12\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1."\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load checkpoint\n";
	print STCH "module load bwa\n";
	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}
	close STCH;
    $last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#################################################
# extract unmapped reads
sub Extract_Host_unmapped_reads{
	$current_job_file1 = "j12_".$sample_name."_Extract_Host_genome_unmapped_reads.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=8G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
#	print STCH "module load samtools hydra\n";
	print STCH "module load samtools bedtools/2.22.1\n";
	print STCH "IN=$sample_dir"."/".$sample_name.".RepeatMasker.goodSeq.RefGenome.mapped.sam \n";
	print STCH "OUTFILE1=$sample_dir"."/".$sample_name."_RepeatMasker_goodSeq_RefGenome_unmapped.bam\n";
	print STCH "OUTFILE2=$sample_dir"."/".$sample_name.".RepeatMasker.goodSeq.RefGenome.unmapped.fastq\n";
	print STCH "OUTFILE3=$sample_dir"."/".$sample_name.".RepeatMasker.goodSeq.RefGenome.unmapped.unmasked.fasta\n";
	print STCH "OUTFILE4=$sample_dir"."/".$sample_name.".RepeatMasker.goodSeq.RefGenome.unmapped.masked.fasta\n";

	print STCH "if [ ! -e $status_log/j12_Extract_Host_unmapped_reads_finished ]\n";
	print STCH "then\n";
 	# extract unmapped
	print STCH "	$samtools view -b -S -f 4  \${IN} > \${OUTFILE1} \n";
	print STCH "	OUT1=\$?\n";
	# convert bam to fastq
	print STCH "	$bamToFastq  -i  \${OUTFILE1}  -fq \${OUTFILE2}  \n";
	print STCH "	OUT2=\$?\n";
	# convert fastq to fasta
	print STCH "	".$run_script_path."FASTQ_to_FASTA_after_BWA.pl \${OUTFILE2} > \${OUTFILE3}\n"; 
	print STCH "	OUT3=\$?\n";

	# change unmasked sequence to masked sequence for BLAST
	print STCH "	".$run_script_path."Unmasked_seq_to_masked_seq.pl $sample_dir \n"; 
	print STCH "	OUT4=\$?\n";
	# did not finish successfully
	print STCH '	if [ ${OUT1} -ne 0 ] || [ ${OUT2} -ne 0 ] || [ ${OUT3} -ne 0 ] || [ ${OUT4} -ne 0 ] ',"\n"; 
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error trying to extract host unmapped reads.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	else\n";
	print STCH "		touch $status_log/j12_Extract_Host_unmapped_reads_finished\n";
	print STCH "	fi\n";
	print STCH "fi\n\n";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
 }

#######################################################################################
# Split files to small ones for MegaBlast against reference genome to identify more divergent
# seqeunces that share significant homology with reference genome.  
sub split_for_MegaBLAST_RefGenome{
	$current_job_file1 = "j13_".$sample_name."_split_MegaBLAST.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=8G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "IN=".${sample_name}.".RepeatMasker.goodSeq.RefGenome.unmapped.masked.fasta\n\n";
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "MegaBLAST_DIR=".$MegaBLAST_dir_host."\n";

	print STCH "if [ ! -e $status_log/j13_split_for_MegaBLAST_finished ]\n";
	print STCH "then\n";
	# if the directory already exist, remove it.
	print STCH "	if [ -d \${MegaBLAST_DIR} ]\n";
	print STCH "	then\n";
	print STCH "		rm -rf \${MegaBLAST_DIR} \n";
	print STCH "	fi\n";

	# make the directory
	print STCH "	mkdir \${MegaBLAST_DIR} \n\n";

	#split fasta 
	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d $sample_dir -i \${IN}  -o \${MegaBLAST_DIR}  -f $file_number_of_MegaBLAST_host -n $num_seq_per_file_MegaBLAST_host \n";
	print STCH '	CHECK1=$?',"\n"; # If correctly completed, value should be 0.
	print STCH '	if [ ${CHECK1} -ne 0 ]',"\n"; # Did not complete correctly
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error trying to extract split for MegaBLAST host genome.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	else\n";
	print STCH "		touch $status_log/j13_split_for_MegaBLAST_finished\n";
	print STCH "	fi\n";
	print STCH "fi\n";

	# remove unwanted intermediate files
#	print STCH "OUTFILE2=$sample_dir"."/".$sample_name.".RepeatMasker.goodSeq.RefGenome.unmapped.fastq\n";
#	print STCH "OUTFILE3=$sample_dir"."/".$sample_name.".RepeatMasker.goodSeq.RefGenome.unmapped.unmasked.fasta\n";
#	print STCH " rm \${OUTFILE2}\n";
#	print STCH " rm \${OUTFILE3}\n";

	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#######################################################################################
# Run MegaBlast against reference genome to remove more host sequences.
sub MegaBLAST_RefGenome{
	# this is the job script
	my $job_script = $job_files_dir."/j14_".$sample_name."_MegaBLAST_HOST_script.sh";
	open(STCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print STCH "BlastNOUT=$MegaBLAST_dir_host/${sample_name}.RepeatMasker.goodSeq.RefGenome.unmapped.masked.fasta_file".'${SLURM_ARRAY_TASK_ID}',".MegaBLAST.out\n"; #full path
	print STCH "QUERY=$MegaBLAST_dir_host/${sample_name}.RepeatMasker.goodSeq.RefGenome.unmapped.masked.fasta_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n\n";

	print STCH 'if [ -e $QUERY ]',"\n"; # make sure query file exist
	print STCH "then\n";
	#if the output file does not exist, run and check the completeness of the out file
	print STCH '	if [ ! -e $BlastNOUT ]',"\n"; # no output file
	print STCH "	then\n";
#	print STCH "		$blastn -task megablast -outfmt 7  -db $reference_genome -evalue 1e-3  -word_size 16 -num_descriptions 1 -num_alignments 1  -query \${QUERY} -out \${BlastNOUT} ","\n";

	print STCH "		$blastn -task megablast  -num_threads \$SLURM_CPUS_PER_TASK  -db $reference_genome -perc_identity 85 -evalue 1e-9  -word_size 16 -num_descriptions 1 -num_alignments 1  -query \${QUERY} -out \${BlastNOUT} ","\n";
	print STCH '		CHECK1=$?',"\n"; # If correctly completed, value should be 0.
	print STCH '		if [ ${CHECK1} -ne 0 ]',"\n"; # Did not complete correctly
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to extract MegaBLAST host genome.\"  \n";
	print STCH "			exit 1\n"; 
	print STCH "		fi\n";

	#if the output file exist, check for completeness of output file
	print STCH "	else\n"; 
	print STCH '		tail -5 ${BlastNOUT}|grep Matrix',"\n";
	print STCH '		CHECK1=$?',"\n";
	print STCH '		if [ ${CHECK1} -eq 1 ] ',"\n";
	print STCH "		then\n";
#	print STCH "		$blastn -task megablast -outfmt 7  -db $reference_genome -evalue 1e-3  -word_size 16 -num_descriptions 1 -num_alignments 1  -query \${QUERY} -out \${BlastNOUT} ","\n";
	print STCH "			$blastn -task megablast -num_threads \$SLURM_CPUS_PER_TASK    -db $reference_genome -perc_identity 85  -evalue 1e-9  -word_size 16 -num_descriptions 1 -num_alignments 1  -query \${QUERY} -out \${BlastNOUT} ","\n";
	print STCH '			CHECK1=$?',"\n";
	print STCH '			if [ ${CHECK1} -ne 0 ]',"\n"; # Did not complete correctly
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to MegaBlast MegaBLAST host genome.\"  \n";
	print STCH "				exit 1\n"; 
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;

 	# here to submit
	$current_job_file1 = "j14_".$sample_name."_MegaBLAST_HOST.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem=20G\n";
	print STCH "#SBATCH --cpus-per-task=8\n";
#	print STCH "#SBATCH --mem-per-cpu=16G\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "#SBATCH --array=1-$file_number_of_MegaBLAST_host\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "module load ncbi-blast\n";
	print STCH "module load checkpoint\n";
	print STCH "set -x\n";

	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

###################################################################################
# parsing MegaBLAST reference genome output file
sub parse_MegaBLAST_RefGenome {
	$current_job_file1 = "j15_".$sample_name."_parseMegaBLAST.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=4G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "#SBATCH --array=1-$file_number_of_MegaBLAST_host\n";
	print STCH "module load bio-perl\n";
	print STCH "Ref_DIR=$MegaBLAST_dir_host\n";
	print STCH "MegaBlastOUT=${sample_name}.RepeatMasker.goodSeq.RefGenome.unmapped.masked.fasta_file".'${SLURM_ARRAY_TASK_ID}',".MegaBLAST.out\n"; # name only, not full path
	print STCH "MegaBlastIN=$MegaBLAST_dir_host/${sample_name}.RepeatMasker.goodSeq.RefGenome.unmapped.masked.fasta_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n";

	print STCH "PARSED=$MegaBLAST_dir_host/${sample_name}.RepeatMasker.goodSeq.RefGenome.unmapped.masked.fasta_file".'${SLURM_ARRAY_TASK_ID}',".MegaBLAST.parsed\n";

	#if the parsed file does not exist, run parser and check the completeness of the parsed file
	print STCH 'if [ -e $MegaBlastIN ]',"\n";  
	print STCH "then\n";

	#if the parsed file does not exist, run parser and check the completeness of the parsed file
	print STCH '	if [ ! -e $PARSED ]',"\n";
	print STCH "	then\n";
	print STCH "		".$run_script_path."MegaBLAST_HG_parser.pl \${Ref_DIR} \${MegaBlastOUT} $refrence_genome_taxonomy \n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${MegaBlastIN} \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed, if not completed, exit
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n"; # value is 0 if completed correctly.
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to process \${MegaBlastOUT}.\"  \n";
	print STCH '			exit 1',"\n";
	print STCH "		fi\n";

	#if the parsed file exists, check the completeness of the parsed file
	print STCH "	else\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${Ref_DIR} \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not correctly completed exit.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n"; # value is 0 if completed correctly.
	print STCH "		then\n"; # rerun the parser
	print STCH "		".$run_script_path."MegaBLAST_HG_parser.pl \${Ref_DIR} \${MegaBlastOUT} $refrence_genome_taxonomy \n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${MegaBlastIN} \${PARSED}\n";
	print STCH '			CHECK=$?',"\n";
	#check if parsed file is completed, if not completed, exit
	print STCH '			if [ ${CHECK} -ne 0 ]',"\n"; # value is 0 if completed correctly.
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to process \${MegaBlastOUT}.\"  \n";
	print STCH '				exit 1',"\n";
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#######################################################################################
sub extract_RefGenomeFiltered_reads{
	$current_job_file1 = "j16_".$sample_name."_ExtractRefGenomeFiltered.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=15G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}

	print STCH "if [ ! -e $status_log/j16_extract_RefGenomeFiltered_finished ]\n";
	print STCH "then\n";

	# file contains reads left for further analysis
	print STCH "	SAMPLE_DIR=".$sample_dir."\n";
	print STCH "	RefGenomeFiltered=".$sample_name.".RefGenomeFiltered.fa\n\n";
	print STCH "	QC_Record=".$sample_dir."/".$sample_name.".MegaBLAST_HOST_finish_check.txt\n\n";

	# check to make sure all MegaBLAST finished and parser finished correctly
	print STCH "	".$run_script_path."check_MegaBLAST_HOST_Finished_correctly.pl \${SAMPLE_DIR} > \${QC_Record}  \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : MegaBLAST_HOST did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";

	# pool the MegaBLAST reference genome filtered reads into a single file 	
	print STCH "	if [ -e $sample_dir/\${RefGenomeFiltered} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${RefGenomeFiltered}\n";
	print STCH "	fi\n";
	print STCH "	cat $MegaBLAST_dir_host/*.RefGenomeFiltered.fa >> $sample_dir/\${RefGenomeFiltered}\n\n";
	print STCH "	touch $status_log/j16_extract_RefGenomeFiltered_finished\n";
	print STCH "fi\n\n";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
# blastn against viral-only database to extract potential viral-sequences
sub split_for_BLASTN_VIRUSDB{
	$current_job_file1 = "j17_".$sample_name."_BLASTN_VIRUSDB_split.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=15G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "BLASTN_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$BLASTN_VIRUSDB_DIR_SUFFIX."\n";
	print STCH "IN=".$sample_name.".RefGenomeFiltered.fa\n\n";

	print STCH "if [ ! -e $status_log/j17_BLASTN_VIRUSDB_split_finished ]\n"; # job never finished
	print STCH "then\n";
	print STCH "	if [ -d \${BLASTN_VIRUSDB_DIR} ] \n"; # If the dir already exists, remove.
	print STCH "	then\n";
	print STCH "		rm -rf \${BLASTN_VIRUSDB_DIR}\n";
	print STCH "	fi\n\n";

	# make the directory
	print STCH "	mkdir \${BLASTN_VIRUSDB_DIR}\n";
	
	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d $sample_dir -i \${IN}  -o \${BLASTN_VIRUSDB_DIR} -f $file_number_of_BLASTN_VIRUSDB -n $num_seq_per_file_BLASTN_VIRUSDB \n";
	print STCH '	CHECK=$?',"\n";
	print STCH '	if [ ${CHECK} -ne 0 ]',"\n"; 
	print STCH "	then\n";
	# split to -n number of files, this number should be consistent with the number of 
	# blastn job array submitted bellow
	print STCH "		echo \"Fatal Error trying to split for BLASTN_VIRUSDB\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	else\n";
	print STCH "		touch $status_log/j17_BLASTN_VIRUSDB_split_finished \n";
	print STCH "	fi\n";
	print STCH "fi\n";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
sub submit_job_array_BLASTN_VIRUSDB{
	# this is the job script
	my $job_script = $job_files_dir."/j18_".$sample_name."_BLASTN_VIRUSDB_run_script.sh";
	open(STCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print STCH "BLASTN_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$BLASTN_VIRUSDB_DIR_SUFFIX."\n";
	print STCH "BlastNOUT=",'${BLASTN_VIRUSDB_DIR}',"/",$sample_name.".RefGenomeFiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastnv.out\n"; # full path
	print STCH "QUERY=",'${BLASTN_VIRUSDB_DIR}',"/".$sample_name.".RefGenomeFiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n\n";
   
	print STCH 'if [ -s $QUERY ]',"\n"; #check if the file is empty
	print STCH "then\n";
	#if the output file does not exist, run and check the completeness of the output file
	print STCH '	if [ ! -e $BlastNOUT ]',"\n";
	print STCH "	then\n";
	print STCH "		$blastn -num_threads \$SLURM_CPUS_PER_TASK  -evalue  1e-4 -show_gis -task blastn -query \${QUERY} -out \${BlastNOUT} -db $NT_VIRUS -db_soft_mask 100 ","\n";
	print STCH '		CHECK=$?',"\n";
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n"; 
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error : BLASTn did not finish correctly.\"  \n";
	print STCH "			exit 1\n"; 
	print STCH "		fi\n";

	#if the output file exists, check the completeness of the output file
	print STCH "	else\n";
	print STCH '		tail -5 ${BlastNOUT}|grep Matrix',"\n";
	print STCH '		CHECK1=$?',"\n";
	print STCH '		if [ ${CHECK1} -eq 1 ] ',"\n";
	print STCH "		then\n";
	print STCH "			$blastn  -num_threads \$SLURM_CPUS_PER_TASK  -evalue  1e-4 -show_gis  -task blastn -query \${QUERY} -out \${BlastNOUT} -db $NT_VIRUS","\n";
	print STCH '			CHECK=$?',"\n";
	print STCH '			if [ ${CHECK} -ne 0 ]',"\n"; 
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error : BLASTn did not finish correctly.\"  \n";
	print STCH "				exit 1\n"; 
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;

 	# here to submit
	$current_job_file1 = "j18_".$sample_name."_BLASTN_VIRUSDB_run.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem=16G\n";
	print STCH "#SBATCH --cpus-per-task=4\n";
	print STCH "#SBATCH --array=1-$file_number_of_BLASTN_VIRUSDB\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "module load ncbi-blast\n";
	print STCH "module load checkpoint\n";
#	print STCH "set -x\n";
	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}

	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
sub parse_BLASTN_VIRUSDB {
	$current_job_file1 = "j19_".$sample_name."_Parse_BLASTN_VIRUSDB.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=1G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "#SBATCH --array=1-$file_number_of_BLASTN_VIRUSDB\n";
	print STCH "module load bio-perl\n";
	print STCH "BLASTN_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$BLASTN_VIRUSDB_DIR_SUFFIX."\n";
	print STCH "BlastNOUT=",$sample_name.".RefGenomeFiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastnv.out\n";#name only, not full path
	print STCH "BlastNIN=",'${BLASTN_VIRUSDB_DIR}',"/",$sample_name.".RefGenomeFiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n";#full path
	print STCH "PARSED=",'${BLASTN_VIRUSDB_DIR}',"/".$sample_name.".RefGenomeFiltered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastnv.parsed\n\n";

	print STCH 'if [ -s $BlastNIN ]',"\n"; 
	print STCH "then\n";

	#if the parsed file does not exist, run parser and check the completeness of the parsed file
	print STCH '	if [ ! -f $PARSED ]',"\n";
	print STCH "	then\n";
	print STCH "		".$run_script_path."BLASTn_VIRUSDB_parser.pl.sqlite ".$sample_dir."/".$sample_name.$BLASTN_VIRUSDB_DIR_SUFFIX." \${BlastNOUT}\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${BlastNIN} \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed, if not completed, exit
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n"; # value is 0 if completed correctly.
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to process \${BlastNOUT}.\"  \n";
	print STCH '			exit 1',"\n";
	print STCH "		fi\n";

	#if the parsed file exists, check the completeness of the parsed file
	print STCH "	else\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${BlastNIN}  \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not correctly completed exit.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n"; # value is 0 if completed correctly.
	print STCH "		then\n"; # rerun the parser
	print STCH "			".$run_script_path."BLASTn_VIRUSDB_parser.pl.sqlite ".$sample_dir."/".$sample_name.$BLASTN_VIRUSDB_DIR_SUFFIX." \${BlastNOUT}\n";
	print STCH "			".$run_script_path."check_Blast_parsed_file.pl \${BlastNIN}  \${PARSED}\n";
	print STCH '			CHECK=$?',"\n";
	#check if parsed file is completed, if not completed, exit
	print STCH '			if [ ${CHECK} -ne 0 ]',"\n"; # value is 0 if completed correctly.
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to process \${BlastNOUT}.\"  \n";
	print STCH '				exit 1',"\n";
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

########################################################################
# This step split data into smaller files to BLAST against protein Virus database
sub split_for_BLASTX_VIRUSDB {
	$current_job_file1 = "j20_".$sample_name."_split_for_BX_virusDB.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=10G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "BLASTX_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$BLASTX_VIRUSDB_DIR_SUFFIX."\n";
	print STCH "BLASTN_VIRUSDB_Filtered_fa=".$sample_name.".BLASTN_VIRUSDB_Filtered.fa\n";
	print STCH "BLASTN_VIRUSDB_Hit_fa=".$sample_name.".BLASTN_VIRUSDB_HIT.fa\n";
	print STCH "BLASTN_VIRUSDB_Phage_fa=".$sample_name.".BLASTN_VIRUSDB_Phage.fa\n";
	print STCH "BLASTN_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$BLASTN_VIRUSDB_DIR_SUFFIX."\n\n";
	print STCH "QC_Record=".$sample_dir."/".$sample_name.".BLASTN_VIRUSDB_finish_check.txt\n\n";

	print STCH "if [ ! -e $status_log/j20_finished_split_for_BLASTX_VIRUSDB ]\n";
	print STCH "then\n";
	# check to make sure the number of reads add togetehr equals the total input reads
	print STCH "	".$run_script_path."check_BLASTN_VIRUSDB_Finished_correctly.pl \${SAMPLE_DIR} > \${QC_Record}  \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : BLASTN_VIRUSDB did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";

	# if the BLASTX Virus database directory already exist, remove it
	print STCH '	if [ -d $BLASTX_VIRUSDB_DIR ] ',"\n";
	print STCH "	then\n";
	print STCH "		rm -rf \${BLASTX_VIRUSDB_DIR}\n";
	print STCH "	fi\n";

	# make BLASTX Virus database directory
	print STCH "	mkdir \${BLASTX_VIRUSDB_DIR}\n";

	# pool all the BLASTN VIRUSDB hits into a single file
	print STCH "	if [ -e $sample_dir/\${BLASTN_VIRUSDB_Hit_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${BLASTN_VIRUSDB_Hit_fa}\n";
	print STCH "	fi\n";
	print STCH "	cat \${BLASTN_VIRUSDB_DIR}/*.BNVIRUSDB_hit.fa >> $sample_dir/\${BLASTN_VIRUSDB_Hit_fa}\n\n";

	# pool all the BLASTN VIRUSDB filtered reads into a single file
	print STCH "	if [ -e $sample_dir/\${BLASTN_VIRUSDB_Filtered_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${BLASTN_VIRUSDB_Filtered_fa}\n";
	print STCH "	fi\n";
	print STCH "	cat \${BLASTN_VIRUSDB_DIR}/*.BNVIRUSDB_filtered.fa >> $sample_dir/\${BLASTN_VIRUSDB_Filtered_fa}\n\n";

	# pool all the BLASTN VIRUSDB phage hits into a single file
	print STCH "	if [ -e $sample_dir/\${BLASTN_VIRUSDB_Phage_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${BLASTN_VIRUSDB_Phage_fa}\n";
	print STCH "	fi\n";
	print STCH "	cat \${BLASTN_VIRUSDB_DIR}/*.BNVIRUSDB_phage.fa >> $sample_dir/\${BLASTN_VIRUSDB_Phage_fa}\n\n";

	# split into smaller files
	# split to -n number of files, this number should be consistent with 
	# the number of blast job array submitted bellow
	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d \${SAMPLE_DIR} -i  \${BLASTN_VIRUSDB_Filtered_fa} -o \${BLASTX_VIRUSDB_DIR} -f $file_number_of_BLASTX_VIRUSDB  -n $num_seq_per_file_BLASTX_VIRUSDB \n";
	print STCH '	CHECK1=$?',"\n";
	print STCH "	".$run_script_path."check_split_BLASTX_VIRUSDB.pl \${SAMPLE_DIR}\n";
	print STCH '	CHECK2=$?',"\n";
#	print STCH '	echo ${CHECK1}, ${CHECK2}', "\n";
	print STCH '	if [ ${CHECK1} -ne 0 -o ${CHECK2} -ne 0 ]',"\n"; 
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error trying to split_BLASTX_VIRUSDB.\"  \n";
	print STCH '		exit 1',"\n";
	print STCH "	else\n";
	print STCH "		touch $status_log/j20_finished_split_for_BLASTX_VIRUSDB \n";
	print STCH "	fi\n";
	print STCH "fi\n";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
# run blastx againt virus-only databse
sub submit_job_array_BLASTX_VIRUSDB{
	# this is the job script
	my $job_script = $job_files_dir."/j21_".$sample_name."_BX_VIRUSDB_script.sh";
	open(STCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print STCH "VIRUS_DIR=".$sample_dir."/".$sample_name.$BLASTX_VIRUSDB_DIR_SUFFIX."\n";
	print STCH "BLASTX_VIRUSDB_OUT=",'${VIRUS_DIR}',"/",$sample_name.".BLASTN_VIRUSDB_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".BLASTX_VIRUSDB.out\n";#full path
	print STCH "QUERY=",'${VIRUS_DIR}',"/".$sample_name.".BLASTN_VIRUSDB_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n\n";

	print STCH 'if [ -s $QUERY ]',"\n"; # check if the file is empty
	print STCH "then\n";
	#if the output file does not exist, run and check the completeness of the output file
	print STCH '	if [ ! -f $BLASTX_VIRUSDB_OUT ]',"\n";
	print STCH "	then\n";
	print STCH "		$blastx  -num_threads \$SLURM_CPUS_PER_TASK  -evalue 1e-2 -show_gis   -query \${QUERY} -out \${BLASTX_VIRUSDB_OUT} -db $NR_VIRUS -db_soft_mask 100","\n";
	print STCH '		CHECK1=$?',"\n";
	print STCH '		if [ ${CHECK1} -ne 0 ]',"\n"; 
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to blastx VIRUSDB.\"  \n";
	print STCH '			exit 1',"\n";
	print STCH "		fi\n";

	#if the output file exists, check the completeness of the output file
	print STCH "	else\n";
	print STCH '		tail -5 ${BLASTX_VIRUSDB_OUT}|grep Matrix',"\n";
	print STCH '		CHECK1=$?',"\n";
	print STCH '		if [ ${CHECK1} -ne 0 ]',"\n"; 
	print STCH "			then\n";
#	print STCH "			$blastx -num_threads \$SLURM_CPUS_PER_TASK  -evalue 1e-2 -show_gis  -query \${QUERY} -out \${BLASTX_VIRUSDB_OUT} -db $NR_VIRUS ","\n";
	print STCH "			$blastx  -num_threads \$SLURM_CPUS_PER_TASK  -evalue 1e-2 -show_gis   -query \${QUERY} -out \${BLASTX_VIRUSDB_OUT} -db $NR_VIRUS -db_soft_mask 100","\n";
	print STCH '			CHECK2=$?',"\n";
	print STCH '			if [ ${CHECK2} -ne 0 ]',"\n"; 
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to blastx VIRUSDB.\"  \n";
	print STCH '				exit 1',"\n";
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi\n";
	close STCH;
	
 	# here to submit
	$current_job_file1 = "j21_".$sample_name."_BX_VIRUSDB.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem=16G\n";
	print STCH "#SBATCH --cpus-per-task=4\n"; # request 4 CPU
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	print STCH "#SBATCH --array=1-$file_number_of_BLASTX_VIRUSDB\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load ncbi-blast\n";
	print STCH "module load checkpoint\n";
#	print STCH "set -x\n";
	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}


#####################################################################################
# parse blastx against virus-only database output file
sub parse_BLASTX_VIRUSDB {
	$current_job_file1 = "j22_".$sample_name."_Parse_BLASTX_VIRUSDB.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=2G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "#SBATCH --array=1-$file_number_of_BLASTX_VIRUSDB\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load bio-perl\n";
#	print STCH "set -x\n";
	print STCH "VIRUS_DIR=".$sample_dir."/".$sample_name.$BLASTX_VIRUSDB_DIR_SUFFIX."\n";
	print STCH "BLASTX_VIRUSDB_OUT=",$sample_name.".BLASTN_VIRUSDB_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".BLASTX_VIRUSDB.out\n"; #name only, not full path
	print STCH "BLASTX_VIRUSDB_IN=",'${VIRUS_DIR}',"/",$sample_name.".BLASTN_VIRUSDB_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n"; #full path
	print STCH "PARSED=",'${VIRUS_DIR}',"/".$sample_name.".BLASTN_VIRUSDB_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".BLASTX_VIRUSDB.parsed\n\n";

	print STCH 'if [ -s $BLASTX_VIRUSDB_IN ]',"\n"; # input file exist
	print STCH "then\n";
	#if the parsed file does not exist, run parser and check the completeness of the parsed file
	print STCH '	if [ ! -f $PARSED ]',"\n";
	print STCH "	then\n";
	print STCH "		".$run_script_path."BLASTx_VIRUSDB_parser.pl.sqlite \${VIRUS_DIR} \${BLASTX_VIRUSDB_OUT}\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${BLASTX_VIRUSDB_IN}  \${PARSED}\n";
	print STCH '		CHECK=$?',"\n"; # if finished correctly the value should be 0.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n"; # parsing did not complete correctly  
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to process \${BLASTX_VIRUSDB_OUT}.\"  \n";
	print STCH "			exit 1\n"; # exit
	print STCH "		fi\n";

	#if the parsed file exists, check the completeness of the parsed file
	print STCH "	else\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${BLASTX_VIRUSDB_IN}  \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n"; # parsing did not complete correctly  
	print STCH "		then\n";
	print STCH "			".$run_script_path."BLASTx_VIRUSDB_parser.pl.sqlite \${VIRUS_DIR} \${BLASTX_VIRUSDB_OUT}\n";
	print STCH "			".$run_script_path."check_Blast_parsed_file.pl \${BLASTX_VIRUSDB_IN}  \${PARSED}\n";
	print STCH '			CHECK=$?',"\n"; # if finished correctly the value should be 0.
	print STCH '			if [ ${CHECK} -ne 0 ]',"\n"; # parsing did not complete correctly  
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to process \${BLASTX_VIRUSDB_OUT}.\"  \n";
	print STCH "				exit 1\n"; # exit
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
sub pool_seq_for_mapping_Bacteria {
 	$current_job_file1 = "j23_".$sample_name."_pool_seq_for_mapping_Bacteria.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=6G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "VIRUSDB_HIT_all=".$sample_name."_VIRUSDB_HIT.fa\n";
	print STCH "BLASTN_VIRUSDB_Hit_fa=".$sample_name.".BLASTN_VIRUSDB_HIT.fa\n";
	print STCH "BLASTX_VIRUSDB_Hit_fa=".$sample_name.".BLASTX_VIRUSDB_HIT.fa\n";
	print STCH "BLASTX_VIRUSDB_Phage_fa=".$sample_name.".BLASTX_VIRUSDB_Phage.fa\n";
	print STCH "BLASTX_VIRUSDB_Filtered_fa=".$sample_name.".BLASTX_VIRUSDB_Filtered.fa\n";
	print STCH "BLASTX_VIRUSDB_DIR=".$sample_dir."/".$sample_name.$BLASTX_VIRUSDB_DIR_SUFFIX."\n\n";
	print STCH "QC_Record=".$sample_dir."/".$sample_name.".BLASTX_VIRUSDB_finish_check.txt\n\n";

	print STCH "if [ ! -e $status_log/j23_finished_pool_seq_for_mapping_Bacteria ]\n"; # the job never finished
	print STCH "then\n";

	# check to make sure the number of reads add togetehr equals the total input reads
	print STCH "	".$run_script_path."check_BLASTX_VIRUSDB_Finished_correctly.pl \${SAMPLE_DIR} > \${QC_Record}  \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : BLASTX_VIRUSDB did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";

	# pool all the BLASTX VIRUSDB hits into a single file
	print STCH "	if [ -e $sample_dir/\${BLASTX_VIRUSDB_Hit_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${BLASTX_VIRUSDB_Hit_fa} \n";
	print STCH "	fi\n";
	print STCH "	cat \${BLASTX_VIRUSDB_DIR}/*.BXVIRUSDB_hit.fa >> $sample_dir/\${BLASTX_VIRUSDB_Hit_fa} \n\n";

	# pool all the BLASTX VIRUSDB filtered into a single file
	print STCH "	if [ -e $sample_dir/\${BLASTX_VIRUSDB_Filtered_fa} ] \n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${BLASTX_VIRUSDB_Filtered_fa} \n";
	print STCH "	fi\n";
	print STCH "	cat \${BLASTX_VIRUSDB_DIR}/*.BXVIRUSDB_filtered.fa >> $sample_dir/\${BLASTX_VIRUSDB_Filtered_fa} \n\n";

	# pool all the BLASTX VIRUSDB phage hits into a single file
	print STCH "	if [ -e $sample_dir/\${BLASTX_VIRUSDB_Phage_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${BLASTX_VIRUSDB_Phage_fa} \n";
	print STCH "	fi\n";
	print STCH "	cat \${BLASTX_VIRUSDB_DIR}/*.BXVIRUSDB_phage.fa >> $sample_dir/\${BLASTX_VIRUSDB_Phage_fa} \n\n";

	# pool BLASTX VIRUSDB hits and BLASTN VIRUSDB hits into a single file
	print STCH "	if [ -e $sample_dir/\${VIRUSDB_HIT_all} ]\n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${VIRUSDB_HIT_all} \n";
	print STCH "	fi\n";
	print STCH "	cat  $sample_dir/\${BLASTX_VIRUSDB_Hit_fa} $sample_dir/\${BLASTN_VIRUSDB_Hit_fa} > \${SAMPLE_DIR}/\${VIRUSDB_HIT_all} \n";

	# change masked sequence to unmasked sequence for mapping to Bacteria
	print STCH "	".$run_script_path."Masked_seq_to_unmasked_seq.pl $sample_dir \n"; 
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : Masked_seq_to_unmasked_seq did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
	print STCH "	touch $status_log/j23_finished_pool_seq_for_mapping_Bacteria \n";
	print STCH "fi\n";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}


################################################################################
# map reads to Bacteria reference genome, output unmapped reads in .soap_unmapped.fa file
# mapping using single-end mode
sub map_Bacteria_Genome{
	$current_job_file1 = "j24_".$sample_name."_Map_Bacteria_Genome.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n"; 
	print STCH "#SBATCH --mem=24G\n";
	print STCH "#SBATCH --cpus-per-task=12\n"; # request 4 CPU
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load bwa\n";
	print STCH "IN=$sample_dir"."/".$sample_name."_VIRUSDB_HIT_unmasked.fa\n";
	print STCH "OUTFILE=$sample_dir"."/".$sample_name."_VIRUSDB_HIT_Bacteria_mapped.sam\n";
#	print STCH "UNMAPPED=$sample_dir"."/".$sample_name.".NonPhage.Bacteria.unmapped.fa\n\n";
	print STCH "TIMEFILE=$sample_dir"."/".$sample_name."_VIRUSDB_HIT_Bacteria_mapping_time.txt\n";

	# Check if SOAPaligner finished successfully. If finished successfully, 
	# resubmitting this job will skip SOAPAligner step and will not overwrite 
	# the results.
	print STCH "if [ ! -e $status_log/j24_Map_Bacteria_finished ]\n";
	print STCH "then\n";
	print STCH "	date > \${TIMEFILE}\n";
	# -t 8 request 8 CPUs
	print STCH "	$bwa mem -k 15 -t \$SLURM_CPUS_PER_TASK  $Bacteria_genome  \${IN} >  \${OUTFILE}  \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : bwa map to bacteria genome did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
	print STCH "	date >> \${TIMEFILE}\n";
	print STCH "	touch $status_log/j24_Map_Bacteria_finished\n";
	print STCH "fi\n\n";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#################################################
# extract unmapped reads
sub Extract_Bacteria_unmapped_reads{
	$current_job_file1 = "j25_".$sample_name."_Extract_Bacteria_unmapped_reads.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=6G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
#	print STCH "module load samtools hydra\n";
	print STCH "module load samtools bedtools/2.22.1\n";
	print STCH "IN=$sample_dir"."/".$sample_name."_VIRUSDB_HIT_Bacteria_mapped.sam \n";
	print STCH "OUTFILE1=$sample_dir"."/".$sample_name."_VIRUSDB_HIT_Bacteria_unmapped.bam\n";
	print STCH "OUTFILE2=$sample_dir"."/".$sample_name."_VIRUSDB_HIT_Bacteria_unmapped.fastq\n";
	print STCH "OUTFILE3=$sample_dir"."/".$sample_name."_VIRUSDB_HIT_Bacteria_unmapped.fasta\n";

	print STCH "if [ ! -e $status_log/j25_Extract_Bacteria_unmapped_reads_finished ]\n";
	print STCH "then\n";
	# extract unmapped
	print STCH "	$samtools view -b -S -f 4  \${IN} > \${OUTFILE1} \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : extract unmapped reads did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";


	# convert bam to fastq
	print STCH "	$bamToFastq  -i  \${OUTFILE1}  -fq \${OUTFILE2}  \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : convert bam to fastq did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
    
	# convert fastq to fasta
	print STCH "	".$run_script_path."FASTQ_to_FASTA_after_BWA.pl \${OUTFILE2} > \${OUTFILE3}\n"; 
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : convert fastq to fasta did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
    
	# change unmasked sequence to masked for BLAST
	print STCH "	".$run_script_path."Unmasked_seq_to_masked_seq_After_BacteriaMapping.pl $sample_dir \n"; 
 	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : change unmasked sequence to masked did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
    
	print STCH "	touch $status_log/j25_Extract_Bacteria_unmapped_reads_finished\n";
	print STCH "fi\n\n";

	# remove unwanted intermediate files
#	print STCH "OUTFILE2=$sample_dir"."/".$sample_name.".VIRUSDB_HIT.Bacteria.unmapped.fastq\n";
#	print STCH " rm \${OUTFILE2}\n";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
 }

#########################################################################################
sub split_for_MegaBLAST_NT{
 	$current_job_file1 = "j26_".$sample_name."_MegaBLAST_split.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=6G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "MegaBLAST_DIR=".$MegaBLAST_dir_NT."\n";
	print STCH "IN=".$sample_name."_VIRUSDB_HIT_Bacteria_unmapped_masked.fasta\n\n";

	print STCH "if [ ! -e $status_log/j26_finished_MegaBLAST_NT_split ]\n"; # the job never finished
	print STCH "then\n";

	# if the directory already exist, remove it.
	print STCH "	if [ -d \${MegaBLAST_DIR} ]\n";
	print STCH "	then\n";
	print STCH "		rm -rf \${MegaBLAST_DIR}\n";
	print STCH "	fi\n";

	# make the directory
	print STCH "	mkdir \${MegaBLAST_DIR}\n\n";

	# split into smaller files
	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d $sample_dir -i \${IN} -o \${MegaBLAST_DIR}   -f $file_number_of_MegaBLAST_NT -n $num_seq_per_file_MegaBLAST_NT \n";
	print STCH "	".$run_script_path."check_split_MegaBLAST_NT.pl $sample_dir \n";
 	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : split_for_MegaBLAST_NT did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
	print STCH "	touch $status_log/j26_finished_MegaBLAST_NT_split \n";
	print STCH "fi\n";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#########################################################################################
# MegaBLAST against full NT to quickly remove false positives
sub submit_job_array_MegaBLAST_NT {
 	# this is the job script
	my $job_script = $job_files_dir."/j27_".$sample_name."_MegaBLAST_NT_script.sh";
	open(STCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print STCH "BLAST_DIR=".$MegaBLAST_dir_NT."\n";
	print STCH "QUERY=",'${BLAST_DIR}',"/".$sample_name."_VIRUSDB_HIT_Bacteria_unmapped_masked.fasta_file".'${SLURM_ARRAY_TASK_ID}'.".fasta\n";
	print STCH "BlastOUT=",'${BLAST_DIR}',"/".$sample_name."_VIRUSDB_HIT_Bacteria_unmapped_masked.fasta_file".'${SLURM_ARRAY_TASK_ID}',".megablast.out\n";

	print STCH 'if [ -s $QUERY ]',"\n"; 
	print STCH "then\n"; 
	#if blast output file does not exist, do blast 
	print STCH '	if [ ! -f $BlastOUT ]',"\n";
	print STCH "	then\n"; # Default -task is `megablast'
	# default MegaBLASt word size is 28
	print STCH "		$blastn -task megablast -num_threads \$SLURM_CPUS_PER_TASK  -evalue  1e-8  -word_size 16  -show_gis -num_descriptions 50 -num_alignments 50  -query  \${QUERY}  -out \${BlastOUT} -db $db_BN","\n";
 	print STCH "		OUT=\$?\n"; 
	print STCH '		if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error : MegaBLAST \${QUERY}  did not finish correctly.\"  \n";
	print STCH "			exit 1\n"; 
	print STCH "		fi\n";
	#if blast output file exists, check the completeness of output
	print STCH "	else\n";
	print STCH '		tail -10 ${BlastOUT}|grep Matrix',"\n";
	print STCH '		CHECK=$?',"\n";
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";
	print STCH "		then\n";
	print STCH "			$blastn -task megablast -num_threads \$SLURM_CPUS_PER_TASK  -evalue  1e-8  -word_size 16  -show_gis -num_descriptions 50 -num_alignments 50  -query  \${QUERY}  -out \${BlastOUT} -db $db_BN","\n";
 	print STCH "			OUT=\$?\n"; 
	print STCH '			if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error : MegaBLAST \${QUERY}  did not finish correctly.\"  \n";
	print STCH "				exit 1\n"; 
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;

 	# here to submit
	$current_job_file1 = "j27_".$sample_name."_MegaBLAST_NT.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n"; 
 	print STCH "#SBATCH --mem=48G\n"; # request total memory of 24G
	print STCH "#SBATCH --cpus-per-task=8\n"; # request 8 CPU
#	print STCH "#SBATCH --mem-per-cpu=20G\n";
	print STCH "#SBATCH --array=1-$file_number_of_MegaBLAST_NT\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load ncbi-blast\n";
	print STCH "module load checkpoint\n";
#	print STCH "set -x\n";

	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#########################################################################################
sub parse_MegaBLAST_NT{
	$current_job_file1 = "j28_".$sample_name."_Parse_MegaBLAST.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=4G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "#SBATCH --array=1-$file_number_of_MegaBLAST_NT\n";
	print STCH "module load bio-perl\n";
	print STCH "BLAST_DIR=".$MegaBLAST_dir_NT."\n";
	print STCH "QUERY=\${BLAST_DIR}/".$sample_name."_VIRUSDB_HIT_Bacteria_unmapped_masked.fasta_file".'${SLURM_ARRAY_TASK_ID}'.".fasta\n";
	print STCH "BlastOUT=".$sample_name."_VIRUSDB_HIT_Bacteria_unmapped_masked.fasta_file".'${SLURM_ARRAY_TASK_ID}',".megablast.out\n"; #name only, not full path
	print STCH "PARSED=\${BLAST_DIR}/".$sample_name."_VIRUSDB_HIT_Bacteria_unmapped_masked.fasta_file".'${SLURM_ARRAY_TASK_ID}'.".megablast.parsed\n\n";

	print STCH 'if [ -s $QUERY ]',"\n"; 
	print STCH "then\n";

	#if the parsed file does not exist, run parser and check the completeness of the parsed file
	print STCH '	if [ ! -f $PARSED ]',"\n"; # no parsed file
	print STCH "	then\n";
	print STCH "		".$run_script_path."MegaBLAST_NT_parser.pl.sqlite  \${BLAST_DIR} \${BlastOUT}\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${QUERY}  \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to process \${BlastNOUT}.\"  \n";
	print STCH "			exit 1\n";
	print STCH "   		fi\n";

	#if the parsed file exists, check the completeness of the parsed file
	print STCH "	else\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${QUERY}  \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, rerun parser.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "		then\n";
	print STCH "			".$run_script_path."MegaBLAST_NT_parser.pl.sqlite  \${BLAST_DIR} \${BlastOUT}\n";
	print STCH "			".$run_script_path."check_Blast_parsed_file.pl \${QUERY}  \${PARSED}\n";
	print STCH '			CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '			if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to process \${BlastNOUT}.\"  \n";
	print STCH "				exit 1\n";
	print STCH "   			fi\n";
	print STCH "   		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
sub pool_split_for_BlastN{
	$current_job_file1 = "j29_".$sample_name."_BN_split.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=6G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "BLASTN_NT_DIR=".$BLASTN_NT_dir."\n";
	print STCH "MegaBLAST_DIR=".$MegaBLAST_dir_NT."\n";
	print STCH "MegaBLAST_Filtered_fa=".$sample_name.".VIRUSDB_HIT_MegaBLAST_NT_Filtered.fa\n\n";
	print STCH "QC_Record=".$sample_dir."/".$sample_name.".MegaBLAST_NT_check.txt\n\n";

	print STCH "if [ ! -e $status_log/j29_finished_split_for_BLASTN_NT ]\n"; # job never finished
	print STCH "then\n";
	# check to make sure all parser finished 
	print STCH "	".$run_script_path."check_MegaBLAST_NT_Finished_correctly.pl \${SAMPLE_DIR} > \${QC_Record}  \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : MegaBLAST_NT did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";

	# if the directory already exist, remove it.
	print STCH "	if [ -d \${BLASTN_NT_DIR} ]\n";
	print STCH "	then\n";
	print STCH "		rm -rf \${BLASTN_NT_DIR}\n";
	print STCH "	fi\n";

	# make the directory
	print STCH "	mkdir \${BLASTN_NT_DIR}\n\n";

	# if the files already exist, remove it
	print STCH "	if [ -f $sample_dir/\${MegaBLAST_Filtered_fa} ] \n";
	print STCH "	then\n";
	print STCH "		rm $sample_dir/\${MegaBLAST_Filtered_fa}\n";
	print STCH "	fi\n";
	print STCH "	cat \${MegaBLAST_DIR}/*.MegaBLAST_filtered.fa >> $sample_dir/\${MegaBLAST_Filtered_fa}\n";
	
	# split into smaller files
	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d $sample_dir -i \${MegaBLAST_Filtered_fa} -o \${BLASTN_NT_DIR}   -f $file_number_of_BLASTN -n $num_seq_per_file_BLASTN \n";
	print STCH "	".$run_script_path."check_split_BN.pl $sample_dir \n";
 	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : split_for_BlastN  did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
	print STCH "	touch $status_log/j29_finished_split_for_BLASTN_NT \n";
	print STCH "fi\n";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################

sub submit_job_array_BLASTN{
	# this is the job script
	my $job_script = $job_files_dir."/j30_".$sample_name."_BN_script.sh";
	open(STCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print STCH "BN_DIR=".$BLASTN_NT_dir."\n";
	print STCH "BlastNOUT=",'${BN_DIR}',"/",$sample_name.".VIRUSDB_HIT_MegaBLAST_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastn.out\n";#full path
	print STCH "QUERY=",'${BN_DIR}',"/".$sample_name.".VIRUSDB_HIT_MegaBLAST_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n\n";

	print STCH 'if [ -s $QUERY ]',"\n"; #check if the file is empty
	print STCH "then\n";
	#if the output file does not exist, run and check the completeness of the output file
	print STCH '	if [ ! -e $BlastNOUT ]',"\n";
	print STCH "	then\n";
	print STCH "		$blastn -evalue  1e-8 -show_gis -num_threads \$SLURM_CPUS_PER_TASK   -task blastn -num_descriptions 50 -num_alignments 50  -query \${QUERY} -out \${BlastNOUT} -db $db_BN","\n";
	print STCH "		OUT=\$?\n"; 
	print STCH '		if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error : BLASTn \${QUERY}  did not finish correctly.\"  \n";
	print STCH "			exit 1\n"; 
	print STCH "		fi\n";

	#if the output file exists, check the completeness of the output file
	print STCH "	else\n";
	print STCH '		tail -5 ${BlastNOUT}|grep Matrix',"\n";
	print STCH '		CHECK=$?',"\n";
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";
	print STCH "		then\n";
	print STCH "			$blastn -evalue  1e-8 -show_gis -num_threads \$SLURM_CPUS_PER_TASK  -task blastn  -num_descriptions 50 -num_alignments 50 -query \${QUERY} -out \${BlastNOUT} -db $db_BN","\n";
 	print STCH "			OUT=\$?\n"; 
	print STCH '			if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error : BLAST \${QUERY}  did not finish correctly.\"  \n";
	print STCH "				exit 1\n"; 
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi\n";
	close STCH;

 	# here to submit
	$current_job_file1 = "j30_".$sample_name."_BN.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";         
	print STCH "#SBATCH --mem=48G\n";
	print STCH "#SBATCH --cpus-per-task=6\n";
	print STCH "#SBATCH --array=1-$file_number_of_BLASTN\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error \n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load ncbi-blast\n";
	print STCH "module load checkpoint\n";
#	print STCH "set -x\n";
	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}


#####################################################################################

sub parse_BLASTN{
	$current_job_file1 = "j31_".$sample_name."_PBN.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=2G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "#SBATCH --array=1-$file_number_of_BLASTN\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load bio-perl\n";
	print STCH "BN_DIR=".$BLASTN_NT_dir."\n";
	print STCH "BlastNOUT=",$sample_name.".VIRUSDB_HIT_MegaBLAST_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastn.out\n";#name only, not full path
	print STCH "BlastNIN=",'${BN_DIR}',"/",$sample_name.".VIRUSDB_HIT_MegaBLAST_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n";#full path
	print STCH "PARSED=",'${BN_DIR}',"/".$sample_name.".VIRUSDB_HIT_MegaBLAST_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastn.parsed\n\n";

	print STCH 'if [ -s $BlastNIN ]',"\n";  
	print STCH "then\n";
	#if the parsed file does not exist, run parser and check the completeness of the parsed file
	print STCH '	if [ ! -e $PARSED ]',"\n";
	print STCH "	then\n";
	print STCH "		".$run_script_path."BLASTn_NT_parser.pl.sqlite ".$BLASTN_NT_dir." \${BlastNOUT}\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${BlastNIN}  \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to process \${BlastNOUT}.\"  \n";
	print STCH "			exit 1\n";
	print STCH "   		fi\n";

	#if the parsed file exists, check the completeness of the parsed file
	print STCH "	else\n";
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${BlastNIN}  \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, rerun parser.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "		then\n";
	print STCH "			".$run_script_path."BLASTn_NT_parser.pl.sqlite ".$BLASTN_NT_dir." \${BlastNOUT}\n";
	print STCH "			".$run_script_path."check_Blast_parsed_file.pl \${BlastNIN}  \${PARSED}\n";
	print STCH '			CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '			if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to process \${BlastNOUT}.\"  \n";
	print STCH "				exit 1\n";
	print STCH "   			fi\n";
	print STCH "   		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
sub pool_split_for_BLASTX{
	$current_job_file1 = "j32_".$sample_name."_BX_split.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=10G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "BX_DIR=".$BLASTX_NR_dir."\n";
	print STCH "BNFiltered_fa=".$sample_name."_BLASTN_NT_Filtered.fa\n";
	print STCH "QC_Record=".$sample_dir."/".$sample_name.".BLAST_NT_check.txt\n\n";

	print STCH "if [ ! -e $status_log/j32_BX_split_finished ]\n"; # job never finished
	print STCH "then\n";
	# check to make sure all parser finished 
	print STCH "	".$run_script_path."check_BLAST_NT_Finished_correctly.pl \${SAMPLE_DIR} > \${QC_Record}  \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : BLAST_NT did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";

	print STCH "	if [ -d \${BX_DIR} ]\n"; # If the dir already exists, remove.
	print STCH "	then\n";
	print STCH "		rm -rf \${BX_DIR}\n";
	print STCH "	fi\n\n";
	print STCH "	mkdir \${BX_DIR}\n";

	print STCH "	if [ -e \${BNFiltered_fa} ]\n";
	print STCH "	then\n";
	print STCH "		rm \$BNFiltered_fa\n";
	print STCH "	fi\n\n";
	print STCH "	cat $BLASTN_NT_dir/*.BNfiltered.fa >> $sample_dir/\${BNFiltered_fa}\n";

	print STCH "	".$run_script_path."FASTA_file_split_given_numberOfFiles.pl -d $sample_dir -i \${BNFiltered_fa}  -o \${BX_DIR} -f $file_number_of_BLASTX -n $num_seq_per_file_BLASTX \n";
	print STCH "	".$run_script_path."check_split_BX.pl ".$sample_dir."\n";
 	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : change unmasked sequence to masked did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
	print STCH "	touch $status_log/j32_BX_split_finished  \n";
	print STCH "fi\n";
	close STCH;
	
    $last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
# this is script for job submission with checkpointing
sub submit_job_array_BLASTX{
	# this is the job script
	my $job_script = $job_files_dir."/j33_".$sample_name."_BX_script.sh";
	open(STCH, ">$job_script") or die $!;
	chmod 0755, $job_script;

	print STCH "set -x\n";
	print STCH "BX_DIR=".$BLASTX_NR_dir."\n";
	print STCH "BlastXOUT=",'${BX_DIR}',"/",$sample_name."_BLASTN_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastx.out\n";#full path
	print STCH "QUERY=",'${BX_DIR}',"/".$sample_name."_BLASTN_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n\n";

	print STCH 'if [ -s $QUERY ]',"\n"; # check if the file is not empty
	print STCH "then\n";

	#if the blastx output file does not exist, run blastx and check the completeness of the blastx output file
	print STCH '	if [ ! -e $BlastXOUT ]',"\n";
	print STCH "	then\n";
	print STCH "		$blastx  -num_threads \$SLURM_CPUS_PER_TASK   -evalue 1e-2 -show_gis -num_descriptions 50 -num_alignments 50   -query \${QUERY} -out \${BlastXOUT} -db $db_BX ","\n";
	print STCH "		OUT=\$?\n"; 
	print STCH '		if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error : BLASTx \${QUERY}  did not finish correctly.\"  \n";
	print STCH "			exit 1\n"; 
	print STCH "		fi\n";
	#if the blastx output file exists, check the completeness of the blastx output file
	print STCH "	else\n";
	print STCH '		tail -50 ${BlastXOUT}|grep Matrix ',"\n";
	print STCH '		CHECK=$?',"\n";
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";
	print STCH "		then\n";
	print STCH "			$blastx -num_threads \$SLURM_CPUS_PER_TASK  -evalue 1e-2 -show_gis -num_descriptions 50 -num_alignments 50 -query \${QUERY} -out \${BlastXOUT} -db $db_BX ","\n";
 	print STCH "			OUT=\$?\n"; 
	print STCH '			if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error : BLASTx \${QUERY}  did not finish correctly.\"  \n";
	print STCH "				exit 1\n"; 
	print STCH "			fi\n";
	print STCH "		fi\n";
	print STCH "	fi\n";
	print STCH "fi\n"; 
	close STCH;
  
	# here to submit
	$current_job_file1 = "j33_".$sample_name."_BX.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem=50G\n";
	print STCH "#SBATCH --cpus-per-task=8\n";
	print STCH "#SBATCH --array=1-$file_number_of_BLASTX\n";
	print STCH "#SBATCH --error=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.error\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load ncbi-blast\n";
	print STCH "module load checkpoint\n";
#	print STCH "module load checkpoint/test\n";
	print STCH "set -x\n";

	if ($use_checkpoint) { 
		print STCH "checkpoint bash  $job_script \n";
	}
	else {
		print STCH " bash  $job_script \n";
	}
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
sub parse_BLASTX{
	$current_job_file1 = "j34_".$sample_name."_PBX.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=4G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".%A_%a.out \n";
	print STCH "#SBATCH --array=1-$file_number_of_BLASTX\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load bio-perl\n";
	print STCH "BX_DIR=".$BLASTX_NR_dir."\n";
	print STCH "BlastXOUT=",$sample_name."_BLASTN_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastx.out\n";#name only, not full path
	print STCH "BlastXIN=",'${BX_DIR}',"/".$sample_name."_BLASTN_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".fasta\n";
	print STCH "PARSED=",'${BX_DIR}',"/".$sample_name."_BLASTN_NT_Filtered.fa_file".'${SLURM_ARRAY_TASK_ID}',".blastx.parsed\n";

	print STCH 'if [ -s $BlastXIN ] ',"\n";  
	print STCH "then\n";
	print STCH '	if [ ! -e $PARSED ] ',"\n"; #blastx.parsed file not exist
	print STCH "	then\n";
	print STCH "		".$run_script_path."BLASTx_NR_parser.pl.sqlite \${BX_DIR} \${BlastXOUT}\n"; #run parser
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${BlastXIN}  \${PARSED}\n";
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "		then\n";
	print STCH "			echo \"Fatal Error trying to process \${MegaBlastOUT}.\"  \n";
	print STCH "			exit 1\n";
	print STCH "   		fi\n";
	
	print STCH "	else\n"; #blastx.parsed file exists
	print STCH "		".$run_script_path."check_Blast_parsed_file.pl \${BlastXIN}  \${PARSED}\n"; #check if the blastx.parsed file completed
	print STCH '		CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '		if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "		then\n";
	print STCH "			".$run_script_path."BLASTx_NR_parser.pl.sqlite \${BX_DIR} \${BlastXOUT}\n"; #run parser
	print STCH "			".$run_script_path."check_Blast_parsed_file.pl \${BlastXIN}  \${PARSED}\n";
	print STCH '			CHECK=$?',"\n";
	#check if parsed file is completed. If not completed, exit.
	print STCH '			if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "			then\n";
	print STCH "				echo \"Fatal Error trying to process \${MegaBlastOUT}.\"  \n";
	print STCH "				exit 1\n";
	print STCH "   			fi\n";
	print STCH "   		fi\n";
	print STCH "	fi\n";
	print STCH "fi";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}


#####################################################################################
sub generate_assignment_report{
	$current_job_file1 = "j35_".$sample_name."_generate_report.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=10G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load bio-perl perl-modules\n";
	print STCH "set -x\n";
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "REPORT=".$sample_dir."/".$sample_name.".AssignmentReport\n";
	print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".AssignmentReport_time.txt\n";
	print STCH "QC_Record=".$sample_dir."/".$sample_name.".BLASTX_NR_check.txt\n\n";

	print STCH "if [ ! -e $status_log/j35_AssignmentReport_finished ]\n"; # job never finished
	print STCH "then\n";
	print STCH "	date > \${TIMEFILE}\n";

	# check to make sure all parser finished 

	print STCH "	".$run_script_path."check_BLASTX_NR_Finished_correctly.pl \${SAMPLE_DIR} > \${QC_Record}  \n";
	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : BLASTX_NR did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";

	print STCH "	".$run_script_path."Generate_assignment_report.pl ".$sample_dir." \${INPUT}\n";
	print STCH '	grep "# Finished Assignment Report" ${REPORT}',"\n";
	print STCH '	CHECK=$?',"\n";
	print STCH '	if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error in generating report for \${SAMPLE_DIR}.\"  \n";
	print STCH "		exit 1\n";
	print STCH "   	fi\n";
	print STCH "fi\n";
	print STCH "touch $status_log/j35_AssignmentReport_finished  \n";
	print STCH "date >> \${TIMEFILE}\n";
	close STCH;
	
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
sub generate_assignment_summary {
	$current_job_file1 = "j36_".$sample_name."_generate_summary.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=10G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}

	print STCH "OUTPUT=".$sample_dir."/".$sample_name.".AssignmentSummary\n";
	print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".AssignmentSummary_time.txt\n";

	print STCH "if [ ! -e $status_log/j36_AssignmentSummary_finished ]\n"; # job never finished
	print STCH "then\n";
	print STCH "	date > \${TIMEFILE}\n";
	print STCH "	".$run_script_path."Generate_assignment_summary.pl ".$sample_dir." \${BAD_SEQ}\n";
 	print STCH "	OUT=\$?\n"; 
	print STCH '	if [ ${OUT} -ne 0 ]',"\n"; # did not finish successfully
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error : Generate_assignment_summary did not finish correctly.\"  \n";
	print STCH "		exit 1\n"; 
	print STCH "	fi\n";
	print STCH "fi\n";
	print STCH "date >> \${TIMEFILE}\n";
	print STCH "touch $status_log/j36_AssignmentSummary_finished  \n";
	close STCH;

	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
sub generate_phage_report{
	$current_job_file1 = "j37_".$sample_name."_generate_phage_report.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=6G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}
	print STCH "module load perl-modules\n";
	print STCH "module load bio-perl\n";
	print STCH "set -x\n";
	print STCH "SAMPLE_DIR=".$sample_dir."\n";
	print STCH "REPORT=".$sample_dir."/".$sample_name.".PhageAssignmentReport\n";
	print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".PhageAssignmentReport_time.txt\n";

	print STCH "if [ ! -e $status_log/j37_PhageAssignmentReport_finished ]\n"; # job never finished
	print STCH "then\n";
	print STCH "	date > \${TIMEFILE}\n";
	print STCH "	".$run_script_path."Phage_Generate_assignment_report.pl ".$sample_dir."\n";
	print STCH '	grep "# Finished Assignment Report" ${REPORT}',"\n";
	print STCH '	CHECK=$?',"\n";
	print STCH '	if [ ${CHECK} -ne 0 ]',"\n";   
	print STCH "	then\n";
	print STCH "		echo \"Fatal Error in generating phage report for \${SAMPLE_DIR}.\"  \n";
	print STCH "		exit 1\n";
	print STCH "   	fi\n";
	print STCH "fi\n";
	print STCH "touch $status_log/j37_PhageAssignmentReport_finished \n";
	print STCH "date >> \${TIMEFILE}\n";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}

#####################################################################################
sub generate_phage_summary {
	$current_job_file1 = "j38_".$sample_name."_generate_phage_summary.sh";
	open(STCH, ">$job_files_dir/$current_job_file1") or die $!;
	print STCH "#!/bin/bash\n";  
	print STCH "#SBATCH --mem-per-cpu=2G\n";
	print STCH "#SBATCH --output=".$SLURM_files_dir."/".$current_job_file1.".out\n";
	if (!$step_number) {
		print STCH "#SBATCH --depend=afterok:$last_jobid\n";
	}

	print STCH "OUTPUT=".$sample_dir."/".$sample_name.".PhageAssignmentSummary\n";
	print STCH "TIMEFILE=$sample_dir"."/".$sample_name.".PhageAssignmentSummary_time.txt\n";

	print STCH "if [ ! -e $status_log/j38_PhageAssignmentSummary_finished ]\n"; # job never finished
	print STCH "then\n";
	print STCH "	date > \${TIMEFILE}\n";
	print STCH '	CHECK=1',"\n";
	print STCH '	while [ ${CHECK} -ne 0 ] ',"\n"; 
	print STCH "	do\n";
	print STCH "		".$run_script_path."Phage_Generate_assignment_summary.pl ".$sample_dir."\n";
	print STCH '		grep "# Finished Assignment Summary" ${OUTPUT}',"\n";
	print STCH '		CHECK=$?',"\n";
	print STCH "	done\n";
	print STCH "fi\n";
	print STCH "date >> \${TIMEFILE}\n";
	print STCH "	touch $status_log/j38_PhageAssignmentSummary_finished  \n";
	close STCH;
	$last_jobid = `sbatch $job_files_dir/$current_job_file1 | awk '{ print \$NF }'`;
}





