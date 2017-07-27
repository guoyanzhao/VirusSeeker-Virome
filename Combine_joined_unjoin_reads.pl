################################################################################
##  Author:
##      Guoyan Zhao <gzhao@pathology.wustl.edu>
#*
#* Copyright (C) Washington University in St. Louis Developed by Guoyan Zhao,
#* Washington University in St. Louis, Department of Pathology and Immunology.
#*
###############################################################################

#!/usr/bin/perl -w

use strict;

my $usage = '
Given the fastq-join output un1, un2 and join fastq files, the script will
rename the id in the sequencing data. Joined reads will have a prefix of
"join", un1 reads (unjoined SE1 reads) will have a prefix of "un1" and un2 
reads will have a prefix of "un2".

perl script < full path and the fastq-join output file prefix for un1, un2 
and join fastq files>

';

die $usage unless scalar @ARGV == 1;
my ($inName) = @ARGV;

my $file_join=$inName.".join.fq";
my $file_un1=$inName.".un1.fq"; 
my $file_un2=$inName.".un2.fq";  
my $file_out=$inName.".stitched.fastq";

open (INjoin, "<$file_join") or die "Can't Open file: $file_join";
open (INun1,"<$file_un1") or die "Can't Open file: $file_un1";
open (INun2,"<$file_un2") or die "Can't Open file: $file_un2";
open (OUT,">$file_out") or die "Can't Open file: $file_out";

my $cc=0;
while (my $line = <INjoin>){
	$cc++; 
	if ($cc-4*int($cc/4)==1) {
		#print $line; <STDIN>; 
		my $info = $line;
		chomp $line;
		my @name1 = split(/\s+/,$line);  # get rid of the mate info 1:N:0:1
		my @name2 = split(/\:/, $name1[0]); 
		my $name = $name2[3];
		for(my $i=4; $i <= $#name2; $i++){ 
			$name=join("\:", $name, $name2[$i]); 
		}
 		print OUT "\@join:$name\n";
	}
	else { print OUT $line; }
}

$cc=0;
while (my $line = <INun1>){
	$cc++;
	if ($cc-4*int($cc/4)==1) {
		my $info = $line;
		chomp $line;
		my @name1 = split(/\s+/,$line);  # get rid of the mate info 1:N:0:1
		my @name2 = split(/\:/, $name1[0]); 
		my $name = $name2[3];
		for(my $i=4; $i <= $#name2; $i++){ 
			$name=join("\:", $name, $name2[$i]); 
		}
		print OUT "\@un1:$name\n";
	}
	else { print OUT $line; }
}

$cc=0; 
while (my $line = <INun2>){
	$cc++;
	if ($cc-4*int($cc/4)==1) {
		my $info = $line;
		chomp $line;
		my @name1 = split(/\s+/,$line);  # get rid of the mate info 1:N:0:1
		my @name2 = split(/\:/, $name1[0]); 
		my $name = $name2[3];
		for(my $i=4; $i <= $#name2; $i++){ 
			$name=join("\:", $name, $name2[$i]); 
		}
		print OUT "\@un2:$name\n";
	}
	else { print OUT $line; }
}

close INjoin; 
close INun2; 
close INun1; 

exit;

