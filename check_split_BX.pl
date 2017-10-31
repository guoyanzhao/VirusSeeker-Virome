#!/usr/bin/perl -w
use strict;
my $usage='
perl script <sample dir> <sample name>
<sample dir> = full path of the folder holding files for this sample
               without last "/"
<sample name> = sample name

';
die $usage unless scalar @ARGV == 2;
my ( $dir, $sample_name ) = @ARGV;

if ($dir =~/(.+)\/$/) {
        $dir = $1;
}

my $finished = &check_split_output($dir);
#print $finished; 
exit ($finished);

########################################################################
sub check_split_output {
	my ( $dir ) = @_;
	my $tot_in_seq=0;
	my $tot_out_seq = 0;

	my $BNfiltered_File = $dir."/".$sample_name."_BLASTN_NT_Filtered.fa";
	$tot_in_seq = `grep \">\" $BNfiltered_File |wc -l`;

	my $BLASTX_dir = $dir."/".$sample_name."_BLASTX_NR";
	opendir(DH, $BLASTX_dir) or return 10;
	foreach my $file (readdir DH) {
		if ($file =~ /\.fasta$/ ) { 
			my $faFile = $BLASTX_dir."/".$file;
			my $temp = `grep \">\" $faFile |wc -l`;
			$tot_out_seq += $temp;                    
		}
	}
	close DH; 

	print "total input sequence: $tot_in_seq\n";
	print "total output sequence : $tot_out_seq\n"; 

	if($tot_in_seq==$tot_out_seq) { return 0; }  # use 0 to make it consistent with Linux convention
	else { # remove content of BLASTX directory, to rerun split 
		opendir(DH, $BLASTX_dir) or return 10;
		foreach my $file (readdir DH) {	
			if (!($file eq ".") && !($file eq "..")) {
				my $faFile = $BLASTX_dir."/".$file;
				# print "$faFile\n"; 
				unlink $faFile; 
			}
		}
	  	close DH;  
	  	return 10; 
	}
}
