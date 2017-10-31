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

##############################################################
sub check_split_output {
	my ( $dir ) = @_;
	my $tot_in_seq=0;
	my $tot_out_seq = 0;

	my $in_file = $dir."/".$sample_name.".VIRUSDB_HIT_MegaBLAST_NT_Filtered.fa";
	$tot_in_seq = `grep \">\" $in_file |wc -l`;

	my $BLASTN_dir = $dir."/".$sample_name."_BLASTN_NT";
	opendir(DH, $BLASTN_dir) or return 10;
	foreach my $file (readdir DH) {
		if ($file =~ /\.fasta$/ ) { 
			my $faFile = $BLASTN_dir."/".$file;
			my $temp = `grep \">\" $faFile |wc -l`;
			$tot_out_seq += $temp;                    
		}
	}
	close DH; 

	print "total input sequence: $tot_in_seq\n";
	print "total output sequence : $tot_out_seq\n"; 

	if($tot_in_seq==$tot_out_seq) { return 0; } 
	else { # remove content of BLASTN directory, to rerun split 
		opendir(DH, $BLASTN_dir) or return 10;
		foreach my $file (readdir DH) {	
			if (!($file eq ".") && !($file eq "..")) {
				my $faFile = $BLASTN_dir."/".$file;
        	    # print "$faFile\n"; 
				unlink $faFile; 
			}
		}
	  	close DH;  
	  	return 10; 
	}
}
