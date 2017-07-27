#!/usr/bin/perl -w
use strict;
(my $usage = <<OUT) =~ s/\t+//g;
perl $0 <sample dir>
<sample dir> = full path of the folder holding files for this sample without last "/"

OUT

die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;

if ($dir =~/(.+)\/$/) {
        $dir = $1;
}
my $sample_name = (split(/\//,$dir))[-1];
#print $sample_name,"\n";

my $finished = &check_split_output($dir);
#print $finished; 
exit ($finished);

##############################################################
sub check_split_output {
	my ( $dir ) = @_;
	my $tot_in_seq=0;
	my $tot_out_seq = 0;

	my $in_file = $dir."/".$sample_name.".QCed.cdhit.tantan.goodSeq.fa";
	$tot_in_seq = `grep \">\" $in_file |wc -l`;

	my $RepeatMasker_dir = $dir."/".$sample_name."_RepeatMasker";
	print "output dir is \\$RepeatMasker_dir\\\n";
	opendir(DH, $RepeatMasker_dir) or return 10;
	foreach my $file (readdir DH) {
		if ($file =~ /\.fasta$/ ) { 
			my $faFile = $RepeatMasker_dir."/".$file;
			my $temp = `grep \">\" $faFile |wc -l`;
			$tot_out_seq += $temp;                    
		}
	}
	close DH; 

	print "total input sequence: $tot_in_seq\n";
	print "total output sequence : $tot_out_seq\n"; 

	if($tot_in_seq==$tot_out_seq) { return 0; }  # use 0 to be consistant with Linux convention
	else { # remove content of RepeatMasker directory, to rerun split 
		opendir(DH, $RepeatMasker_dir) or return 10;
		foreach my $file (readdir DH) {	
			if (!($file eq ".") && !($file eq "..")) {
				my $faFile = $RepeatMasker_dir."/".$file;
#				print "to remove $faFile\n"; 
				unlink $faFile; 
			}
		}
	  	close DH;  
	  	return 10; 
	}
}
