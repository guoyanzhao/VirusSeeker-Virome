#!/usr/bin/perl -w
use strict;
my $usage='
perl script <sample dir>
<sample dir> = full path of the folder holding files for a sample

';
die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;

if ($dir =~/(.+)\/$/) {
	$dir = $1;
}
my $sample_name = (split(/\//,$dir))[-1];
#print $sample_name,"\n";

my $finished = &check_split($dir);
#print $finished; 
exit ($finished);

##########################################################################
sub check_split {
	my ( $dir ) = @_;
	my $tot_in_seq=0;
	my $tot_out_seq = 0;

	my $in_file = $dir."/".$sample_name.".BLASTN_VIRUSDB_Filtered.fa";
#	print "in file is $in_file\n";
	$tot_in_seq = `grep \">\" $in_file |wc -l`;

	my $BLASTX_VIRUSDB_dir = $dir."/".$sample_name."_BLASTX_VIRUSDB";
	opendir(DH, $BLASTX_VIRUSDB_dir) or return 10;
	foreach my $file (readdir DH) {
		if ($file =~ /\.fasta$/ ) {
			my $faFile = $BLASTX_VIRUSDB_dir."/".$file;
			my $temp = `grep \">\" $faFile |wc -l`;
			$tot_out_seq += $temp;
		}
	}
	close DH;

	print "total seq in input file: $tot_in_seq\n";
	print "total sequence in splited output dir: $tot_out_seq\n\n"; 

	if($tot_in_seq==$tot_out_seq) { return 0; } 
	else { 
		return 10; 
	}
}
