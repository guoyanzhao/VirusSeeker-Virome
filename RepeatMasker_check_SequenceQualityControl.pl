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

my $finished = &check_QC_read_number($dir);
#print $finished; 
exit ($finished);

##########################################################################
sub check_QC_read_number {
	my ( $dir ) = @_;
	my $tot_in_seq=0;
	my $tot_out_seq = 0;


	my $in_file = $dir."/".$sample_name.".QCed.cdhit.tantan.goodSeq.fa";
	$tot_in_seq = `grep \">\" $in_file |wc -l`;

	print "in file is $in_file\n";
	print "total sequence in input file: $tot_in_seq\n";

	my $badSeq_file = $dir."/".$sample_name.".RepeatMasker.badSeq.fa";;
	my $temp = `grep \"total seq\" $badSeq_file`;
	if ($temp =~ /total seq = (\d+)/) {
		$tot_out_seq = $1;
	}

	print "RepeatMasker QC output file: $badSeq_file\n";
	print "total sequence in RepeatMasker QC output file: $tot_out_seq\n\n";

	if($tot_in_seq==$tot_out_seq) { return 0; } # use 0 to be consistent with linux convention 
	else { 
		return 10; 
	}
}

