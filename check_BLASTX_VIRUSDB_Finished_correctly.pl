#!/usr/bin/perl -w
use strict;
my $usage='
perl script <sample dir>
<sample dir> = full path of the folder holding files for this sample
               without last "/"
';
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
	my @file_not_finished_blast = ();
	my @file_not_finished_parsed = ();

	# total input sequence
	my $BLAST_dir = $dir."/".$sample_name."_BLASTX_VIRUSDB";
	opendir(DH, $BLAST_dir) or return 10;
	foreach my $file (readdir DH) {
		if ($file =~ /\.fasta$/ ) { 
			my $faFile = $BLAST_dir."/".$file;
			my $temp = `grep \">\" $faFile |wc -l`;
			$tot_in_seq += $temp;                    

			# get parsed file not finished
			my $file_name = $file;
			$file_name =~ s/fasta$/BLASTX_VIRUSDB\.parsed/;
			my $parsedFile = $BLAST_dir."/".$file_name;
#			print "file name is $parsedFile\n";
			if (-e $parsedFile ) {# have parsed file, check for completeness
				my $temp = `grep \"Summary\" $parsedFile `;
				if (!($temp =~ /Summary/)) {
					push @file_not_finished_parsed, $file;
				}
			}
			else { # parsed file does not exist
				push @file_not_finished_parsed, $file;
			}

			# get blast output file not finished
			$file_name = $file;
			$file_name =~ s/fasta$/BLASTX_VIRUSDB\.out/;
			my $BlastOutFile = $BLAST_dir."/".$file_name;
#			print "file name is $BlastOutFile\n";
			if (-e $BlastOutFile ) {# have out put file, check for completeness
				my $temp = `grep \"Matrix:\" $BlastOutFile `;
				if (!($temp =~ /Matrix:/)) {
					push @file_not_finished_blast, $file;
				}
			}
			else { # blast output file does not exist
				push @file_not_finished_blast, $file;
			}

		}
	}
	close DH; 

	# total processed sequence
	opendir(DH, $BLAST_dir) or return 10;
	foreach my $file (readdir DH) {
		if ($file =~ /\.fa$/ ) { 
			my $faFile = $BLAST_dir."/".$file;
			my $temp = `grep \">\" $faFile |wc -l`;
			$tot_out_seq += $temp;                    
		}
	}
	close DH; 

	print "total input sequence: $tot_in_seq\n";
	print "total output sequence : $tot_out_seq\n\n"; 

	if($tot_in_seq==$tot_out_seq) { return 0; } 
	else {
 		my @sorted = sort (@file_not_finished_blast);
		print "following files did not finish blast : \n";
		foreach my $file ( @sorted) {
			print "$file\n";
		}
		print"\n";
 
 		@sorted = sort (@file_not_finished_parsed);
		print "following files parsed file is not correct : \n";
		foreach my $file ( @sorted) {
			print "$file\n";
		}
 
	  	return 10; 
	}
}
