################################################################################
##  Author:
##      Guoyan Zhao <gzhao@pathology.wustl.edu>
#*
#* Copyright (C) Washington University in St. Louis Developed by Guoyan Zhao,
#* Washington University in St. Louis, Department of Pathology and Immunology.
#*
###############################################################################

#!/usr/bin/perl
use strict;
use warnings;

my $usage = "
This script will check all .blastn.parsed or .blastx.parsed files 
 to make sure parsing blast output file is finished for each file.

perl $0  <blast input file> <blast parsed file>
";

exit( 10 ) unless scalar @ARGV == 2;
my ( $fa_file, $PARSED ) = @ARGV;
my $HOME = $ENV{HOME};

my $finished = &check_blastnParsed_output($fa_file, $PARSED);

exit ($finished);
################################################################
sub check_blastnParsed_output {
	my ($fa_file,  $in_file ) = @_;
	my $have_summary_line = 0;
	my $line_count = 0;
	my $saved_seq = 0;
	my $total_seq = 0;
	my $num_undefined_taxon = 0;
	my $num_other_taxon = 0; 

	# total number of sequences in the input fasta file
	my $total_input_seq = `grep \">\" $fa_file |wc -l`;

	open (TEMP, "<$in_file") or return 10;
	while (my $line = <TEMP>) {
		$line_count++;
		if ($line =~ /# Summary: (\d+) out of (\d+)/) {
			$saved_seq = $1; 
			$total_seq = $2; # total record the parser processed
			$have_summary_line = 1;
		}
		if ($line =~ /undefined taxon/) {
			$num_undefined_taxon++;
		}
		if ($line =~ /other/) {
#			if (!($line =~ /cellular organisms|Viruses|other sequences|unclassified sequences/)) { # sometimes the lineage retrieval did not work and the sequence was assigned to "other"
			if ($line =~ /;{2,}/) { # sometimes the lineage retrieval did not work and the sequence was assigned to "other"

				print "sequence assigned to \"other\" but lineage is not correctively retrieved.\n";
				return 10;;
			}
		}
	
	}
	close TEMP;

	if (!$have_summary_line) {
		print "no summary line\n";
		return 10;
	}

	# check whether total input sequence is equal to total sequence been parsed
	if ($total_input_seq > $total_seq ) {
		print "BLAST output is not complete\n";
		return 10;
	}


	# taxonomy record has to be equal or greater than the number of sequences get 
	# successful phylotyped because some sequence could be assigned multiple taxonomy
	# categories. Should have at least $num_phylotyped + 1 lines
	my $num_phylotyped = $total_seq - $saved_seq;	
	if ( $num_phylotyped == 0 ) { # every sequence is unassigned
		print "num_phylotyped = $num_phylotyped, every sequence is unassigned\n";
		return 0;
	}
	# deal with situation where all records showed as undefined taxon and relative 
	# to number of phylotyped sequences
	elsif ( $num_phylotyped <= $num_undefined_taxon) { 
		print "num_undefined_taxon = $num_undefined_taxon, every sequence is undefined taxon\n";
		return 10; #changed from 0 to 10, the system default $? is 0, avoid the same value, the same reason below
	}

	if ( ($line_count - 1) == $num_undefined_taxon) { # deal with situation where all records showed as undefined taxon
		print "every sequence is un defined taxon\n";
		return 10;
	}

	# deal with old situation where some reads were not recorded because of no 
	# entry of gi-taxon record in the database 
	if ($num_phylotyped > ($line_count -1 ) ) {
		print "record number less than num phylotyped\n";
		return 10;
	}


	print "everything correct.\n";
	return 0; # use 0 to be consistent with Linux convention
}

