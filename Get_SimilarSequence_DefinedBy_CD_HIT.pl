
#!/usr/bin/perl -w

use strict;
use Tie::IxHash;

my $usage = '
Given a library directory, a file with unique sequences from this library and 
the number of similar sequences, this script will get read ID for all the reads 
that are similar to the given reads defined by CD-Hit and output given number 
of similar reads (maximum given number, if the cluster has less number of similar
reads it will output all the reads in the cluster) .

perl script 
<input total read file> = full path to the file with all reads 
<input unique file>     = full path to the file contains unique reads
<clustering file>       = cd-hit output .clstr file
<number of similar reads to output> = maximum number of similar reads 

output to standard output
';

die $usage unless @ARGV == 4;

my ($lib_fa, $unique_file, $clustering_file, $output_num_reads) = @ARGV;

# get all the reads
my %sequences = (); # all sequences in this library
&read_FASTA_data($lib_fa, \%sequences);

# get the unique sequences
my %unique_seq = (); 
tie %unique_seq, "Tie::IxHash";
&read_FASTA_data($unique_file, \%unique_seq );

# go to CD-HIT output file and get all the reads that are in the same
# cluster as the viral reads.
my %cluster = (); # read_ID => [read_IDs]
&read_CD_HIT_clstr_file($clustering_file, \%cluster);

foreach my $read (keys %unique_seq) {
	# output the representative sequence
	print ">$read\n";
	print $sequences{$read}, "\n";
	if (defined $cluster{$read}) { # have similar sequence
		my $cluster_member_hash_ref = $cluster{$read};

#		foreach my $key (sort {$cluster_member_hash_ref->{$b} <=> $cluster_member_hash_ref->{$a} } keys %{$cluster{$read}}) {
#			print "cluser: $read, $key, $cluster{$read}->{$key}\n";
#		}

		my $i = 1;
		foreach my $ele (sort {$cluster_member_hash_ref->{$b} <=> $cluster_member_hash_ref->{$a}} keys %{$cluster_member_hash_ref}) {
			if ($i <= $output_num_reads) {
				print ">$ele\n";
				print $sequences{$ele}, "\n";
				$i++;
			}
			else {
				last;
			}
		}
	}

#		print "******************************************\n\n";
}

exit;


######################################################
# in this version, sequence length is saved
sub read_CD_HIT_clstr_file {
	my ($in_file, $hash_ref) = @_;


	open(IN, "<$in_file") or die "can not open file $in_file!\n";
	my $content = do {local $/; <IN>}; # read the whole content
	close IN;
#	print "content is \\$content\\\n";
	while ($content =~ /(^>Cluster \d+\n(\d+.*\n)+)/gm) {
		my %seqID_to_length = (); # sequence ID => sequence length
		my $rep = "";

		my $temp = $1;
#		print "\nfind the match, \\$temp\\ \n";
		my @lines = split(/\n/, $temp);
		shift @lines; # remove >Cluster 0 etc 
		foreach my $line (@lines) {
#			print "line is : \\$line\\\n";
			my @temp = split(/\s+/, $line);
#			print "num = \\$temp[0]\\, len = \\$temp[1]\\, name = \\$temp[2]\\\n";
			my $len = $temp[1];
			$len =~ s/nt,//;
			my $name = $temp[2];
 			$name =~ s/>//;
			$name =~ s/\.\.\.//;
#			print "now length = $len, name = \\$name\\\n";
			if ($line =~ /\*/) {
				$rep = $name;
			}
			else {
				$seqID_to_length{$name} = $len;
			}
		}

=head1
		print "rep = $rep, \n";
		foreach my $key (sort {$seqID_to_length{$b} <=> $seqID_to_length{$a}} keys %seqID_to_length) {
			print "$key, $seqID_to_length{$key}\n";

		}
print "******************************************\n\n";
=cut
		$hash_ref->{$rep} = \%seqID_to_length;
	}
}

#####################################################################
sub read_FASTA_data () {
    my ($fastaFile, $hash_ref) =  @_;

    #keep old read seperator and set new read seperator to ">"
    my $oldseperator = $/;
    $/ = ">";

    open (FastaFile, $fastaFile) or die "Can't Open FASTA file: $fastaFile";
    while (my $line = <FastaFile>){
		# Discard blank lines
        if ($line =~ /^\s*$/) {
			next;
		}
		# discard comment lines
		elsif ($line =~ /^\s*#/) {
	       next;
		}
		# discard the first line which only has ">", keep the rest
		elsif ($line ne ">") {
		    chomp $line;
		    my @rows = ();
		    @rows = split (/\n/m, $line);	
		    my $firstLine = shift @rows;
			my @temp_arr = split(/\s/, $firstLine);
			my $Name = shift @temp_arr;
		    my $Seq = join("", @rows);
		    $Seq =~ s/\s//g; #remove white space
		    $hash_ref->{$Name} = $Seq;
#			print " name = \\$Name\\, seq  = \\$Seq\\\n\n";
		}
    }

    # execute right
#	 foreach my $key (keys %$hash_ref){
#	      print "Here is the key for fasta seq: $key \t $hash_ref->{$key}\n";
#	 }

    #reset the read seperator
    $/ = $oldseperator;
	close FastaFile;
}


