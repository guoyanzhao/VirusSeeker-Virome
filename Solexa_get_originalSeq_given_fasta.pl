
#!/usr/bin/perl -w

use strict;
use warnings;

my $usage = '
This script will accept a fasta file, get the original sequences for those reads, 
output to an output file.

perl script <input FASTAQ file> <fasta file>

<input FASTAQ file1> = full path to FASTAQ file with original data
<fasta file> = full path to the the fasta file

output to a FASTAQ files with the same name except end with .fastq in the same 
directory as the input FASTAQ files
';

die $usage unless scalar @ARGV == 2;
my ($in_FASTAQ, $unmapped_reads_file) = @ARGV;

my $outfile = $unmapped_reads_file."stq";
#print "out file is $outfile\n";
open (OUT, ">$outfile") or die $!;

my %read_names = (); # name of the read => 1
open (IN, "<$unmapped_reads_file") or die $!;
while (my $line = <IN>) {
	if ($line =~ />/) {
		chomp $line;
		my @fields = split (/\t/, $line);
		my $read_name = $fields[0];
		$read_name =~ s/>//;
#		print "read = \\$read_name\\\n";
		$read_names{$read_name} = 1;
	}
}
close IN;

my $keys = keys %read_names;
print "total  reads = $keys\n\n";

open(IN, $in_FASTAQ)||die $!;
while (my $line = <IN>) {
#	print "line = \\$line\\\n";
	chomp $line; # the first line with @seq-name
	my @temp = split(/\s+/, $line);
	my $temp_seq_name = $temp[0];
	my $seq_name = substr($temp_seq_name, 1);
#	print ">$seq_name<\n";
	if ( exists $read_names{$seq_name}) {
		print OUT "@", $seq_name, "\n";
		$line = <IN>; # seq line
		print OUT $line;
		$line = <IN>; # + line
		print OUT $line;
		$line = <IN>; # quality line
		print OUT $line;

	}
}
close IN;
close OUT;

exit;

