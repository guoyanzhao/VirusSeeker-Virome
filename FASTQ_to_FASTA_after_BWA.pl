#########################################################
# in this version the read name ID is trimmed for "/2" for the purpose 
# to fit the pipeline because I used bamToFastq to generate the 
# fastq from bam mapping file and the read ID were postfixed with "/2" 
#!/usr/bin/perl -w

use strict;

my $usage = '
Given a fastq file, this script will conver it to a fasta format file.

perl script <fastq file>
<fastq file> = full path to the fastq file

output to standard output.

';
die $usage unless scalar @ARGV == 1;
my ($inFile) = @ARGV;

open (IN, "<$inFile") or die "Can't Open file: $inFile";

while (my $line = <IN>){
	my $info = $line; # first line
	chomp $line;
	my $name = $line;
	$name =~ s/@//;
	print ">$name\n";
	$line = <IN>; # second line: sequence
	print "$line";
	$line = <IN>; # third line: + 
	$line = <IN>; # fourth line: quality line
}

close IN;

exit;

