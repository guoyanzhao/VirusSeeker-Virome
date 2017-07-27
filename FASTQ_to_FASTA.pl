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

# add extra sequence to avoid pipeline stopping
# these are unassigned sequence from I3886_17803_Project618_Stool_Botswana_Kwon_Stool_2_XE28202
print ">un2:1:2102:16457:16451\n";
print "ATGAGTGTCTAGTAATATTCTGTAATTACTTCTTAGATGAATATTTTACTATAGCATTATATTCATCTAAGAAGTAATTACAGAATATATCTAGAGTTTCATTAGAGAATCTACGAGGAATTATAGTGGCTTCCATAGAACCTTCCGCATTAGGAGAGTACTTTATAGTGCCTAAAAAGTCTTTTGGTAATTTAATACATTCTATACTATACCAAGTCTTTTTATCGGACATA\n";
print ">join:1:1101:8673:17533\n";
print "CGAATGCCTTCTAAGACAGATAGCTGTCTATAGTATTATCGCTATCTATTTCCAGCAACTCATTAAACAAGTCTTTAGCCAATGCTTTCCACTGCTCTCCCCAATCACGGAGATTCTCGACCTTTGACTGTATGTCTTCGAAATAAGAATCTACGTCTGATTTGATTGATTTTGAATAGTATTTAACATCCTTCTCATCCCCATCCATCATATAATCACATTGTGCCCTGATATCTTTTATATGACTGTCTATATCACTTCACATATAATCAACAGGTTTACGTATATTGAATATAGCTTCTGACGTAAGACCGGTTATATCTTGTATGTCTTTTAAATTACCCATGATTTAATCAATTAAATGCCAACCATCCCCCCTGCGC\n";
close IN;

exit;

