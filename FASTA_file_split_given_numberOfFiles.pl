
#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $usage = "
perl $0 <-d input dir><-i input fasta file><-o output dir><-f number of files to split into><-n minimum number of reads per file><-h>
<-d> = full path to the directory where the fasta file reside without last \"/\"
<-i> = name of the fasta file
<-o> = output directory path
<-f> = number of files to split into
<-n> = minimum number of reads in each file 

Given a fasta file, number of files to be splited into -f, this script will split 
input fasta file to maximum -f number of files. Each file will contain minimum -n 
number of sequences. 

Generated files have the same name as 
the given file with numbered suffix .file0.fa .file1.fa ... etc.  

";

my %opts;
GetOptions(\%opts, "d:s", "i:s", "o:s", "f:i", "n=i", "h");

die $usage if(defined($opts{h}) || !defined($opts{d}) || !defined($opts{i}) || !defined($opts{o}) || !defined($opts{f}) || !defined($opts{n}) );
$| = 1;

my $dir = $opts{d};
my $fasta = $opts{i};
my $out_dir = $opts{o};
my $numFile = $opts{f};
my $minSeq = $opts{n}; 

my $inFile = $dir."/".$fasta;
print "input is $inFile\n";
my $oldSeperator = $/;

# get total number of sequence in input file
my $total_seq = `grep ">" $inFile |wc -l`;

# calculate how many sequences in each file
my $size = int ($total_seq/$opts{f}) + 1;
#print "$count_seq $size\n";

my  $numSeq = 0; # number of sequences to be in each file 
if ($size > $minSeq) {
	$numSeq = $size; # has to have $size number of sequence each file in order to keep the number of files <= $numFile
}
else {
	$numSeq =  $minSeq;
}

# start spliting
my $geneCount = 0;
my $fileCount = 1;
my $geneCount_in_last_file = 0;
my $last_file = ""; 
my $outFile = $out_dir."/".$fasta."_file".$fileCount.".fasta";
open (OUT, ">$outFile") or die "can not open file $outFile!\n";

open (IN, "<$inFile") or die "Can't Open file: $inFile";
$/=">";
<IN>; # remove the first ">"
while (my $line = <IN>){
#	print "$line\n";
	chomp $line;
	print OUT ">", $line; # first read

	$geneCount++;
	$geneCount_in_last_file++;

	if (!($geneCount%$numSeq)) {
		close OUT;
		$fileCount++;
		$geneCount_in_last_file = 0;
		$outFile = $out_dir."/".$fasta."_file".$fileCount.".fasta";
		open (OUT, ">$outFile") or die "can not open file $outFile!\n";
		$last_file = $outFile;
	}
}

close IN;
close OUT;

if (!$geneCount_in_last_file) {
	unlink $last_file;
}

$/ = $oldSeperator;

exit;

