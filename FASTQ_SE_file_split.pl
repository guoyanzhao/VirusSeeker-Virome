
#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $usage = "
Given a fastq file this script will split it to a number of files. Each file will 
contain given number of sequences. Generated files have the same name as the given file 
with numbered suffix .file0.fa .file1.fa ... etc. This script make sure paired reads
remain paired after spliting.

Each sequence should occupy 4 lines: name ID, sequence, + or ID, quality score.

perl $0 <-d input dir><-i input fastq file><-o output dir><-n number of reads per file><-f max number of files><-h>
<-d> = full path to the directory where the fastq file reside without last \"/\"
<-i> = name of the fastq file
<-o> = output directory path
<-f> = maximum number of files. This is for SGE job submission.
<-n> = minimum number of reads in each file, 500000 for Newbler assembly

";

my %opts;
GetOptions(\%opts, "d:s", "i:s", "o:s", "n=i", "f=i", "h");

die $usage if(defined($opts{h}) || !defined($opts{d}) || !defined($opts{i}) || !defined($opts{o})  || !defined($opts{n}) || !defined($opts{f}));
$| = 1;

my $dir = $opts{d};
my $fastq = $opts{i};
my $out_dir = $opts{o};
my $numFile = $opts{f};
my $minSeq = $opts{n};

my $inFile = $dir."/".$fastq;
print "input is $inFile\n";

# get total number of sequence in input file
my $total_seq = `wc -l < $inFile `;
print "total line = $total_seq\n";
$total_seq = $total_seq/4; # because this is for SE sequence, each seq has 4 lines.
print "total # seq = $total_seq\n";

# calculate how many sequences in each file
my $size = int ($total_seq/$numFile) + 1;
print "$total_seq, # seq per file:  $size\n";

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
my $outFile = $out_dir."/".$fastq."_file".$fileCount.".fastq";
open (OUT, ">$outFile") or die "can not open file $outFile!\n";

open (IN, "<$inFile") or die "Can't Open file: $inFile";
while (my $line = <IN>){
	print OUT "$line"; # first line of first read
	$line = <IN>; 
	print OUT "$line"; # 2nd line of first read
	$line = <IN>; 
	print OUT "$line"; # 3rd line of first read
	$line = <IN>; 
	print OUT "$line"; # 4th line of first read

	$geneCount++;
	$geneCount_in_last_file++;

	if (!($geneCount%$numSeq)) {
		close OUT;
		$fileCount++;
		$geneCount_in_last_file = 0;
		$outFile = $out_dir."/".$fastq."_file".$fileCount.".fastq";
		open (OUT, ">$outFile") or die "can not open file $outFile!\n";
		$last_file = $outFile;
	}
}

close IN;

if (!($geneCount_in_last_file)) {
	unlink $last_file;
}

exit;

