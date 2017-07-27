In this version:

Modified BLASTn and BLASTx parser to bin sequences that hit both
virus and other species to "Ambiguous" bin. Sequences only hit virus with
significant e value are counted as true "viral".
Finished modification 2014/02/04.
Updated to keep sequences that hit both virus and vector sequences.

Modified BLASTn nt, BLASTx NR and BLASTX VIRUSDB parser to generate
best hit alignment information in the parsed file. 

Replace Soap mapping to human genome with BWA-mem
Use CD-HIT 100% ID with 100% coverage for shorter sequence.
Adapted to htcf using The Simple Linux Utility for Resource Management (Slurm) 
as the job scheduling system.

3/23/2015
TANTAN QC step, do not classify sequences into tantan.RepeatLowComplexSeq
All those sequences will be kept in the analysis

V0.06
3/25/2015
After RepeatMasker, use original reads for mapping to host genome instead of
using masked.
After mapping to host genome use masked reads for MegaBLAST.
After pooling BLASTn VIRUSDB and BLASTx VIRUSDB hits, use original (Unmasked) reads 
for mapping to bacteria genome
After mapping to bacteria genome use masked reads for following BLAST

4/30/2015
MegaBLAST NT e value cutoff change to 1e-10

V0.061
10/08/2015
modified following scripts to check completeness at the corresponding steps. 
check to make sure all blast input files have corresponding output files and 
the output file finished correctly. If not will report which files did not
finish correctly.

check_BLAST_NT_Finished_correctly.pl
check_BLASTN_VIRUSDB_Finished_correctly.pl
check_BLASTX_NR_Finished_correctly.pl
check_BLASTX_VIRUSDB_Finished_correctly.pl
check_MegaBLAST_HOST_Finished_correctly.pl
check_MegaBLAST_NT_Finished_correctly.pl

modified following file to make sure the blast input file query sequence number
match the number of records in the parsed summary line. If the blast output file
did not fihish you can have a summary line but the total records is less than the 
total number of input sequences in the input file. The modified script will check
to make sure it's correct.
check_Blast_parsed_file.pl

modified following scripts so that if a seqence hit both virus and phage it will 
be classified as phage instead of "Ambiguous" as in the old version.
BLASTn_VIRUSDB_parser.pl
BLASTx_VIRUSDB_parser.pl

modified following script so that sequences assigned to phiX174 are saved to a 
seperate file *phiX174.fasta . The .PhageReads_All.fa file contains all the phage
sequences obtain by all blast steps except the phiX174 sequences.
Phage_Generate_assignment_report.pl

# V0.062 06/01/2016
use SQLite to replace MySQL database 
removed all while loop

# 2016/10/17
add extra sequences at FASTQ_to_FASTA.pl step to the file ".RemoveAdapter.stitched.prinseq.fasta"
 to prevent premature stop of pipeline
because no sequence left after sequencial blast analysis
e are unassigned sequence from I3886_17803_Project618_Stool_Botswana_Kwon_Stool_2_XE28202
print ">un2:1:2102:16457:16451\n";
print "ATGAGTGTCTAGTAATATTCTGTAATTACTTCTTAGATGAATATTTTACTATAGCATTATATTCATCTAAGAAGTAATTACAGAATATATCTAGAGTTTCATTAGAGAATCTACGAGGAATTATAGTGGCTTCCATAGAACCTTCCGCATTAGGAGAGTACTTTATAGTGCCTAAAAAGTCTTTTGGTAATTTAATACATTCTATACTATACCAAGTCTTTTTATCGGACATA\n";
print ">join:1:1101:8673:17533\n";
print "CGAATGCCTTCTAAGACAGATAGCTGTCTATAGTATTATCGCTATCTATTTCCAGCAACTCATTAAACAAGTCTTTAGCCAATGCTTTCCACTGCTCTCCCCAATCACGGAGATTCTCGACCTTTGACTGTATGTCTTCGAAATAAGAATCTACGTCTGATTTGATTGATTTTGAATAGTATTTAACATCCTTCTCATCCCCATCCATCATATAATCACATTGTGCCCTGATATCTTTTATATGACTGTCTATATCACTTCACATATAATCAACAGGTTTACGTATATTGAATATAGCTTCTGACGTAAGACCGGTTATATCTTGTATGTCTTTTAAATTACCCATGATTTAATCAATTAAATGCCAACCATCCCCCCTGCGC\n";

# 2016/10/19
Added steps to generate Phage summary files to the fully automated pipeline.

# 2016/10/24
Change to use bamToFastq from Hydra pacakge to bedtools package. It's from the same lab
Aaron Quinlan
http://quinlanlab.org/

# 2016/11/16
changed taxonomy temp directory to a common directory  /scratch/dwlab/tmp_taxonomy
