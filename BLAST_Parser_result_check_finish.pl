
#----------------------------------------------------------#
# Copyright (C) 2005 Washington University, St. Louis, MO. #
# All Rights Reserved.                                     #
#                                                          #
# Author: Guoyan Zhao                                        #
# Send all comments to gzhao@ural.wustl.edu                #
#                                                          #
# DISCLAIMER: THIS SOFTWARE IS PROVIDED "AS IS"            #
#             WITHOUT WARRANTY OF ANY KIND.                #
#----------------------------------------------------------#

use strict;
#use constant TOOLS=> "/home2/gzhao/tools/";
my $usage = '
usage: script <directory> 
<directory> = dir
';
die $usage unless scalar @ARGV  == 1;

my ( $dir) = @ARGV;

opendir (DH, $dir) or die "can not open $dir!\n";

foreach my $file ( readdir DH ){
	if ($file =~ /parsed$/) {
		my $file_path = $dir."/".$file;	
		open(File, "<$file_path") or die "can not open file $file_path!\n";	
		my $found = 0;
		while (my $line = <File>) {
			if ($line =~ /Summary/) {
				$found = 1;
				last;
			}   
		}
		if (!$found){
			print " $file is not finished!\n";
		}
	}
}

exit;

