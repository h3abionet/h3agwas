#!/usr/bin/env perl

use strict;

open IN, "<$missing" or die "Cannot open missing file \n";
open OUT, ">$failed";
while(<IN>){
	s/^\\s+//;
	my @fields = split /\\s+/, \$_;
	unless(\$fields[0] eq 'CHR'){
	    if(\$fields[$probcol] < "$cut_diff_miss") {
		print "\$fields[1] \$fields[4]\n";
			print OUT "\$fields[1]\n";
		}
	}
}
