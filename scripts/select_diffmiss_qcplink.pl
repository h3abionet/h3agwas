#!/usr/bin/env perl

use strict;

open IN, '<', $ARGV[1] or die "Cannot open missing file \n";
open OUT, '>', $ARGV[2];
while(<IN>){
	s/^\s+//;
	my @fields = split /\s+/, $_;
	unless($fields[0] eq 'CHR'){
		if($fields[4] < $ARGV[0]){
			print OUT "$fields[1]\n";
		}
	}
}
