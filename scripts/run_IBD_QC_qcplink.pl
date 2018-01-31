#!/usr/bin/env perl

use strict;

my %imiss;
my %removed;

open IMISS, '<', $ARGV[0]
    or die "Cannot open missing file (".$ARGV[0]."): $!\n";

print "Reading PLINK .imiss file ".$ARGV[0]."\n";

while(<IMISS>){
    s/^\s+//;
    my @fields = split /\s+/, $_;
    $imiss{$fields[0]}{$fields[1]} = $fields[5];
}

open GENOME, '<', $ARGV[1]
    or die "Cannot open genotypes file (".$ARGV[1]."): $!\n";
open OUT, '>', $ARGV[2];
print "Reading PLINK .genome file ".$ARGV[1]."\n";
while(<GENOME>){
    s/^\s+//;
    my @fields = split /\s+/, $_;
    if($fields[9] > 0.185){
	if($imiss{$fields[0]}{$fields[1]}>$imiss{$fields[2]}{$fields[3]}){
	    unless($removed{$fields[0]}{$fields[1]}){
		print OUT "$fields[0] $fields[1]\n";
		$removed{$fields[0]}{$fields[1]} = 1;
	    }
	}
	elsif($imiss{$fields[0]}{$fields[1]}<$imiss{$fields[2]}{$fields[3]}){
	    unless($removed{$fields[2]}{$fields[3]}){
		print OUT "$fields[2] $fields[3]\n";
		$removed{$fields[2]}{$fields[3]} = 1;
	    }
	}
	else{
	    unless($removed{$fields[0]}{$fields[1]}){
		print OUT "$fields[0] $fields[1]\n";
		$removed{$fields[0]}{$fields[1]} = 1;
	    }
	}
    }
}
