#!/usr/bin/env perl

use warnings;
use strict;
use autodie;
use File::Basename;
use IO::Zlib;

my $inputFile = $ARGV[0];
my $inputBase = basename($ARGV[0]);
my $numberOfHeaderLines = $ARGV[1]
my $gtReport = new IO::Zlib;
$gtReport->open($inputFile, 'rb');

my $outFileName = $inputBase =~ s/.csv/.lgen/gr;
open(outputFile,">",$outFileName);

my $lineNumber = 0;

while(my $line = <$gtReport>) 
{
  next if $lineNumber++ < $numberOfHeaderLines+1;
  my @column = split (/,/,$line);
  my $selectedColumns = "$column[1] $column[1] $column[0] $column[2] $column[3]\n";
  print outputFile $selectedColumns;
}
$gtReport->close;
close(outputFile);
