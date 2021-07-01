#!/usr/bin/env perl

use warnings;
use strict;
use File::Basename;
use IO::Zlib;

my \$inputFile = "${genotypeReport}";
my \$numberOfHeaderLines = ${params.numberOfGtReportHeaderLines};
my \$gtReport = new IO::Zlib;
\$gtReport->open(\$inputFile, 'rb');

my \$outFileName = \$inputFile =~ s/.csv/.lgen/gr;
open(outputFile,">",\$outFileName) or die \$!;

my \$lineNumber = 0;

while(my \$line = <\$gtReport>) 
{
  next if \$lineNumber++ < \$numberOfHeaderLines+1;
  my (\$column1, \$column2, \$column3, \$column4) = (split /,/, \$line)[1, 0, 2, 3];
  my %replace = ( 
        "-" => "0", 
        "I" => "0", 
        "D" => "0" 
  );
  \$column3 = \$column3 =~ s/(-|I|D)/\$replace{\$1}/gr;
  \$column4 = \$column4 =~ s/(-|I|D)/\$replace{\$1}/gr;
  #\$column3 = \$column3 =~ s/(-)|(I)|(D)/0/gr;
  #\$column4 = \$column4 =~ s/-/0/gr;
  print outputFile "\$column1 \$column1 \$column2 \$column3 \$column4\n";
}
\$gtReport->close;
close(outputFile);
