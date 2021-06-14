#!/usr/bin/perl

# Script to check plink .bim files against HRC/1000G/TOPMed for strand, id names, positions, alleles, ref/alt assignment
# W.Rayner 2015 - 2020
# wrayner@well.ox.ac.uk
#
# Version 4.2
#
#  -v4.0
#  - removes SNPs not in the reference panel
#  - removes all A/T G/C SNPs with MAF >40% in the reference data set
#  - removed all SNPs with an AF difference >0.2, between reference and data set frequency file, frequency file is 
#    expected to be a plink frequency file with the same number of SNPs as the bim file
#
#  -v4.1
#  - removes duplicates that may be introduced with the position update
#  - fixed problem with AF vs MAF not correctly identifying AF diff > 0.2
#  - fixed issue with checks introduced as fix for previous bug that resulted in duplicates being introduced to the final file
#  -v4.2
#  - added ability to use multi ethnic 1000G reference panel
#  - added command line options to specific multiple reference panels/ethnicities
#  - added reference panel used to output filenames
#  - reworked all output filenames to more user friendly format
#  - added verbose option for toggling amount written to log file
#  - Indels are now added to the exclude file
#  -v4.2.1
#  - Fixed bug causing all SNPs with same position but differing name to be labelled duplicates and removed
#  -v4.2.2  
#  - added ability to set allele frequency difference threshold, was fixed at 0.2
#  - added ability to not exclude based on allele frequency difference
#  - changed force allele file to now include all variants regardless of if they need a change to fix a bug where variants could be missed 
#  -v4.2.3
#  - fixed bug whereby a SNP with an incorrect but identical position to the HRC would not be moved chromosome, even though the position was updated
#  -v4.2.4
#  - Added chrX support from r1.1 of HRC
#  -v4.2.5
#  - Updated, adding Chr X to plink command file
#  -v4.2.6
#  - Added gzip reading for reference panels
#  - Added a check of variant #'s between freq and bim file
#  -v4.2.7
#  - Updated code for determining window width to allow better compatibility with Windows
#  -v4.2.8
#  - Added -c flag to allow checking of a subset or individual chromosomes
#  -v4.2.9
#  - Added ability to read .bim files with 1,2,3,4 allele coding
#  - Fixed bug which would give a null entry if no differences were found between the reference and bim file
#  -v4.2.10
#  - Added checking to ensure bim and frq filenames and paths are valid
#  -v4.2.11
#  - Changed the commands to update and retain the Ref/Alt alleles in the plink conversion commands
#  -v4.2.12
#  - Added -a flag to disable the automatic removal of palindromic SNPs with MAF > 0.4
#  -v4.2.13
#  - Added ability to read the bgzipped reference panels, as well as plain gzipped files
#  - Added -l flag to set path to preferred plink executable   
#  - Added better support for the paths to the plink files in the shell script
#  - Added -o flag to allow the final output path for all the plink and VCF files to be specified
#  -v4.3
#  - Added TOPMed and rebuilt code to use less memory
#  - Removed the interactive terminal size check, this version will work both interactively and on a cluster
#
#
#
#
#
# NOTES:
# Script is based on release 1 of the HRC, filename HRC.r1.GRCh37.autosomes.mac5.sites.tab, can be overridden with the -r flag
# # 1000G has X and is also checked
# At the moment indels are not checked, beyond counting them in the bim file, they are exluded from the bim files
# May add -i flag to keep or even check these in future but tricky due to the nature of the labelling (I/D vs actual alleles)
# Script needs ~20Gb RAM to run for HRC, 32GB recommended, and 64GB for 1000G.
# 
#

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
#use Term::ReadKey   qw/ GetTerminalSize /;
use File::Spec;

$| = 1;

my $columns = getwidth(); 
my $mid = int($columns/2+0.5);
my $version = '4.3';

print "\n\n";
printf("%*s", $mid+24, "Script to check plink .bim files against HRC/1000G for\n");
printf("%*s", $mid+25, "strand, id names, positions, alleles, ref/alt assignment\n");
printf("%*s", $mid+5, "William Rayner 2015-2020\n");
printf("%*s", $mid+6, "wrayner\@well.ox.ac.uk\n");
print "\n";
printf("%*s", $mid+5, "Version $version\n\n\n");

# default input filenames (HRC or 1000G file name)
my $hrc_file = 'HRC.r1.GRCh37.autosomes.mac5.sites.tab';
my $kg_file =  '1000GP_Phase3_combined.legend';

my $in_file;
my $bim_file;
my $frq_file;
my $hrcflag = 0;
my $kgflag = 0;
my $verbose;
my %id;
my %refalt;
my %rs;
my $indel = 0;
my $mismatchpos = 0;
my $nomatch = 0;
my $unchanged = 0;
my $strand = 0;
my $nostrand = 0;
my $nomatchalleles = 0;
my $idmismatch = 0;
my $idmatch = 0;
my $nothing = 0;
my $idallelemismatch = 0;
my $hrcdot = 0;
my $total = 0;
my $altchr = 0;
my %AltAf;
my %af;
my $palin = 0;
my $allelediff = 0;
my %seen;
my $duplicate = 0;
my @alleles;
my @tr_alleles;
my $population;
my $referenceused = 'HRC';
my $indelflag = 0;
my $plotflag = 0;
my $threshold = 0.2; #allele frequency difference threshold, set as default 0.2
my $noexclude = 0;
my $chrflag = 0;
my $palinFlag = 0;
my $volume; 
my $directories;
my $parsed_file;
my $abs_path;
my $plink;
my $plink_run_script = 'Run-plink.sh';
my $outpath;
my $abs_outpath;
my %chromosomes;
my %chrPosBim;
my %rsBim;


GetOptions
 (
 "f|frequency=s" => \$frq_file,   # input frequency filename (from plink)
 "b|bim=s"       => \$bim_file,   # input plink bim filename
 "h|hrc"         => \$hrcflag,    # flag to set HRC check
 "r|ref=s"       => \$in_file,    # input reference file (1000G or HRC) 
 "g|1000g"       => \$kgflag,     # flag to set 1000G check
 "p|pop=s"       => \$population, # 1000G population frequency
 "v|verbose"     => \$verbose,    # set verbose logging
 "t|threshold=s" => \$threshold,  # set the allele frequency difference threshold
 "n|noexclude"   => \$noexclude,  # sets flag to keep all SNPs regardless of allele frequency differences  
 "c|chromosome"  => \$chrflag,    # sets flag to say using a smaller than expected reference panel
 "i|indels"      => \$indelflag,  # sets flag for keeping/checking indels in the bim file
 "a|acgt"        => \$palinFlag,  # sets flag to keep the palindromic SNPs
 "l|plink=s"     => \$plink,      # sets the path to the plink executable
 "o|output=s"    => \$outpath,    # sets the output path for the final files
 "x|xyplot"      => \$plotflag    # sets flag for invoking frequency plots at the end of the comparison, requires GD or R
  );

#$plotflag = 0; # temporary setting to disable plotting until the scripts are sorted.

if (!$hrcflag and !$kgflag)
 {
 print "ERROR: you must specify HRC or 1000G reference panel using -g or -h\n";
 usage();
 die "exiting\n";
 }

if (!$plink) 
 {
 $plink = 'plink';
 }
elsif (!-e $plink)
 {
 $plink = 'plink';
 }
else
 {
 }

if ($in_file)
 {
 # reference file has been set via command line so no need to use defaults
 }
elsif ($hrcflag)
 {
 $in_file = $hrc_file;
 }
elsif ($kgflag)
 {
 $in_file = $kg_file;
 }
else
 {
 print "ERROR: no frequency file specified\n";
 usage();
 die "exiting\n";
 }

# set text for filenames and logging
if ($hrcflag)
 {
 $referenceused = 'HRC';
 }
elsif ($kgflag)
 {
 $referenceused = '1000G';
 }

if (!$frq_file)
 {
 print "ERROR: no .frq file specified\n";
 usage();
 die "exiting\n";
 }

if (!$bim_file)
 {
 print "ERROR: no bim file specified\n";
 usage();
 die "exiting\n";
 }

if ($kgflag and !$population)
 {
 print "WARNING: 1000G panel specified and no population selected will default to ALL\n";
 $population = 'ALL';
 }

print "Options Set:\n";
print "Reference Panel:             $referenceused\n";
print "Bim filename:                $bim_file\n";
print "Reference filename:          $in_file\n";
print "Allele frequencies filename: $frq_file\n";
print "Plink executable to use:     $plink\n";
print "\n";

if ($chrflag)
 {
 print "Chromosome flag set:         Yes\n";
 }
else
 {
 print "Chromosome flag set:         No\n";
 }

if (!$noexclude)
 {
 print "Allele frequency threshold:  $threshold\n";
 }
else
 {
 print "Will not exclude any SNPs based on allele frequency differences\n";
 }

# optional flags to be highlighted if set 
if ($kgflag)
 {
 print "Population for 1000G:        $population\n";
 }

if ($verbose)
 {
 print "Verbose logging flag set\n";
 }

if ($palinFlag)
 {
 print "Palindromic Flag set:\n";
 print "   A/T G/C SNPs with MAF > 0.4 will be included in the checks\n";
 }

print "\n";

my $bim_count = 0;
my $frq_count = 0;

if (-e $bim_file)
 {
 read_bim_sites($bim_file);
 $bim_count = get_counts($bim_file);
 ($volume,$directories,$parsed_file) = File::Spec->splitpath($bim_file);
 # print "Volume $volume\nDirectories $directories\nFile $parsed_file\n";
 $abs_path = File::Spec->rel2abs($directories);
 print "Path to plink bim file: $abs_path\n";
 }
else
 {
 print "ERROR: Unable to open specified bim file: $bim_file\n";
 print "Please check path and filename and try again\n";
 exit;
 }
 
if (-e $frq_file)
 {
 $frq_count = get_counts($frq_file);
 }
else
 {
 print "ERROR: Unable to open specified frequency file: $frq_file\n";
 print "Please check path and filename and try again\n";
 exit;
 }

if ($bim_count != ($frq_count-1))
 {
 print "WARNING: The number of variants in the bim and frq files are different\n";
 }

my $allele_coding_flag = check_allele_coding($bim_file);
if ($allele_coding_flag == 2)
 {
 print "ERROR: Alleles coded as 1,2 this coding is not supported\n";
 exit;
 }
 
if ($outpath)
 {
 print "Writing output files to: $outpath\n\n";
 }

print "\nReading $in_file\n";

if ($hrcflag)
 {
 read_hrc($in_file)
 }
elsif ($kgflag)
 {
 read_kg($in_file, $population)
 }
else
 {
 #should never end up here, but just in case
 print "ERROR: you must specify HRC or 1000G type reference panel using -g or -h\n";
 usage();
 die "exiting\n";
 }

open IN, "$bim_file" or die $!; # bim file
open FRQ, "$frq_file" or die $!; # frequency file 

while (<FRQ>)
 {
 chomp;
 s/^\s+//;
 my @temp = split/\s+/;
 $af{$temp[1]} = $temp[4];
 }
close FRQ;

#open all the output files for the different plink update lists
#my $filename = fileparse($bim_file);
#$filename =~ /(.*)\.bim$/;
$parsed_file =~ /(.*)\.bim$/;
my $file_stem = $1;
my $plink_file_path = File::Spec->catfile($abs_path, $file_stem);

if ($outpath)
 {
 my @created = make_path($outpath);
 $abs_outpath = File::Spec->rel2abs($outpath);
 }

my $logfilename = File::Spec->catfile($abs_path, 'LOG-'.$file_stem.'-'.$referenceused.'.txt');
if ($outpath)
 {
 $logfilename = File::Spec->catfile($abs_outpath, 'LOG-'.$file_stem.'-'.$referenceused.'.txt');
 }
print "\n\nDetails written to log file: $logfilename\n";
#print "$logfilename\n";
open L, ">$logfilename" or die $!;

print "\nCreating variant lists\n";

my $forcefile = File::Spec->catfile($abs_path, 'Force-Allele1-'.$file_stem.'-'.$referenceused.'.txt');
if ($outpath)
 {
 $forcefile = File::Spec->catfile($abs_outpath, 'Force-Allele1-'.$file_stem.'-'.$referenceused.'.txt');
 }
print "$forcefile\n";
open F, ">$forcefile" or die $!;

my $strandfile = File::Spec->catfile($abs_path, 'Strand-Flip-'.$file_stem.'-'.$referenceused.'.txt');
if ($outpath)
 {
 $strandfile = File::Spec->catfile($abs_outpath, 'Strand-Flip-'.$file_stem.'-'.$referenceused.'.txt');
 }
print "$strandfile\n";
open S, ">$strandfile" or die $!;

my $idfile = File::Spec->catfile($abs_path, 'ID-'.$file_stem.'-'.$referenceused.'.txt');
if ($outpath)
 {
 $idfile = File::Spec->catfile($abs_outpath, 'ID-'.$file_stem.'-'.$referenceused.'.txt');
 }
print "$idfile\n";
open I, ">$idfile" or die $!;

my $posfile = File::Spec->catfile($abs_path, 'Position-'.$file_stem.'-'.$referenceused.'.txt');
if ($outpath)
 {
 $posfile = File::Spec->catfile($abs_outpath, 'Position-'.$file_stem.'-'.$referenceused.'.txt');
 }
print "$posfile\n";
open P, ">$posfile" or die $!;

my $chrfile = File::Spec->catfile($abs_path, 'Chromosome-'.$file_stem.'-'.$referenceused.'.txt');
if ($outpath)
 {
 $chrfile = File::Spec->catfile($abs_outpath, 'Chromosome-'.$file_stem.'-'.$referenceused.'.txt');
 }
print "$chrfile\n";
open C, ">$chrfile" or die $!;

my $excludefile = File::Spec->catfile($abs_path, 'Exclude-'.$file_stem.'-'.$referenceused.'.txt'); 
if ($outpath)
 {
 $excludefile = File::Spec->catfile($abs_outpath, 'Exclude-'.$file_stem.'-'.$referenceused.'.txt'); 
 }
print "$excludefile\n";
open E, ">$excludefile" or die $!;

my $plotfile = File::Spec->catfile($abs_path, 'FreqPlot-'.$file_stem.'-'.$referenceused.'.txt');
if ($outpath)
 {
 $plotfile = File::Spec->catfile($abs_outpath, 'FreqPlot-'.$file_stem.'-'.$referenceused.'.txt');
 }
print "$plotfile\n";
open PL, ">$plotfile" or dir $!;

print L "Version $version\n\n";
print L "Options Set:\n";
print L "Reference Panel:             $referenceused\n";
print L "Bim filename:                $bim_file\n";
print L "Reference filename:          $in_file\n";
print L "Allele frequencies filename: $frq_file\n";
print L "Plink executable to use:     $plink\n";

if ($chrflag)
 {
 print L "Chromosome flag set:         Yes\n";
 }
else
 {
 print L "Chromosome flag set:         No\n";
 }

if (!$noexclude)
 {
 print L "Allele frequency threshold:  $threshold\n";
 }
else
 {
 print L "Will not exclude any SNPs based on allele frequency differences\n";
 }

if ($kgflag)
 {
 print L "Population for 1000G:        $population\n";
 }
if ($verbose)
 {
 print L "Verbose logging flag set\n";
 } 

print L "\nPath to plink bim file: $abs_path\n";

if ($outpath)
 {
 print L "Writing output files to: $outpath\n\n";
 }

print L "\n\nDetails written to log file: $logfilename\n";
print L "\nCreating variant lists\n";
print L "$forcefile\n";
print L "$strandfile\n";
print L "$idfile\n";
print L "$posfile\n";
print L "$chrfile\n";
print L "$excludefile\n";
print L "$plotfile\n";
print L "\n\n";

while (<IN>)
 {
 chomp;
 my $indelflag = 0;
 $total++;
 
 #split line
 my @temp = split/\s+/;
 if ($temp[0] <= 23) 
  {
  #set chr-position id for checks 
  my $chrpos = $temp[0].'-'.$temp[3];
  
  #chromosomes in the file, for final update and split
  $chromosomes{$temp[0]}++;

  #set alleles for strand and ref/alt checks
  my $allele1 = $temp[4];
  my $allele2 = $temp[5];
  if ($allele_coding_flag == 1)
   {
   $allele1 =~ s/1234/ACGT/;
   $allele2 =~ s/1234/ACGT/;
   }
  $alleles[0] = $allele1;
  $alleles[1] = $allele2;
  my $bim_alleles = $allele1.':'.$allele2;
  my @sorted_alleles = sort {$a cmp $b} @alleles;
  my $sort_alleles = $sorted_alleles[0].':'.$sorted_alleles[1];
  
  $allele1 =~ tr/ACGT/TGCA/;
  $allele2 =~ tr/ACGT/TGCA/;
  
  #create alternate strand alleles for this SNP
  $tr_alleles[0] = $allele1;
  $tr_alleles[1] = $allele2;
  my @sorted_tr_alleles = sort {$a cmp $b} @tr_alleles;
  my $sort_tr_alleles = $sorted_tr_alleles[0].':'.$sorted_tr_alleles[1];
  my $ChrPosTrAlleles =  $chrpos.'-'.$sort_tr_alleles; 
  
  # if indel, adjust position by -1 before checking
  if ($temp[4] eq '-' or $temp[5] eq '-' or $temp[4] eq 'I' or $temp[5] eq 'I' or $temp[4] eq 'D' or $temp[5] eq 'D')
   {
   $temp[3] = $temp[3] - 1;
   $indel++;
   $indelflag = 1;
   print E "$temp[1]\n";
   }
  # no indels in r1 of HRC so skip for now, indels will flag up at the same pos/different alleles stage 
  # due to the way Illumina represent as -/A but in 1000G/HRC represented as T/TA
  elsif ($id{$chrpos}) # position Match
   {
   my $ChrPosAlleles =  $chrpos.'-'.$sort_alleles; # set flag for duplicate removal, based on chr pos alleles
   if ($seen{$ChrPosAlleles}) # chr-position has been seen before remove from the data set, 
    {
    print E "$temp[1]\n";
    if ($verbose)
     {
     print L "Duplicate $temp[1]\t$chrpos\n";
     }
    $duplicate++;
    }
   else
    {
    $seen{$ChrPosAlleles} = 1;
    $seen{$ChrPosTrAlleles} = 1;

    my $idmismatching = 0;
    if ($id{$chrpos} eq $temp[1]) # id match
     {
     $idmatch++;
     }
    else # positions the same but ids are not, is this a genuine difference of name, or a mismapped SNP
     {
     if ($rs{$temp[1]} and $id{$chrpos} ne '.' and $chrpos ne $rs{$temp[1]}) # checks for mismapped SNP
      {
      # SNP rs exists in the reference dataset, and the HRC name is not = '.' (to exclude spurious mismatches where HRC name is not assigned),
      # and position in the reference, based on rs id, is different from the current position
      # check to see if the id is used elsewhere, if it is the wrong SNP has been chosen by position, 
      # happens when there are adjacent SNPs and bim file SNP position is out, correct here by changing
      # bim file location to that of the rs id in HRC
      
      if ($verbose)
       {
       print L "$temp[1]\t$id{$chrpos}\t$chrpos\t$rs{$temp[1]}\t$bim_alleles\t$refalt{$chrpos}\t$refalt{$rs{$temp[1]}}\n";
       }
       
      # set exclusions for this SNP based on the new position to remove duplicates
      $ChrPosTrAlleles =  $rs{$temp[1]}.'-'.$sort_tr_alleles;
      $ChrPosAlleles =  $rs{$temp[1]}.'-'.$sort_alleles;
     
      $chrpos = $rs{$temp[1]};
      
      #check whether this SNP exists or not already in the data set at this new position as could create duplicate in the data set otherwise
      if ($seen{$ChrPosAlleles} or $seen{$ChrPosTrAlleles})
       {
       print E "$temp[1]\n";
       if ($verbose)
        {
        print L "Duplicate $temp[1]\t$rs{$temp[1]}\n";
        }
       $duplicate++;
       }
      else
       {
       $chrpos =~ /(.*)\-(.*)/;
       print P "$temp[1]\t$2\n";
       print C "$temp[1]\t$1\n";
       $mismatchpos++;
       }
      
      # update hash to reflect that this SNP has been seen in the data set before
      $seen{$ChrPosAlleles} = 1;
      $seen{$ChrPosTrAlleles} = 1; 
      }
     else
      {
      $idmismatch++;
      print I "$temp[1]\t$id{$chrpos}\n"; #update ID
      if ($verbose)
       {
       print L "$temp[1]\t$id{$chrpos}\t$chrpos\t$bim_alleles\t$refalt{$chrpos}\n";
       }
      $idmismatching = 1;
      }
     }
     
    my $checking = check_strand($refalt{$chrpos}, $bim_alleles, $temp[1], $AltAf{$chrpos}, $af{$temp[1]});
    if (!$checking)
     {
     if ($idmismatching)
      {
      #alleles and ids don't match
      $idallelemismatch++;
      if ($id{$chrpos} eq '.')
       {
       $hrcdot++;
       }
      }
     
     if ($verbose)
      {
      print L "nomatch $temp[1]\t$id{$chrpos}\t$refalt{$chrpos}\t$bim_alleles\n";
      if ($refalt{$chrpos} eq 'N:N')
       {
       print L "$temp[1] MultiAllelic\n";
       }
      }
     $nomatchalleles++;   
    
     #print to an exclusion file
     print E "$temp[1]\n";
     }
    }
   }
  elsif ($rs{$temp[1]}) #match on id, check why position did not match, set position to reference
   {
   my $ChrPosAlleles = $rs{$temp[1]}.'-'.$sort_alleles; # set flag for duplicate removal, based on chr pos alleles
   if ($seen{$ChrPosAlleles}) # chr-position has been seen before remove from the data set, 
    {
    print E "$temp[1]\n";
    if ($verbose)
     {
     print L "Duplicate $temp[1]\t$rs{$temp[1]}\n";
     }
    $duplicate++;
    }
   else
    {
    my @ChrPosRef = split(/-/, $rs{$temp[1]});
    $seen{$ChrPosAlleles} = 1;
    $seen{$ChrPosTrAlleles} = 1;

    print C "$temp[1]\t$ChrPosRef[0]\n"; #print element [0] chromosome
    print P "$temp[1]\t$ChrPosRef[1]\n"; #print element [1] position
     
    $mismatchpos++;
    my $checking = check_strand($refalt{$rs{$temp[1]}}, $bim_alleles, $temp[1], $AltAf{$rs{$temp[1]}}, $af{$temp[1]});
    
    if (!$checking)
     {
     if ($verbose)
      {
      print L "nomatch $temp[1]\t$id{$rs{$temp[1]}}\t$refalt{$rs{$temp[1]}}\t$bim_alleles\n";
      }
     $nomatchalleles++;   
    
     #print to an exclusion file
     print E "$temp[1]\n";
     }
    } 
   }
  else # no match on position or variant id, check +/- 1???
   {
   $nothing++;
   if ($verbose)
    {
    print L "Not in $referenceused\t$temp[1]\n";
    }
   print E "$temp[1]\n";
   }
  } #end of Chromosome check section
 else
  {# total all the skipped lines here
  $altchr++;
  print E "$temp[1]\n";
  }
 }
 
#print "Total bim File Rows $total\n";

my $check_total = $idmatch + $idmismatch + $mismatchpos + $nothing + $altchr ;
my $check_total1 = $unchanged + $nomatch;
my $pos_check = $idmatch + $idmismatch;
my $worked_check = $idmatch + $idmismatch + $mismatchpos;
my $worked_check1 = $strand + $nostrand;

print "\n\nMatching to $referenceused\n";
print "\nPosition Matches\n ID matches $referenceused $idmatch\n ID Doesn't match $referenceused $idmismatch\n Total Position Matches $pos_check\nID Match\n Position different from $referenceused $mismatchpos\nNo Match to $referenceused $nothing\nSkipped (MT) $altchr\nTotal in bim file $total\nTotal processed $check_total\n\n"; 
print "Indels $indel\n\n";
print "SNPs not changed $unchanged\nSNPs to change ref alt $nomatch\nStrand ok $strand\nTotal Strand ok $check_total1\n\n";
print "Strand to change $nostrand\nTotal checked $worked_check\nTotal checked Strand $worked_check1\n";
if (!$noexclude)
 {
 print "Total removed for allele Frequency diff > $threshold $allelediff\n";
 }
print "Palindromic SNPs with Freq > 0.4 $palin\n\n";
print "\nNon Matching alleles $nomatchalleles\n";
print "ID and allele mismatching $idallelemismatch; where $referenceused is . $hrcdot\n";
print "Duplicates removed $duplicate\n";

#print L "Total bim File Rows $total\n";
print L "\nMatching to $referenceused\n";
print L "\nPosition Matches\n ID matches $referenceused $idmatch\n ID Doesn't match $referenceused $idmismatch\n Total Position Matches $pos_check\nID Match\n Position different from $referenceused $mismatchpos\nNo Match to $referenceused $nothing\nSkipped (MT) $altchr\nTotal in bim file $total\nTotal processed $check_total\n\n"; 
print L "Indels $indel\n\n";
print L "SNPs not changed $unchanged\nSNPs to change ref alt $nomatch\nStrand ok $strand\nTotal Strand ok $check_total1\n\n";
print L "Strand to change $nostrand\nTotal checked $worked_check\nTotal checked Strand $worked_check1\n";
if (!$noexclude)
 {
 print L "Total removed for allele Frequency diff > $threshold $allelediff\n";
 }
print L "Palindromic SNPs with Freq > 0.4 $palin\n\n";
print L "\nNon Matching alleles $nomatchalleles\n";
print L "ID and allele mismatching $idallelemismatch; where $referenceused is . $hrcdot\n";
print L "Duplicates removed $duplicate\n";

print "\n\nWriting plink commands to: $plink_run_script\n";
write_shell_script();

#close the file with the allele frequencies
close PL;
close L;
close I;
close P;
close C;
close E;
close S;

#clear the HRC hashes
#%id = ();
#%rs = ();
#%AltAf = ();
#%refalt = ();

#my $checkinc = inc_check();
#my $checkr = r_check();

#if ($plotflag )
# {
# if (-e 'plotting.pl' and $checkinc)
#  { # check includes can be satisfied
#  my $arg = "perl plotting.pl $plotfile";
#  exec $arg;
#  }
# elsif (-e 'run-r.pl' and $checkr)
#  {
#  # if R available then run that insted of Perl script
#  my $arg = "run-r.pl $plotfile";
#  exec $arg;
#  }
# else
#  {
#  print "Plotting selected but neither R or GD::Graph appear installed on your system\n";
#  print "Plotting script (plotting.pl) not found in directory with this script\n";
#  print "manually use \'perl plotting.pl $plotfile\' to invoke\n";
#  print "To get the Perl plotting to run install GD by running the following commands\n";
#  print "sudo apt-get install cpanminus\nsudo apt-get -y install libgd2-xpm-dev build-essential\nsudo cpanm GD::Graph\n";
#  }
# }
#
#sub r_check
# {
# my $check = 0;
# my $return = `whereis R`;
# #print "$return\n";
# if ($return =~ /^R\:\n/)
#  {
#  print "R executable not found\n";
#  }
# else
#  {
#  $check = 1;
#  }
# return $check;
# }
#
#sub inc_check
# {
# my $check = 0;
# eval {use GD::Graph::points; };
# if (!$@)
#  {
#  $check = 1; 
#  }
# return $check; 
# }

sub write_shell_script
 {
 # shell script for running plink
 my $run_script = File::Spec->catfile($abs_path, $plink_run_script);
 if ($outpath)
  {
  $run_script = File::Spec->catfile($abs_outpath, $plink_run_script);
  }
 open SH, ">$run_script" or die $!;

 my $tempcount = 1;
 my $tempfile = File::Spec->catfile($abs_path,'TEMP'.$tempcount);
 #remove SNPs
 print SH "$plink --bfile $plink_file_path --exclude $excludefile --make-bed --out $tempfile\n"; 

 #change chromosome
 print SH "$plink --bfile $tempfile --update-map $chrfile --update-chr --make-bed --out ";
 $tempcount++;
 $tempfile = File::Spec->catfile($abs_path, 'TEMP'.$tempcount);
 print SH "$tempfile\n"; 

 #change positions
 print SH "$plink --bfile $tempfile --update-map $posfile --make-bed --out ";
 $tempcount++;
 $tempfile = File::Spec->catfile($abs_path, 'TEMP'.$tempcount);
 print SH "$tempfile\n";

 #flip strand
 print SH "$plink --bfile $tempfile --flip $strandfile --make-bed --out ";
 $tempcount++;
 $tempfile = File::Spec->catfile($abs_path, 'TEMP'.$tempcount);
 print SH "$tempfile\n";

 #update ids
 #remove the following 4 lines if you want don't want to update the SNP identifiers to match the HRC
 #print SH "$plink --bfile $tempfile --update-map $idfile --update-name --make-bed --out ";
 #$tempcount++;
 #$tempfile = 'TEMP'.$tempcount;
 #print SH "$tempfile\n";

 #force alleles
 my $newfile = File::Spec->catfile($abs_path, $file_stem.'-updated');
 if ($outpath)
  {
  $newfile = File::Spec->catfile($abs_outpath, $file_stem.'-updated');
  } 
 print SH "$plink --bfile $tempfile --a2-allele $forcefile --make-bed --out $newfile\n";

 #split into per chromosome files
 for (my $i = 1; $i <= 23; $i++)
  {
  if ($chromosomes{$i})
   {
   my $perchrfile = $newfile.'-chr'.$i;
   #my $perchrfile = File::Spec->catfile($abs_path, $newfile.'-chr'.$i);
   #if ($outpath)
   # {
   # $perchrfile = File::Spec->catfile($abs_outpath, $newfile.'-chr'.$i);
   # }
   print SH "$plink --bfile $newfile --real-ref-alleles --make-bed --chr $i --out $perchrfile\n";
   print SH "$plink --bfile $newfile --real-ref-alleles --recode vcf --chr $i --out $perchrfile\n";
   }
  }
 print SH "rm TEMP*\n";
 close SH;
 }


sub check_strand
 {
 my $check = 0;
 my $a1 = $_[0]; #HRC alleles
 my $a2 = $_[1]; #bimfile alleles
 my $id = $_[2];
 my $altaf = $_[3]; #HRC alt allele frequency
 my $bimaf = $_[4]; #bim file allele frequency
 my $diff = 0;
 my $maf = 0;
 
 my @alleles1 = split(/\:/, $a1);
 my @alleles2 = split(/\:/, $a2);
 #set the ref allele
 my $ref = $alleles1[0];
 
 if ($verbose)
  {
  print L "$id\t$a1\t$a2";
  }
  
 # flip one set and check if they match opposite strand
 $a2 =~ tr/ACGTN/TGCAN/;
 
 if ($verbose)
  {
  print L "\t$a2\n";
  }
  
 my @allelesflip = split(/\:/, $a2); 
 
 if ($altaf > 0.5)
  {
  $maf = 1 - $altaf;
  }
 else
  {
  $maf = $altaf
  }
  
 if (!$palinFlag)
  {
  #check MAFs for palindromic SNPs first, this is an absolute failure, so return here if conditions not met
  if ($maf > 0.4 and ($a1 eq 'A:T' or $a1 eq 'T:A' or $a1 eq 'G:C' or $a1 eq 'C:G')) 
   {
   print E "$id\n";
   
   if ($verbose)
    {
    print L "$id\t$maf\t$a1\n";
    }
    
   $check = 5;
   $palin++;
   return $check;
   }
  }
 #print PL "$id\t$refaf\t$af\t";
 #$diff = $refaf - $bimaf;

 
 # check alleles/strand are the same and print frequencies to file
 if ($alleles1[0] eq $alleles2[0] and $alleles1[1] eq $alleles2[1])
  { # strand ok, ref/alt ok
  $check = 1;
  $strand++;
  $unchanged++;
  #print PL "1 $alleles1[0]\t$alleles2[0]\t$alleles1[1]\t$alleles2[1]\t";
  my $RefAf = 1 - $altaf;
  $diff = $RefAf - $bimaf;
  print PL "$id\t$RefAf\t$bimaf\t$diff\t1\n";
  }
 elsif ($alleles1[0] eq $alleles2[1] and $alleles1[1] eq $alleles2[0])
  { # strand ok, ref alt swapped
  $check = 2;
  $strand++;
  $nomatch++;
  #print F "$id\t$ref\n";
  #print PL "2 $alleles1[0]\t$alleles2[0]\t$alleles1[1]\t$alleles2[1]\t";
  #my $newaf = 1-$refaf;
  $diff = $altaf - $bimaf;
  print PL "$id\t$altaf\t$bimaf\t$diff\t2\n";
  }
 elsif ($alleles1[0] eq $allelesflip[0] and $alleles1[1] eq $allelesflip[1])
  { # strand flipped, ref alt ok
  $check = 3;
  $nostrand++;
  print S "$id\n";
  #print PL "3 $alleles1[0]\t$allelesflip[0]\t$alleles1[1]\t$allelesflip[1]\t";
  my $RefAf = 1 - $altaf;
  $diff = $RefAf - $bimaf;
  print PL "$id\t$RefAf\t$bimaf\t$diff\t3\n";
  }
 elsif ($alleles1[0] eq $allelesflip[1] and $alleles1[1] eq $allelesflip[0])
  { # strand flipped, ref alt swapped
  $check = 4;
  $nostrand++;
  $nomatch++;
  print S "$id\n";
  #print F "$id\t$ref\n"; 
  #print PL "4 $alleles1[0]\t$allelesflip[0]\t$alleles1[1]\t$allelesflip[1]\t";
  #my $af = 1-$altaf;
  $diff = $altaf - $bimaf;
  print PL "$id\t$altaf\t$bimaf\t$diff\t4\n";
  } 
 else
  {
  $check = 0;
  #print PL "7\n";
  }

 print F "$id\t$ref\n"; #print all variants to Force Allele file so as to ensure none are missed

 if ($diff < 0)
  {
  $diff = $diff * -1;
  }
 
  if ($diff > $threshold)# removed the following for A/T G/C: and ($a1 eq 'A:T' or $a1 eq 'T:A' or $a1 eq 'G:C' or $a1 eq 'C:G'))
   {
   if (!$noexclude) #only add to exclusion if the noexclude flag is not set
    {
    print E "$id\n"; 
    }
    
   if ($verbose)
    {
    print L "$id\t$bimaf\t$altaf\t$diff\n";
    }
    
   $check = 6;
   $allelediff++;
   }
  
 return $check;
 } 
 
sub getwidth
 {
 #my $output = `stty size`;
 #my @rowcols = split(/\s/, $output);
 #my $cols = $rowcols[1];
 
 #my @winsize = &GetTerminalSize(\*STDOUT);
 #my ($cols, $rows, $xpix, $ypix) = @winsize;
 my $cols = 80;
 if (!$cols)
  {
  $cols = 80;
  }
 return $cols;
 } 

sub usage
 {
 print "\nUsage:\nFor HRC:\nperl HRC-1000G-check-bim.pl -b <bim file> -f <Frequency file> -r <Reference panel> -h [-v -t <allele frequency threshold -n]\n";
 print "\nFor 1000G:\nperl HRC-1000G-check-bim.pl -b <bim file> -f <Frequency file> -r <Reference panel> -g -p <population> [-v -t <allele frequency threshold -n]\n";
 print "\n\n";
 printusage("-b --bim", "bim file", "Plink format .bim file");
 printusage("-f --frequency", "Frequency file", "Plink format .frq allele frequency file, from plink --freq command");
 printusage("-r --ref", "Reference panel", "Reference Panel file, either 1000G or HRC");
 printusage("-h --hrc", "", "Flag to indicate Reference panel file given is HRC");
 printusage("-g --1000g", "", "Flag to indicate Reference panel file given is 1000G");
 printusage("-p --pop", "Population", "Population to check frequency against, 1000G only. Default ALL, options ALL, EUR, AFR, AMR, SAS, EAS");
 printusage("-v --verbose", "", "Optional flag to increase verbosity in the log file");
 printusage("-t --threshold", "Freq threshold", "Frequency difference to use when checking allele frequency of data set versus reference; default: 0.2; range: 0-1");
 printusage("-c --chromosome", "Chromsome flag", "Optional flag to indicate to the program you are checking a subset of chromosomes and so to expect a smaller than normal reference panel");
 printusage("-l --plink", "Path to the plink executable to add to the shell script", "Optional flag to indicate which plink executable you want to use in the Run-plink.sh shell script");
 printusage("-o --output", "Path for the output files", "Optional flag to indicate the directory to use for the output files");
 #printusage("-i --indels", "", "Optional flag to keep or exclude indels in the output, default: exclude");
 #printusage("-x --xyplot", "", "Optional flag to invoke plotting of the data set vs reference allele frequencies");
 printusage("-n --noexclude", "", "Optional flag to include all SNPs regardless of allele frequency differences, default is exclude based on -t threshold, overrides -t");
 print "\n\n";
 }
 
sub printusage
 {
 my $option = $_[0];
 my $file = $_[1];
 my $text = $_[2];
 
 my $cols = getwidth(); 

 printf("%-18s", $option);
 printf("%-15s", $file);
 my $textwidth = $cols - 36; #33 plus 3 spaces
 my $newtext = substr ($text, 0, $textwidth);
 
 if (length($text) > $textwidth)
  {
  printf("  %*s", $textwidth, $newtext);
  for (my $i = $textwidth; $i < length($text); $i+= $textwidth)
   {
   my $newtext = substr ($text, $i, $textwidth);
   $newtext =~ s/^\s//;
   printf("%33s", " ");
   printf("   %*s", -$textwidth, $newtext);
   }
  }
 else
  {
  printf("  %*s", -$textwidth, $newtext);
  }
 print "\n"; 
 } 
 
sub read_hrc
 {
 my $file = $_[0];
 #check if file has .gz suffix
 my $zipped = checkgz($file);
 my $z;
 
 if ($zipped)
  {
  open $z, '-|', '/bin/gunzip', '-c', "$file";
  if (!$z)
   {
   $z = new IO::Uncompress::Gunzip "$file" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
   }
  }
 else
  {
  open $z, "$file" or die $!;
  }
  
 while (<$z>)
  {
  chomp;
  if (!/\#.*/)
   {
   if ($. % 2000000 == 0)
    {
    print "$.\n";
    }
   my @temp = split/\s+/;
   if ($temp[0] eq 'X')
    {
    $temp[0] = 23;
    }
   my $chrpos = $temp[0].'-'.$temp[1];
   if ($chrPosBim{$chrpos} or $rsBim{$temp[2]}) #only read in variant if this position or rs number exists in the supplied bim file
    {
    $id{$chrpos} = $temp[2];
    $rs{$temp[2]} = $chrpos;
    $refalt{$chrpos} = $temp[3].':'.$temp[4];
    $AltAf{$chrpos} = $temp[7];
    }
   }
  }
 if ($. < 10000000 and !$chrflag)
  {
  print "\nERROR: Reference appears smaller than expected this can be due to bgzip\nor by using a subset of the genome without the -c flag\nBgzip is not supported with the current Perl Library, please use the unzipped version\n";
  exit;
  }
 print "$.\n ...Done\n";
 close IN;
 }
 
sub read_kg
 {
 my $file = $_[0];
 my $pop = $_[1];
 my $freqcol;
 my $typecol = 0;
 
 my $zipped = checkgz($file);
 my $z;
  
 if ($zipped)
  {
  print "Reference Panel is zipped\n";
  open $z, '-|', '/bin/gunzip', '-c', "$file";
  if (!$z) # if cannot find gunzip then try uising the library in case not bgzipped
   {
   $z = new IO::Uncompress::Gunzip "$file" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
   }
  }
 else
  {
  open $z, "$file" or die $!;
  }
 
 my $header = <$z>;
 chomp $header;
 my @titles = split(/\s+/, $header);
 
 for (my $i = 0; $i <= $#titles; $i++)
  {
  if ($titles[$i] eq $pop)
   {
   $freqcol = $i;
   }
  if ($titles[$i] eq 'TYPE')
   {
   $typecol = $i;
   }
  }
 
 if (!$freqcol)
  {
  print "ERROR: Population specified, $pop not found in $file\n";
  die;
  }
 
 while (<$z>)
  {
  chomp;
  if ($. % 2000000 == 0)
   {
   print "$.\n";
   }
  my @temp = split/\s+/;
  my $rsID = 'rs';

  my $chrpos = $temp[1].'-'.$temp[2];
  if ($temp[0] =~ /^rs.*/)
   {
   my @tempids = split(/\:/, $temp[0]); #take first element = rs number
   $rsID = $tempids[0];
   $rs{$rsID} = $chrpos;
   }
  if ($chrPosBim{$chrpos} or $rsBim{$rsID}) #only read in variant if this position or rs number exists in the supplied bim file
   {
   $id{$chrpos} = $temp[0];
   $refalt{$chrpos} = $temp[3].':'.$temp[4];
   $AltAf{$chrpos} = $temp[$freqcol];
  
   if ($typecol) # if column with SNP type exists, check for Multiallelic
    {
    if ($temp[$typecol] =~ /^Multiallelic.*/)
     { # set multiallelic alleles to N so will always fail allele check
     $refalt{$chrpos} = 'N:N';
     }
    }
   }
  }
 print "$.\n ...Done\n";
 
 close IN;

 }
 
 
sub checkgz
 {
 my $file = $_[0];
 
 my @filecomponents = split(/\./, $file);
 my $zipped = 0;
 
 if ($filecomponents[$#filecomponents] eq 'gz')
  {
  print "Reference Panel is zipped\n";
  $zipped = 1;
  }
 elsif ($filecomponents[$#filecomponents] eq 'zip')
  {
  print "WARNING: .zip format is not supported, skipping $file\n";
  $zipped = -1;
  } 
 else
  {
  $zipped = 0;
  }
 return $zipped; 
 } 
 
sub get_counts
 {
 my $file = $_[0];
 my $counts = 0;
 open I, "$file" or die $!;
 while (<I>)
  {
  $counts++;
  }
 close I; 
 return $counts;
 }
 
 
sub check_allele_coding
 {
 my $file = $_[0];
 my $total = 0;
 my $code = 0;
 my $oneTwo = 0;
 my $oneTwoThreeFour = 0;
 
 open I, "$file" or die $!;
 while (<I>)
  {
  chomp;
  my @temp = split/\s+/;  
  my $allele1 = $temp[4];
  my $allele2 = $temp[5];
  my $alleles = $allele1.$allele2;
  
  if ($allele1 ne '0' and $allele2 ne '0')
   {
   if ($alleles eq '12' or $alleles eq '21')
    {
    $oneTwo++;
    }
   if ($alleles eq '23' or $alleles eq '34')
    {
    $oneTwoThreeFour++;
    }
   $total++;
   }
  }
 
 if ($oneTwo == $total)
  {
  $code = 2;
  }
 elsif ($oneTwoThreeFour) 
  {
  $code = 1;
  }
  
 return $code; 
 }

sub read_bim_sites
 {
 my $file = $_[0];
 open I, "$file" or die $!;
 while (<I>)
  {
  chomp;
  my @temp = split/\s+/;  
  my $chrPos = $temp[0].'-'.$temp[3];
  $chrPosBim{$chrPos} = 1;
  $rsBim{$temp[1]} = $chrPos; # add an allele check
  }
 close I; 
 }
