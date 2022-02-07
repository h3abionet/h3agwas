#!/usr/bin/env perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use File::Spec;
use Cwd;

our $REVISION = '$Revision: 6b5656c3a8e640eea2984b2037ef4226001dcc94 $';
our $DATE =	'$Date: 2020-06-07 23:56:37 -0400 (Sun,  7 Jun 2020) $';  
our $AUTHOR =	'$Author: Kai Wang <kaichop@gmail.com> $';

our ($verbose, $help, $man);
our ($queryfile, $dbloc);
our ($outfile, $separate, $batchsize, $dbtype, $neargene, $genomebinsize, $geneanno, $regionanno, $filter, $downdb, $buildver, $score_threshold, $normscore_threshold, $minqueryfrac, $expandbin, $splicing_threshold,
	$maf_threshold, $chromosome, $zerostart, $rawscore, $memfree, $memtotal, $sift_threshold, $gff3dbfile, $genericdbfile, $vcfdbfile, $time, $wget, $precedence,
	$webfrom, $colsWanted, $comment, $scorecolumn, $poscolumn, $transfun, $exonsort, $avcolumn, $bedfile, $hgvs, $reverse, $indexfilter_threshold, $otherinfo, $seq_padding, $indel_splicing_threshold, $infoasscore,
	$firstcodondel, $aamatrixfile, $gff3attr, $exonicsplicing, $infosep, $dbm, $idasscore, $thread, $maxgenethread, $mingenelinecount);

our (%valichr, $dbtype1);		# valid chromosome name to process, internal alias of dbtype
our (@precedence, @colsWanted, @avcolumn);	# precedence of functional categories to print, desired columns to be printed out, redefine standard ANNOVAR 5-column input
our $aamatrix;				# amino acid substitution matrix
my ($pad_fh, $cDNA_pad); 		# Write seq pad here
sub printerr;				# declare a subroutine

# codon table
our %codon1 = (TTT=>"F", TTC=>"F", TCT=>"S", TCC=>"S", TAT=>"Y", TAC=>"Y", TGT=>"C", TGC=>"C", TTA=>"L", TCA=>"S", TAA=>"*", TGA=>"*", TTG=>"L", TCG=>"S", TAG=>"*", TGG=>"W", CTT=>"L", CTC=>"L", CCT=>"P", CCC=>"P", CAT=>"H", CAC=>"H", CGT=>"R", CGC=>"R", CTA=>"L", CTG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", ATT=>"I", ATC=>"I", ACT=>"T", ACC=>"T", AAT=>"N", AAC=>"N", AGT=>"S", AGC=>"S", ATA=>"I", ACA=>"T", AAA=>"K", AGA=>"R", ATG=>"M", ACG=>"T", AAG=>"K", AGG=>"R", GTT=>"V", GTC=>"V", GCT=>"A", GCC=>"A", GAT=>"D", GAC=>"D", GGT=>"G", GGC=>"G", GTA=>"V", GTG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");
our %codon3 = (TTT=>"Phe", TTC=>"Phe", TCT=>"Ser", TCC=>"Ser", TAT=>"Tyr", TAC=>"Tyr", TGT=>"Cys", TGC=>"Cys", TTA=>"Leu", TCA=>"Ser", TAA=>"*", TGA=>"*", TTG=>"Leu", TCG=>"Ser", TAG=>"*", TGG=>"Trp", CTT=>"Leu", CTC=>"Leu", CCT=>"Pro", CCC=>"Pro", CAT=>"His", CAC=>"His", CGT=>"Arg", CGC=>"Arg", CTA=>"Leu", CTG=>"Leu", CCA=>"Pro", CCG=>"Pro", CAA=>"Gln", CAG=>"Gln", CGA=>"Arg", CGG=>"Arg", ATT=>"Ile", ATC=>"Ile", ACT=>"Thr", ACC=>"Thr", AAT=>"Asn", AAC=>"Asn", AGT=>"Ser", AGC=>"Ser", ATA=>"Ile", ACA=>"Thr", AAA=>"Lys", AGA=>"Arg", ATG=>"Met", ACG=>"Thr", AAG=>"Lys", AGG=>"Arg", GTT=>"Val", GTC=>"Val", GCT=>"Ala", GCC=>"Ala", GAT=>"Asp", GAC=>"Asp", GGT=>"Gly", GGC=>"Gly", GTA=>"Val", GTG=>"Val", GCA=>"Ala", GCG=>"Ala", GAA=>"Glu", GAG=>"Glu", GGA=>"Gly", GGG=>"Gly");
our %codonfull = (TTT=>"Phenylalanine", TTC=>"Phenylalanine", TCT=>"Serine", TCC=>"Serine", TAT=>"Tyrosine", TAC=>"Tyrosine", TGT=>"Cysteine", TGC=>"Cysteine", TTA=>"Leucine", TCA=>"Serine", TAA=>"Stop", TGA=>"Stop", TTG=>"Leucine", TCG=>"Serine", TAG=>"Stop", TGG=>"Tryptophan", CTT=>"Leucine", CTC=>"Leucine", CCT=>"Proline", CCC=>"Proline", CAT=>"Histidine", CAC=>"Histidine", CGT=>"Arginine", CGC=>"Arginine", CTA=>"Leucine", CTG=>"Leucine", CCA=>"Proline", CCG=>"Proline", CAA=>"Glutamine", CAG=>"Glutamine", CGA=>"Arginine", CGG=>"Arginine", ATT=>"Isoleucine", ATC=>"Isoleucine", ACT=>"Threonine", ACC=>"Threonine", AAT=>"Asparagine", AAC=>"Asparagine", AGT=>"Serine", AGC=>"Serine", ATA=>"Isoleucine", ACA=>"Threonine", AAA=>"Lysine", AGA=>"Arginine", ATG=>"Methionine", ACG=>"Threonine", AAG=>"Lysine", AGG=>"Arginine", GTT=>"Valine", GTC=>"Valine", GCT=>"Alanine", GCC=>"Alanine", GAT=>"Aspartic acid", GAC=>"Aspartic acid", GGT=>"Glycine", GGC=>"Glycine", GTA=>"Valine", GTG=>"Valine", GCA=>"Alanine", GCG=>"Alanine", GAA=>"Glutamic acid", GAG=>"Glutamic acid", GGA=>"Glycine", GGG=>"Glycine");
our %codonr1 = (UUU=>"F", UUC=>"F", UCU=>"S", UCC=>"S", UAU=>"Y", UAC=>"Y", UGU=>"C", UGC=>"C", UUA=>"L", UCA=>"S", UAA=>"*", UGA=>"*", UUG=>"L", UCG=>"S", UAG=>"*", UGG=>"W", CUU=>"L", CUC=>"L", CCU=>"P", CCC=>"P", CAU=>"H", CAC=>"H", CGU=>"R", CGC=>"R", CUA=>"L", CUG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", AUU=>"I", AUC=>"I", ACU=>"T", ACC=>"T", AAU=>"N", AAC=>"N", AGU=>"S", AGC=>"S", AUA=>"I", ACA=>"T", AAA=>"K", AGA=>"R", AUG=>"M", ACG=>"T", AAG=>"K", AGG=>"R", GUU=>"V", GUC=>"V", GCU=>"A", GCC=>"A", GAU=>"D", GAC=>"D", GGU=>"G", GGC=>"G", GUA=>"V", GUG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");
our %codonr3 = (UUU=>"Phe", UUC=>"Phe", UCU=>"Ser", UCC=>"Ser", UAU=>"Tyr", UAC=>"Tyr", UGU=>"Cys", UGC=>"Cys", UUA=>"Leu", UCA=>"Ser", UAA=>"*", UGA=>"*", UUG=>"Leu", UCG=>"Ser", UAG=>"*", UGG=>"Trp", CUU=>"Leu", CUC=>"Leu", CCU=>"Pro", CCC=>"Pro", CAU=>"His", CAC=>"His", CGU=>"Arg", CGC=>"Arg", CUA=>"Leu", CUG=>"Leu", CCA=>"Pro", CCG=>"Pro", CAA=>"Gln", CAG=>"Gln", CGA=>"Arg", CGG=>"Arg", AUU=>"Ile", AUC=>"Ile", ACU=>"Thr", ACC=>"Thr", AAU=>"Asn", AAC=>"Asn", AGU=>"Ser", AGC=>"Ser", AUA=>"Ile", ACA=>"Thr", AAA=>"Lys", AGA=>"Arg", AUG=>"Met", ACG=>"Thr", AAG=>"Lys", AGG=>"Arg", GUU=>"Val", GUC=>"Val", GCU=>"Ala", GCC=>"Ala", GAU=>"Asp", GAC=>"Asp", GGU=>"Gly", GGC=>"Gly", GUA=>"Val", GUG=>"Val", GCA=>"Ala", GCG=>"Ala", GAA=>"Glu", GAG=>"Glu", GGA=>"Gly", GGG=>"Gly");
our %codonrfull = (UUU=>"Phenylalanine", UUC=>"Phenylalanine", UCU=>"Serine", UCC=>"Serine", UAU=>"Tyrosine", UAC=>"Tyrosine", UGU=>"Cysteine", UGC=>"Cysteine", UUA=>"Leucine", UCA=>"Serine", UAA=>"Stop", UGA=>"Stop", UUG=>"Leucine", UCG=>"Serine", UAG=>"Stop", UGG=>"Tryptophan", CUU=>"Leucine", CUC=>"Leucine", CCU=>"Proline", CCC=>"Proline", CAU=>"Histidine", CAC=>"Histidine", CGU=>"Arginine", CGC=>"Arginine", CUA=>"Leucine", CUG=>"Leucine", CCA=>"Proline", CCG=>"Proline", CAA=>"Glutamine", CAG=>"Glutamine", CGA=>"Arginine", CGG=>"Arginine", AUU=>"Isoleucine", AUC=>"Isoleucine", ACU=>"Threonine", ACC=>"Threonine", AAU=>"Asparagine", AAC=>"Asparagine", AGU=>"Serine", AGC=>"Serine", AUA=>"Isoleucine", ACA=>"Threonine", AAA=>"Lysine", AGA=>"Arginine", AUG=>"Methionine", ACG=>"Threonine", AAG=>"Lysine", AGG=>"Arginine", GUU=>"Valine", GUC=>"Valine", GCU=>"Alanine", GCC=>"Alanine", GAU=>"Aspartic acid", GAC=>"Aspartic acid", GGU=>"Glycine", GGC=>"Glycine", GUA=>"Valine", GUG=>"Valine", GCA=>"Alanine", GCG=>"Alanine", GAA=>"Glutamic acid", GAG=>"Glutamic acid", GGA=>"Glycine", GGG=>"Glycine");

# mitochondria codon
our %codon1m = (TTT=>"F", TTC=>"F", TCT=>"S", TCC=>"S", TAT=>"Y", TAC=>"Y", TGT=>"C", TGC=>"C", TTA=>"L", TCA=>"S", TAA=>"*", TGA=>"W", TTG=>"L", TCG=>"S", TAG=>"*", TGG=>"W", CTT=>"L", CTC=>"L", CCT=>"P", CCC=>"P", CAT=>"H", CAC=>"H", CGT=>"R", CGC=>"R", CTA=>"L", CTG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", ATT=>"I", ATC=>"I", ACT=>"T", ACC=>"T", AAT=>"N", AAC=>"N", AGT=>"S", AGC=>"S", ATA=>"M", ACA=>"T", AAA=>"K", AGA=>"*", ATG=>"M", ACG=>"T", AAG=>"K", AGG=>"*", GTT=>"V", GTC=>"V", GCT=>"A", GCC=>"A", GAT=>"D", GAC=>"D", GGT=>"G", GGC=>"G", GTA=>"V", GTG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");
our %codon3m = (TTT=>"Phe", TTC=>"Phe", TCT=>"Ser", TCC=>"Ser", TAT=>"Tyr", TAC=>"Tyr", TGT=>"Cys", TGC=>"Cys", TTA=>"Leu", TCA=>"Ser", TAA=>"*", TGA=>"Trp", TTG=>"Leu", TCG=>"Ser", TAG=>"*", TGG=>"Trp", CTT=>"Leu", CTC=>"Leu", CCT=>"Pro", CCC=>"Pro", CAT=>"His", CAC=>"His", CGT=>"Arg", CGC=>"Arg", CTA=>"Leu", CTG=>"Leu", CCA=>"Pro", CCG=>"Pro", CAA=>"Gln", CAG=>"Gln", CGA=>"Arg", CGG=>"Arg", ATT=>"Ile", ATC=>"Ile", ACT=>"Thr", ACC=>"Thr", AAT=>"Asn", AAC=>"Asn", AGT=>"Ser", AGC=>"Ser", ATA=>"Met", ACA=>"Thr", AAA=>"Lys", AGA=>"*", ATG=>"Met", ACG=>"Thr", AAG=>"Lys", AGG=>"*", GTT=>"Val", GTC=>"Val", GCT=>"Ala", GCC=>"Ala", GAT=>"Asp", GAC=>"Asp", GGT=>"Gly", GGC=>"Gly", GTA=>"Val", GTG=>"Val", GCA=>"Ala", GCG=>"Ala", GAA=>"Glu", GAG=>"Glu", GGA=>"Gly", GGG=>"Gly");
our %codonfullm = (TTT=>"Phenylalanine", TTC=>"Phenylalanine", TCT=>"Serine", TCC=>"Serine", TAT=>"Tyrosine", TAC=>"Tyrosine", TGT=>"Cysteine", TGC=>"Cysteine", TTA=>"Leucine", TCA=>"Serine", TAA=>"Stop", TGA=>"Tryptophan", TTG=>"Leucine", TCG=>"Serine", TAG=>"Stop", TGG=>"Tryptophan", CTT=>"Leucine", CTC=>"Leucine", CCT=>"Proline", CCC=>"Proline", CAT=>"Histidine", CAC=>"Histidine", CGT=>"Arginine", CGC=>"Arginine", CTA=>"Leucine", CTG=>"Leucine", CCA=>"Proline", CCG=>"Proline", CAA=>"Glutamine", CAG=>"Glutamine", CGA=>"Arginine", CGG=>"Arginine", ATT=>"Isoleucine", ATC=>"Isoleucine", ACT=>"Threonine", ACC=>"Threonine", AAT=>"Asparagine", AAC=>"Asparagine", AGT=>"Serine", AGC=>"Serine", ATA=>"Methionine", ACA=>"Threonine", AAA=>"Lysine", AGA=>"Stop", ATG=>"Methionine", ACG=>"Threonine", AAG=>"Lysine", AGG=>"Stop", GTT=>"Valine", GTC=>"Valine", GCT=>"Alanine", GCC=>"Alanine", GAT=>"Aspartic acid", GAC=>"Aspartic acid", GGT=>"Glycine", GGC=>"Glycine", GTA=>"Valine", GTG=>"Valine", GCA=>"Alanine", GCG=>"Alanine", GAA=>"Glutamic acid", GAG=>"Glutamic acid", GGA=>"Glycine", GGG=>"Glycine");
our %codonr1m = (UUU=>"F", UUC=>"F", UCU=>"S", UCC=>"S", UAU=>"Y", UAC=>"Y", UGU=>"C", UGC=>"C", UUA=>"L", UCA=>"S", UAA=>"*", UGA=>"W", UUG=>"L", UCG=>"S", UAG=>"*", UGG=>"W", CUU=>"L", CUC=>"L", CCU=>"P", CCC=>"P", CAU=>"H", CAC=>"H", CGU=>"R", CGC=>"R", CUA=>"L", CUG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", AUU=>"I", AUC=>"I", ACU=>"T", ACC=>"T", AAU=>"N", AAC=>"N", AGU=>"S", AGC=>"S", AUA=>"M", ACA=>"T", AAA=>"K", AGA=>"*", AUG=>"M", ACG=>"T", AAG=>"K", AGG=>"*", GUU=>"V", GUC=>"V", GCU=>"A", GCC=>"A", GAU=>"D", GAC=>"D", GGU=>"G", GGC=>"G", GUA=>"V", GUG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");
our %codonr3m = (UUU=>"Phe", UUC=>"Phe", UCU=>"Ser", UCC=>"Ser", UAU=>"Tyr", UAC=>"Tyr", UGU=>"Cys", UGC=>"Cys", UUA=>"Leu", UCA=>"Ser", UAA=>"*", UGA=>"Trp", UUG=>"Leu", UCG=>"Ser", UAG=>"*", UGG=>"Trp", CUU=>"Leu", CUC=>"Leu", CCU=>"Pro", CCC=>"Pro", CAU=>"His", CAC=>"His", CGU=>"Arg", CGC=>"Arg", CUA=>"Leu", CUG=>"Leu", CCA=>"Pro", CCG=>"Pro", CAA=>"Gln", CAG=>"Gln", CGA=>"Arg", CGG=>"Arg", AUU=>"Ile", AUC=>"Ile", ACU=>"Thr", ACC=>"Thr", AAU=>"Asn", AAC=>"Asn", AGU=>"Ser", AGC=>"Ser", AUA=>"Met", ACA=>"Thr", AAA=>"Lys", AGA=>"*", AUG=>"Met", ACG=>"Thr", AAG=>"Lys", AGG=>"*", GUU=>"Val", GUC=>"Val", GCU=>"Ala", GCC=>"Ala", GAU=>"Asp", GAC=>"Asp", GGU=>"Gly", GGC=>"Gly", GUA=>"Val", GUG=>"Val", GCA=>"Ala", GCG=>"Ala", GAA=>"Glu", GAG=>"Glu", GGA=>"Gly", GGG=>"Gly");
our %codonrfullm = (UUU=>"Phenylalanine", UUC=>"Phenylalanine", UCU=>"Serine", UCC=>"Serine", UAU=>"Tyrosine", UAC=>"Tyrosine", UGU=>"Cysteine", UGC=>"Cysteine", UUA=>"Leucine", UCA=>"Serine", UAA=>"Stop", UGA=>"Tryptophan", UUG=>"Leucine", UCG=>"Serine", UAG=>"Stop", UGG=>"Tryptophan", CUU=>"Leucine", CUC=>"Leucine", CCU=>"Proline", CCC=>"Proline", CAU=>"Histidine", CAC=>"Histidine", CGU=>"Arginine", CGC=>"Arginine", CUA=>"Leucine", CUG=>"Leucine", CCA=>"Proline", CCG=>"Proline", CAA=>"Glutamine", CAG=>"Glutamine", CGA=>"Arginine", CGG=>"Arginine", AUU=>"Isoleucine", AUC=>"Isoleucine", ACU=>"Threonine", ACC=>"Threonine", AAU=>"Asparagine", AAC=>"Asparagine", AGU=>"Serine", AGC=>"Serine", AUA=>"Methionine", ACA=>"Threonine", AAA=>"Lysine", AGA=>"Stop", AUG=>"Methionine", ACG=>"Threonine", AAG=>"Lysine", AGG=>"Stop", GUU=>"Valine", GUC=>"Valine", GCU=>"Alanine", GCC=>"Alanine", GAU=>"Aspartic acid", GAC=>"Aspartic acid", GGU=>"Glycine", GGC=>"Glycine", GUA=>"Valine", GUG=>"Valine", GCA=>"Alanine", GCG=>"Alanine", GAA=>"Glutamic acid", GAG=>"Glutamic acid", GGA=>"Glycine", GGG=>"Glycine"); 		#"

our %iupac = (R=>'AG', Y=>'CT', S=>'GC', W=>'AT', K=>'GT', M=>'AC', A=>'AA', C=>'CC', G=>'GG', T=>'TT', B=>'CGT', D=>'AGT', H=>'ACT', V=>'ACG', N=>'ACGT', '.'=>'-', '-'=>'-');


# process program arguments, set up default values, check for errors, check for existence of db files
processArguments ();

$time and printerr "NOTICE: Current time (before execution) is ", scalar (localtime), "\n";

# print out output file names
if ($geneanno) {
	printerr "NOTICE: Output files are written to $outfile.variant_function, $outfile.exonic_variant_function\n";
} elsif ($regionanno) {
	printerr "NOTICE: Output file is written to $outfile.${buildver}_$dbtype1\n";
} elsif ($filter) {
	printerr "NOTICE: Output file with variants matching filtering criteria is written to $outfile.${buildver}_${dbtype1}_dropped, and output file with other variants is written to $outfile.${buildver}_${dbtype1}_filtered\n";
}

# check number of input lines in query, and adjust number of threads accordingly for gene-based annotation (when the user specifies -thread argument for multi-threading)
# my general observation is that when input line is less than 1 million, multi-threading does not really improve performance and should be disabled (--mingenelinecount controls this behavior)
# another observation is that when number of threads are too high, the program runs slower (--maxgenethread controls this); this is dependent on the I/O of the specific computing platform
my ($queryfile_line_count, $chunk_line_count);
if ($thread) {
	($queryfile_line_count, $chunk_line_count) = calculateChunkLine ($queryfile, $thread);
	if ($geneanno) {
		if ($queryfile_line_count < $mingenelinecount) {		#for gene-based annotation, only use threading if more than 1 million variants are present
			printerr ("NOTICE: threading is disabled for gene-based annotation on file with less than $mingenelinecount input lines\n");
			$thread = undef;
		}
		if ($thread and $thread > $maxgenethread) {				#do not use too many threads for gene-based annotation
			printerr ("NOTICE: number of threads is reduced to $maxgenethread (use --maxgenethread to change this behavior)\n");
			$thread = $maxgenethread;
			($queryfile_line_count, $chunk_line_count) = calculateChunkLine ($queryfile, $thread);	#re-calculate the two measures after adjusting thread
		}
	}
}

# begin the main program to process gene, region or filter-based annotation
if ($thread) {	
	my (@thr, $thr);
	for my $i (0 .. $thread-1) {
		my $start_line = $i * $chunk_line_count + 1;	#the $. linecount starts with 1
		my $end_line = $start_line + $chunk_line_count - 1;
		if ($end_line > $queryfile_line_count) {
			$end_line = $queryfile_line_count;
		}
		printerr ("NOTICE: Creating new threads for query line $start_line to $end_line\n");
		if ($geneanno) {
			$thr = threads->create(\&annotateQueryByGeneThread, "$outfile.variant_function.$i", "$outfile.exonic_variant_function.$i", 
			"$outfile.invalid_input.$i", $start_line, $end_line, $i);
		} elsif ($regionanno) {
			$thr = threads->create(\&annotateQueryByRegionThread, "$outfile.${buildver}_$dbtype1.$i", "$outfile.invalid_input.$i",
			$start_line, $end_line, $i);
		} elsif ($filter) {
			$thr = threads->create(\&filterQueryThread, "$outfile.${buildver}_${dbtype1}_filtered.$i", "$outfile.${buildver}_${dbtype1}_dropped.$i", 
				"$outfile.invalid_input.$i", $start_line, $end_line, $i);
		}
		push @thr, $thr;
	}
	for my $i (0 .. $thread-1) {
		$thr[$i]->join();		#wait for the corresponding thread to complete its execution
	}
	if ($geneanno) {
		combineFile ($outfile, "variant_function", "exonic_variant_function", "invalid_input");
	} elsif ($regionanno) {
		combineFile ($outfile, "${buildver}_$dbtype1", "invalid_input");
	} elsif ($filter) {
		combineFile ($outfile, "${buildver}_${dbtype1}_filtered", "${buildver}_${dbtype1}_dropped", "invalid_input");
	}
} else {
	if ($geneanno) {
		annotateQueryByGeneThread ("$outfile.variant_function", "$outfile.exonic_variant_function", "$outfile.invalid_input");
	} elsif ($regionanno) {
		annotateQueryByRegionThread ("$outfile.${buildver}_$dbtype1", "$outfile.invalid_input");
	} elsif ($filter) {
		filterQueryThread ("$outfile.${buildver}_${dbtype1}_filtered", "$outfile.${buildver}_${dbtype1}_dropped", "$outfile.invalid_input");
	} elsif ($downdb) {
		downloadDB ();
	}
}

# delete invlid_input file if it has zero size, otherwise print out the number of invalid lines
if (not $downdb) {
	if (not -z "$outfile.invalid_input") {
		printerr "NOTICE: Variants with invalid input format are written to $outfile.invalid_input\n";
	} else {
		unlink ("$outfile.invalid_input");
	}
}
$time and printerr "NOTICE: Current time (after execution) is ", scalar (localtime), "\n";


# combine serveral files produced by multiple threads into one single file, then delete the individual files
sub combineFile {
	my ($outfile, @suffix) = @_;
	for my $j (0 .. @suffix-1) {
		if (not rename ("$outfile.$suffix[$j].0", "$outfile.$suffix[$j]")) {
			die "Error: cannot rename $outfile.$suffix[$j].0 to $outfile.$suffix[$j]\n";
		}
	}
	for my $i (1 .. $thread-1) {
		my $msg = qx/cat --help 2>&1/ || '';		#collect the output of the system command
		if ($msg =~ m/^Usage/) {
			for my $j (0 .. @suffix-1) {
				system ("cat $outfile.$suffix[$j].$i >> $outfile.$suffix[$j]");		#if `cat` is available (in Linux/Unix, Apple Mac OS, etc), use it
			}
		} else {
			for my $j (0 .. @suffix-1) {
				appendFile ("$outfile.$suffix[$j].$i", "$outfile.$suffix[$j]");		#if `cat` is not available (in Windows), use appendFile() subroutine instead, which is slower
			}
		}
		for my $j (0 .. @suffix-1) {
			unlink ("$outfile.$suffix[$j].$i");
		}
	}
}

# Calculate how many query lines should be in each chunk, given the number of threads
sub calculateChunkLine {
	my ($queryfile, $thread) = @_;
	my $chunk_line_count;
	my $queryfile_line_count = `cat $queryfile | wc -l 2>&1`;		#if `cat` and `wc` are available (in Linux/Unix, Apple Mac OS, etc), use it
	chomp $queryfile_line_count;
	if ($queryfile_line_count =~ m/^\d+$/) {	#integer should be returned if the command finishes successfully
		printerr ("NOTICE: the queryfile $queryfile contains $queryfile_line_count lines\n");
	} else {									#if `cat` and `wc` are not available (in Windows), open the file and count the number of lines
		open (QUERY, $queryfile);
		$queryfile_line_count++ while (<QUERY>);
		close (QUERY);
		printerr ("NOTICE: the queryfile $queryfile contains $queryfile_line_count lines\n");
	}
	if ($queryfile_line_count % $thread == 0) {
		$chunk_line_count = $queryfile_line_count/$thread;
	} else {
		$chunk_line_count = int($queryfile_line_count/$thread)+1;	#max number of lines to be read from the input
	}
	return ($queryfile_line_count, $chunk_line_count);
}

# Append the content of one file to another, in case "cat" system command is not availalbe in the current operating system
sub appendFile {
	my ($file1, $file2) = @_;
	open (FH1, $file1) or die "Error: cannot read from file $file1: $!\n";
	open (FH2, ">>$file2") or die "Error: cannot append to file $file2: $!\n";
	while (<FH1>) {
		print FH2 $_;
	}
	close (FH1);
	close (FH2);
}

# Parse the arguments to this program and identify any potential errors in the arguments and set up appropriate default values
sub processArguments {
	my @command_line = @ARGV;		#command line argument
	GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'separate'=>\$separate,
	'batchsize=s'=>\$batchsize, 'dbtype=s'=>\$dbtype, 'neargene=i'=>\$neargene, 'genomebinsize=s'=>\$genomebinsize,
	'geneanno'=>\$geneanno, 'regionanno'=>\$regionanno, , 'filter'=>\$filter, 'downdb'=>\$downdb, 'buildver=s'=>\$buildver, 'score_threshold=f'=>\$score_threshold, 
	'normscore_threshold=i'=>\$normscore_threshold,	'minqueryfrac=f'=>\$minqueryfrac, 'expandbin=i'=>\$expandbin, 'splicing_threshold=i'=>\$splicing_threshold,
	'maf_threshold=f'=>\$maf_threshold, 'chromosome=s'=>\$chromosome, 'zerostart'=>\$zerostart, 'rawscore'=>\$rawscore, 'memfree=i'=>\$memfree, 
	'memtotal=i'=>\$memtotal, 'sift_threshold=f'=>\$sift_threshold, 'gff3dbfile=s'=>\$gff3dbfile, 'genericdbfile=s'=>\$genericdbfile, 'vcfdbfile=s'=>\$vcfdbfile,
	'time'=>\$time, 'wget!'=>\$wget, 'precedence=s'=>\$precedence, 'webfrom=s'=>\$webfrom, 'colsWanted=s'=>\$colsWanted, 'comment'=>\$comment,
	'scorecolumn=i'=>\$scorecolumn, 'poscolumn=s'=>\$poscolumn, 'transcript_function'=>\$transfun, 'exonsort'=>\$exonsort, 'avcolumn=s'=>\$avcolumn, 'bedfile=s'=>\$bedfile,
	'hgvs'=>\$hgvs, 'reverse'=>\$reverse, 'indexfilter_threshold=f'=>\$indexfilter_threshold, 'otherinfo'=>\$otherinfo, 
	'seq_padding=i'=>\$seq_padding, 'indel_splicing_threshold=i'=>\$indel_splicing_threshold, 'infoasscore'=>\$infoasscore, 'firstcodondel!'=>\$firstcodondel,
	'aamatrixfile=s'=>\$aamatrixfile, 'gff3attribute'=>\$gff3attr, 'exonicsplicing'=>\$exonicsplicing, 'infosep'=>\$infosep, 'dbm'=>\$dbm, 'idasscore'=>\$idasscore,
	'thread=i'=>\$thread, 'maxgenethread=i'=>\$maxgenethread) or pod2usage ();
	
	$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
	$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
	@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
	@ARGV == 2 or pod2usage ("Syntax error");

	($queryfile, $dbloc) = @ARGV;

	$dbloc =~ s/[\\\/]+$//;			#delete the trailing / or \ sign as part of the directory name
	$dbloc =~ s{^~([^/]*)} { $1 ? (getpwnam($1))[7] : ( $ENV{HOME} || $ENV{LOGDIR} || (getpwuid($>))[7])}ex;		#expand the ~ tilde to full path name (sometimes, annotate_variation.pl is not executed by the shell but excuted by other environments that do not recognize ~, so it is important to expand ~ to increase system compatibility)

	
	if (defined $batchsize) {
		$batchsize =~ s/k$/000/;
		$batchsize =~ s/m$/000000/;
		$batchsize =~ m/^\d+$/ or pod2usage ("Error: the --batchsize argument must be a positive integer (suffix of k or m is okay)");
	} else {
		$batchsize = 5_000_000;
	}
	if (defined $genomebinsize) {
		$genomebinsize =~ s/k$/000/;
		$genomebinsize =~ s/m$/000000/;
		$genomebinsize =~ m/^\d+$/ or pod2usage ("Error: the --genomebinsize argument must be a positive integer (suffix of k or m is okay)");
		$genomebinsize > 1000 or pod2suage ("Error: the --genomebinsize argument must be larger than 1000");
	} else {
		if ($geneanno) {
			$genomebinsize = 100_000;		#gene usually span large genomic regions
		} else {
			$genomebinsize = 10_000;		#MCE, TFBS, miRNA, etc are small genomic regions
		}
	}

	$verbose ||= 0;				#when it is not specified, it is zero
	$neargene ||= 1_000;		#for upstream/downstream annotation of variants, specify the distance threshold between variants and genes
	$expandbin ||= int(2_000_000/$genomebinsize);		#for gene-based annotations, when intergenic variants are found, expand to specified number of nearby bins to find closest genes
	$outfile ||= $queryfile;	#specify the prefix of output file names

	#set up log file
	if ($downdb) {
		if (not -d $dbloc) {
			mkdir ($dbloc) or die "Error: the directory $dbloc does not exist and cannot be created\n";
		}
		my $errfile = File::Spec->catfile ($dbloc, "annovar_downdb.log");
		open (LOG, ">$errfile") or die "Error: cannot write LOG information to log file $errfile: $!\n";
	} else {
		open (LOG, ">$outfile.log") or die "Error: cannot write LOG information to log file $outfile.log: $!\n";
	}
	print LOG "ANNOVAR Version:\n\t", q/$Date: 2020-06-07 23:56:37 -0400 (Sun,  7 Jun 2020) $/, "\n";
	print LOG "ANNOVAR Information:\n\tFor questions, comments, documentation, bug reports and program update, please visit http://www.openbioinformatics.org/annovar/\n";
	print LOG "ANNOVAR Command:\n\t$0 @command_line\n";
	print LOG "ANNOVAR Started:\n\t", scalar (localtime), "\n";
	
	my $num = 0;
	$geneanno and $num++;
	$downdb and $num++;
	$filter and $num++;
	$regionanno and $num++;
	$num <= 1 or pod2usage ("Error in argument: please specify only one of --geneanno, -regionanno, --downdb, --filter");
	if (not $num) {
		$geneanno++;
		printerr "NOTICE: The --geneanno operation is set to ON by default\n";
	}
		
	my %dbtype1 = ('gene'=>'refGene', 'refgene'=>'refGene', 'knowngene'=>'knownGene', 'ensgene'=>'ensGene', 'band'=>'cytoBand', 'cytoband'=>'cytoBand', 'tfbs'=>'tfbsConsSites', 'mirna'=>'wgRna',
			'mirnatarget'=>'targetScanS', 'segdup'=>'genomicSuperDups', 'omimgene'=>'omimGene', 'gwascatalog'=>'gwasCatalog',
			'1000g_ceu'=>'CEU.sites.2009_04', '1000g_yri'=>'YRI.sites.2009_04', '1000g_jptchb'=>'JPTCHB.sites.2009_04', 
			'1000g2010_ceu'=>'CEU.sites.2010_03', '1000g2010_yri'=>'YRI.sites.2010_03', '1000g2010_jptchb'=>'JPTCHB.sites.2010_03',			
			);		#for backward compatibility when these lower-case keywords were used in previous versions of ANNOVAR (the /^1000g(20\d\d)([a-z]{3})_([a-z]+)$/ pattern is handled below)

		
	if ($geneanno) {
		$dbtype ||= 'refGene';
		$dbtype1 = $dbtype1{$dbtype} || $dbtype;
	} elsif ($regionanno) {
		defined $dbtype or pod2usage ("Error in argument: please specify --dbtype (required for the --regionanno operation)");
		$dbtype1 = $dbtype1{$dbtype} || $dbtype;
		if ($dbtype =~ m/^mce(\d+)way/) {			#added 2010Feb16
			$dbtype1 = "phastConsElements$1way";
		}
		if ($dbtype1 eq 'gff3') {
			defined $gff3dbfile or pod2usage ("Error in argument: please specify --gff3dbfile for the --dbtype of 'gff3'");
		}
		if ($dbtype1 eq 'generic') {
			defined $genericdbfile or pod2usage ("Error in argument: please specify --genericdbfile for the --dbtype of 'generic'");
		}
		if ($dbtype1 eq 'bed') {
			defined $bedfile or pod2usage ("Error in argument: please specify --bedfile for the --dbtype of 'bed'");
		}
	} elsif ($filter) {
		defined $dbtype or pod2usage ("Error in argument: please specify --dbtype (required for the --filter operation)");
		#as of Feb 2012, I no longer check the validity of the database name for -filter operation, to give users the maximum amount of flexibility in designing and using their own favorite databases
		$dbtype =~ m/^avsift|1000g_(ceu|yri|jptchb)|1000g2010_(ceu|yri|jptchb)|1000g20\d\d[a-z]{3}_[a-z]+|[A-Z][A-Z][A-Z]\.sites.\d{4}_\d{2}|snp\d+\w+?|vcf|(ljb_[\w\+]+)|esp\d+_[\w]+|gnomad\w+|exac\w+|dbnsfp\w+$/ or print STDERR "NOTICE: the --dbtype $dbtype is assumed to be in generic ANNOVAR database format\n";
		
		$dbtype1 = $dbtype1{$dbtype} || $dbtype;
		
		if ($dbtype1 =~ m/^1000g(20\d\d)([a-z]{3})_([a-z]+)$/) {
			my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
			$dbtype1 = uc ($3) . '.sites.' . $1 . '_' . $monthhash{$2};
		}
		
		if ($dbtype1 eq 'generic') {
			defined $genericdbfile or pod2usage ("Error in argument: please specify --genericdbfile for the --dbtype of 'generic'");
		}
		if ($dbtype eq 'vcf') {
			defined $vcfdbfile or pod2usage ("Error in argument: please specify --vcfdbfile for the --dbtype of 'vcf'");
		}
	} elsif ($downdb) {
		defined $dbtype and pod2usage ("Error in argument: please do not specify --dbtype for the --downdb operation");
		$dbtype1 = $dbtype1{$queryfile} || $queryfile;
		if ($queryfile =~ m/^mce(\d+)way/) {			#added 2013may08
			$dbtype1 = "phastConsElements$1way";
		}
	}
	
	if (not $buildver) {
		$buildver = 'hg18';
		printerr "NOTICE: The --buildver is set as 'hg18' by default\n";
	}
	
	if (defined $score_threshold) {
		#$score_threshold >= 0 or pod2usage ("Error in argument: the --score_threshold must be a positive number or zero (you specified $score_threshold)");	#20130208: score_threshold can be anything now as long as it is a number
		$geneanno || $downdb and pod2usage ("Error in argument: the --score_threshold is not useful for --geneanno or --downdb operations");
	}
	if ($normscore_threshold) {
		$normscore_threshold >= 0 and $normscore_threshold <= 1000 or pod2usage ("Error in argument: the --normscore_threshold must be between 0 and 1000 (you specified $normscore_threshold)");
		$regionanno or pod2usage ("Error in argument: the --normscore_threshold is supported only for the --regionanno operation");
	}
	
	if ($zerostart) {
		pod2usage ("Error: the -zerostart argument is now obselete and will no longer be supported in ANNOVAR");
	}
	
	if (defined $sift_threshold) {
		$filter or pod2usage ("Error in argument: the --sift_threshold is supported only for the --filter operation");
		$dbtype1 eq 'avsift' or pod2usage ("Error in argument: the --sift_threshold argument can be used only if '--dbtype avsift' is used");
		$sift_threshold >= 0 and $sift_threshold <= 1 or pod2usage ("Error in argument: the --sift_threshold must be between 0 and 1 inclusive");
	} else {
		$dbtype1 eq 'avsift' and printerr "NOTICE: The --sift_threshold is set as 0.05 by default\n";
		$sift_threshold = 0.05;
	}
	
	if (defined $indexfilter_threshold) {
		$filter or pod2usage ("Error in argument: the --indexfilter_threshold is supported only for the --filter operation");
		$indexfilter_threshold >= 0 and $indexfilter_threshold <= 1 or pod2usage ("Error in argument: the --indexfilter_threshold must be between 0 and 1 inclusive");
	} else {
		$indexfilter_threshold = 0.9;
	}
	
	#operation-specific argument
	if (defined $splicing_threshold) {
		$geneanno or pod2usage ("Error in argument: the --splicing_threshold is supported only for the --geneanno operation");
	} else {
		$splicing_threshold = 2;	#for splicing annotation, specify the distance threshold between variants and exon/intron boundaries
	}
	if (defined $indel_splicing_threshold) {
		$geneanno or pod2usage ("Error: the --indel_splicing_threshold is supported only for the --geneanno operation");
	}
	else {
		$indel_splicing_threshold = $splicing_threshold;    #if not set, preserve original behavior
	}
	if (defined $maf_threshold) {
		$filter or pod2usage ("Error in argument: the --maf_threshold is supported only for the --filter operation");
		$dbtype =~ m/^1000g/ or pod2usage ("Error in argument: the --maf_threshold is supported only for 1000 Genomes Project data set (try -score_threshold instead)");
	} else {
		$maf_threshold = 0;		#for filter-based annotations on 1000 Genomes Project data, specify the MAF threshold to be used in filtering
	}
	if (defined $minqueryfrac) {
		$regionanno or pod2usage ("Error in argument: the --minqueryfrac is supported only for the --regionanno operation");
	} else {
		$minqueryfrac = 0;		#minimum query overlap to declare a "match" with database records
	}
	if (defined $gff3dbfile) {
		$dbtype eq 'gff3' or pod2usage ("Error in argument: the --gff3dbfile argument can be used only if '--dbtype gff3' is used");
		$geneanno or $regionanno or pod2usage ("Error in argument: the --gff3dbfile argument is supported only for the --geneanno or --regionanno operation");
	}
	if (defined $bedfile) {
		$dbtype eq 'bed' or pod2usage ("Error in argument: the --bedfile argument can be used only if '--dbtype bed' is used");
		$regionanno or pod2usage ("Error in argument: the --bedfile argument is supported only for the --regionanno operation");
	}
	if (defined $genericdbfile) {
		$filter or $regionanno or pod2usage ("Error in argument: the --genericdbfile argument is supported only for the --filter and -region operation");
	}
	if (defined $wget) {
		$downdb or pod2usage ("Error in argument: the --wget argument is supported only for the --downdb operation");
	} else {
		$wget = 1;			#by default, use wget for downloading files from Internet
	}
	if (defined $precedence) {
		$geneanno or pod2usage ("Error in argument: the --precedence argument is supported only for the --geneanno operation");
		if ($precedence =~ m/#/) {
			@precedence = split (/#/, $precedence);
		} else {
			@precedence = split (/,/, $precedence);
		}
		@precedence >= 2 or pod2usage ("Error in argument: the --precedence argument should be comma delimited");
		for my $i (0 .. @precedence-1) {
			$precedence[$i] =~ m/^(exonic|intronic|splicing|utr5|utr3|upstream|downstream|splicing|ncrna)$/ or pod2usage ("Error in argument: the --precedence argument contains invalid keywords (valid ones are exonic|intronic|splicing|utr5|utr3|upstream|downstream|splicing)");
		}
	}
	
	if (defined $colsWanted) {
		$regionanno or $filter or pod2usage ("Error in argument: the --colWanted argument is supported only for the --regionanno and --filter operation");
		if (lc $colsWanted eq 'all') {
			@colsWanted = ('all');
		} elsif (lc $colsWanted eq 'none') {
			@colsWanted = ('none');
		} else {
			if ($colsWanted =~ m/#/) {		#by default, -colsWanted are comma separated, but occasionally (for example, when use in table_annovar.pl, the user may want to use # to separate the columns)
				@colsWanted = split (/#/, $colsWanted);
			} else {
				@colsWanted = split (/,/, $colsWanted);
			}
			for my $i (0 .. @colsWanted-1) {
				$colsWanted[$i]=~m/^\d+$/ or pod2usage ("Error in argument: the --colsWanted argument ($colsWanted) must be a list of comma delimited numbers or be 'all' or be 'none'");
			}
		}
	}
	
	if (defined $scorecolumn) {
		$regionanno or pod2usage ("Error in argument: the --scorecolumn argument is supported only for the --regionanno operation");
	}
	
	if (defined $poscolumn) {
		$regionanno or pod2usage ("Error in argument: the --poscolumn argument is supported only for the --regionanno operation");
		$poscolumn =~ m/^\d+,\d+,\d+$/ or pod2usage ("Error in argument: the --poscolumn must be three integers separated by comma");
	}
	
	if ($exonsort) {
		$geneanno or pod2usage ("Error in argument: the --exonsort argument is supported only for the --geneanno operation");
	}
	
	if (defined $avcolumn) {
		if ($avcolumn =~ m/^\d+,\d+,\d+,\d+,\d+$/) {
			@avcolumn = split (/,/, $avcolumn);
		} elsif ($avcolumn =~ m/^\d+#\d+#\d+#\d+#\d+$/) {
			@avcolumn = split (/#/, $avcolumn);
		} else {
			pod2usage ("Error in argument: the --avcolumn argument must be five integer numbers separated by comma or #");
		}
		
		@avcolumn = map {$_-1} @avcolumn;
	} else {
		@avcolumn = (0..4);		#by default, the first five columns are the required AVINPUT information
	}
	
	if (defined $webfrom) {
		if ($webfrom ne 'ucsc' and $webfrom ne 'annovar') {
			$webfrom =~ m#^(http://|https://|ftp://)# or pod2usage ("Error: the --webfrom argument needs to be 'ucsc', 'annovar', or a URL");	#20191010: add https as an option
		}
	} else {
		$webfrom = 'ucsc';
	}
	
	# Padded output
	if ($seq_padding) {
		open $pad_fh, ">$outfile.seqpad" or die "Error: cannot write to output file $outfile.seqpad: $!\n";
		$cDNA_pad = $seq_padding * 3;
	}

	if ($infoasscore and $idasscore) {
		pod2usage ("Error in argument: you can specify either -infoasscore or -idasscore but not both");
	}
	
	$maf_threshold >= 0 and $maf_threshold <= 0.5 or pod2usage ("Error in argument: the --maf_threshold must be between 0 and 0.5 (you specified $maf_threshold)");
	$minqueryfrac >= 0 and $minqueryfrac <= 1 or pod2usage ("Error in argument: the --minqueryfrac must be between 0 and 1 (you specified $minqueryfrac)");
	$memfree and $memfree >= 100_000 || pod2usage ("Error in argument: the --memfree argument must be at least 100000 (in the order of kilobytes)");
	$memtotal and $memtotal >= 100_000 || pod2usage ("Error in argument: the --memtotal argument must be at least 100000 (in the order of kilobytes)");
	
	if ($chromosome) {
		my @chr;
		if ($chromosome =~ m/#/) {
			@chr = split (/#/, $chromosome);
		} else {
			@chr = split (/,/, $chromosome);
		}
		for my $i (0 .. @chr-1) {
			if ($chr[$i] =~ m/^(\d+)-(\d+)$/) {
				for my $j ($1 .. $2) {
					$valichr{$j}++;
				}
			} else {
				$valichr{$chr[$i]}++;
			}
		}
		printerr "NOTICE: These chromosomes in database will be examined: ", join (",", sort keys %valichr), "\n";
	}
	
	if (not defined $firstcodondel) {
		$firstcodondel = 1;
	}
	
	if (defined $aamatrixfile) {
		$geneanno or pod2usage ("Error in argument: the --aamatrix argument can be used only for gene-based annotation");
		$aamatrix = readAAMatrixFile ($aamatrixfile);
	}
	
	if ($gff3attr) {
		$dbtype eq 'gff3' or pod2usage ("Error in argument: the --gff3attr argument can be used only if '--dbtype gff3' is used");
	}

	if ($thread) {		#enable multi-threaded analysis (currently only filter-based annotation is supported)
		$downdb and pod2usage ("Error in argument: the --thread argument is not supported for --downdb operation");
		$thread > 0 or pod2usage ("Error in argument: the --thread argument must be a positive integer");
		if (not eval 'use threads; 1') {
			die "Error: your system does not support multi-threaded analysis. Please remove the --thread argument\n";
		}
	}
		
	if (not defined $maxgenethread) {
		$maxgenethread = 6;
	} else {
		$thread or pod2usage ("Error in argument: the --maxgenethread argument is supported only if --thread is set");
	}
	if (not defined $mingenelinecount) {
		$mingenelinecount = 1_000_000;
	} else {
		$thread or pod2usage ("Error in argument: the --mingenelinecount argument is supported only if --thread is set");
	}
}

# Read the AA matrxi file that scores the functional deleteriousness of a change from one amino acid to another. This is a simple tab-delimited file, with each row/column represent one amino acid
sub readAAMatrixFile {
	my ($aamatrixfile) = @_;
	my %aamatrix;
	open (MATRIX, $aamatrixfile) or die "Error: cannot read from aamatrixfile $aamatrixfile: $!\n";
	$_ = <MATRIX>;
	s/[\r\n]+$//;
	my @aa1 = split (/\t/, uc $_);
	shift @aa1;
	@aa1 == 20 or die "Error: invalid first line found in aamatrixfile (21 tab-delimited fields expected): <$_>\n";
	for my $i (0 .. @aa1-1) {
		$aa1[$i] =~ m/^[SRLPTAVGIFYCHQNKDEMW]$/ or die "Error: invalid amino acid identifier found in aamatrixfile (SRLPTAVGIFYCHQNKDEMW expected): <$aa1[$i]>\n";
	}
	while (<MATRIX>) {
		s/[\r\n]+$//;
		my @aa2 = split (/\t/, uc $_);
		my $aa2 = shift @aa2;
		@aa2 == 20 or die "Error: invalid line found in aamatrixfile (21 tab-delimited fields expected): <$_>\n";
		$aa2 =~ m/^[SRLPTAVGIFYCHQNKDEMW]$/ or die "Error: invalid amino acid identifier found in aamatrixfile (SRLPTAVGIFYCHQNKDEMW expected): <$_>\n";
		for my $j (0 .. @aa2-1) {
			$aamatrix{$aa2.$aa1[$j]} = $aa2[$j];
		}
	}
	close (MATRIX);
	return (\%aamatrix);
}

# Parse ANNOVAR input lines and detect whether invalid input is present, if not, returns the five ANNOVAR input columns
sub detectInvalidInput {
	my ($line) = @_;
	my $invalid = 0;
	my @nextline = split (/\s+/, $line);
	my ($chr, $start, $end, $ref, $obs) = @nextline[@avcolumn];
	if ( not (defined $chr and defined $start and defined $end and defined $ref and defined $obs)) {
		$invalid++;
	} else {
		($ref, $obs) = (uc $ref, uc $obs);
		$zerostart and $start++;
		$chr =~ s/^chr//;
		$ref =~ s/^\-+$/-/;		#sometimes users use -- as the reference allele for whatever reason
		$obs =~ s/^\-+$/-/;		#sometimes users use -- as the alternative allele for whatever reason
		if ($chr =~ m/[^\w\.\-]/ or $start =~ m/[^\d]/ or $end =~ m/[^\d]/) {		#chr name could contain . (example: GL000212.1, or Zv9_NA###, or ERCC-3324
			$invalid++;
		} elsif ($ref eq '-' and $obs eq '-' 		#both are empty allele
			or $ref =~ m/[^ACTG0\-]/ 		#non-standard nucleotide code
			or $obs =~ m/[^ACGT0\-]/ 		#non-standard nucleotide code
			or $start =~ m/[^\d]/ 			#start is not a number
			or $end =~ m/[^\d]/ 			#end is not a number
			or $start > $end			#start is more than end
			or $ref ne '0' and $end-$start+1 != length ($ref) 	#length mismatch with ref
			or $ref eq '-' and $start != $end	#length mismatch for insertion
			) {
			$invalid++;
		}
	}
	return ($invalid, $chr, $start, $end, $ref, $obs);
}

# Control gene-based annotation by submitting a batch request to processNextQueryBatchByGeneThread(), with the support for multi-threading if $start_line, $end_line and $cur_thread is specified
sub annotateQueryByGeneThread {
	my ($varfun_file, $exvarfun_file, $invalid_file, $start_line, $end_line, $cur_thread) = @_;
	my ($VARFUN, $EXVARFUN, $INVALID, $QUERY);
	open ($VARFUN, ">$varfun_file") or die "Error: cannot write to output file $varfun_file: $!\n";
	open ($EXVARFUN, ">$exvarfun_file") or die "Error: cannot write to output file $exvarfun_file: $!\n";
	open ($INVALID, ">$invalid_file") or die "Error: cannot write to output file $invalid_file: $!\n";	
	
	#my ($totalquerycount, $totalinvalidcount, $batchcount) = qw/0 0 1/;
	open ($QUERY, $queryfile) or die "Error: cannot read from --queryfile ($queryfile): $!\n";

	my ($genedb, $geneidmap, $cdslen, $mrnalen) = readUCSCGeneAnnotation ($dbloc);
	
	my (@variant, $filedone, $batchdone);
	my ($linecount, $batchlinecount, $invalid, $invalidcount) = (0, 0);
	my ($chr, $start, $end, $ref, $obs, $info);
	
	if ($thread and $start_line) {
		for (1 .. $start_line-1) {
			$_ = <$QUERY>;
		}
	}
	while (1) {
		$_ = <$QUERY>;
		if ($thread) {
			if ($. > $end_line) {
				$_ = undef;		#do not read any line any more
			}
		}
		
		if (not defined $_) {
			$filedone++;
		} else {
			$_ =~ s/[\r\n]+$//;
		
			if ($_ =~ m/^#/ and $comment) {			#comment line start with #, do not include this is $linecount
				print $VARFUN "#comment\t$comment\t$_\n";
				next;
			}
			
			($invalid, $chr, $start, $end, $ref, $obs) = detectInvalidInput ($_);
			
			if ($invalid) {
				print $INVALID $_, "\n";			#invalid record found
				$invalidcount++;
				next;
			}
			push @variant, [$chr, $start, $end, $ref, $obs, $., $_];
			
			$linecount++;		#does not include comment line or invalid line
			$batchlinecount++;
			if ($batchlinecount == $batchsize) {
				$batchdone++;
			}
		}
		
		if ($filedone or $batchdone) {
			@variant and printerr "NOTICE: Processing next batch with ${\(scalar @variant)} unique variants in $batchlinecount input lines\n";
			processNextQueryBatchByGeneThread (\@variant, $VARFUN, $EXVARFUN, $INVALID, $genedb, $geneidmap, $cdslen, $mrnalen);
			@variant = ();
			$batchlinecount = 0;
			$batchdone = 0;
		}
		if ($filedone) {
			close ($VARFUN);
			close ($EXVARFUN);
			close ($INVALID);
			close ($QUERY);
			last;
		}
	}
}

# Perform the actual gene-based annotation on a batch of input variants
sub processNextQueryBatchByGeneThread {
	my ($variant, $VARFUN, $EXVARFUN, $INVALID, $genedb, $geneidmap, $cdslen, $mrnalen) = @_;
	my (%refseqvar);
	
	my ($name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend, $name2);
	
	for my $i (0 .. @$variant-1) {
		my ($chr, $start, $end, $ref, $obs, $curlinecount, $nextline) = @{$variant->[$i]};
	
	
		
		my (%intronic, %utr5, %utr3, %exonic, %upstream, %downstream, %ncrna, %intergenic, %splicing, %splicing_anno, %newutr5, %newutr3);
		my $foundgenic;						#variant found in genic region (between start and end position of a gene in genome)
		my ($distl, $distr, $genel, $gener);			#for intergenic variant, the distance and gene name to the left and right side of gene
		my ($distup, $distdown);				#for upstream/downstream variants, the distance to the transcript
		my $bin1 = int ($start/$genomebinsize)-1;		#start bin
		$bin1 < 0 and $bin1=0;
		my $bin2 = int ($end/$genomebinsize)+1;			#end bin (usually same as start bin, unless the query is really big that spans multiple megabases)
		
		while (not exists $genedb->{$chr, $bin1} and $bin1 > int ($start/$genomebinsize)-$expandbin) {		#examine at least 5 bins (by default 5Mb) to the left to make sure that a gene is found in the bin
			$bin1 > 0 or last;
			$bin1--;
		}
		
		while (not exists $genedb->{$chr, $bin2} and $bin2 < int ($end/$genomebinsize)+$expandbin) {		#examine at least 5 bins (by default 5Mb) to the right to make sure that a gene is found in the bin
			$bin2++;
		}

		my (%seen);
		for my $nextbin ($bin1 .. $bin2) {
			exists $genedb->{$chr, $nextbin} or next;		#this genome bin has no annotated gene (a complete intergenic region)
			for my $nextgene (@{$genedb->{$chr, $nextbin}}) {	#when $genedb->{$chr, $nextbin} is undefined, this automatically create an array!!!
				($name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend, $name2) = @$nextgene;
				defined $name2 or printerr "WARNING: name2 field is not provided for transcript $name (start=$txstart end=$txend)\n" and $name2='';
				$seen{$name, $txstart} and next;		#name and txstart uniquely identify a transcript and chromosome position (sometimes same transcript may map to two nearby positions, such as nearby segmental duplications)
				$seen{$name, $txstart}++;			#a transcript may be in two adjacent bins, so if one is already scanned, there is no need to work on it again
				my $current_ncRNA;
				
				if ($transfun) {					#variant_function output contains transcript name, rather than gene name
					$name2 = $name;
					$name2 =~ s/#\w+#\d+//;			#20140711 (issue:multimap) multimap issue
				}
				
				if (not $foundgenic or $separate) {				#this variant has not hit a genic region yet (and -separate is not specified, 20170712)
					if ($start > $txend) {
						defined $distl or $distl = $start-$txend and $genel=$name2;
						$distl > $start-$txend and $distl = $start-$txend and $genel=$name2;	#identify left closest gene
					}
		
					if ($end < $txstart) {
						defined $distr or $distr = $txstart-$end and $gener=$name2;
						$distr > $txstart-$end and $distr = $txstart-$end and $gener=$name2;	#identify right closest gene
					}
				}
				
				if ($end < $txstart) {
					#query ---
					#gene		<-*----*->
					not $separate and $foundgenic and last;			#if found a genic annotation already, end the search of the bins (however, if -separate is specified, we should not end the search)
					if ($end > $txstart - $neargene) {
						if ($dbstrand eq '+') {
							$upstream{$name2}++;
							$distup = $distr;
						} else {
							$downstream{$name2}++;
							$distdown = $distr;
						}
					} else {
						last;						#if transcript is too far away from end, end the search of the bins
					}
				} elsif ($start > $txend) {
					#query            ---
					#gene  <-*----*->
					if ($separate || not $foundgenic and $start < $txend + $neargene) {	#similar to above, if separate is specified, we still record up/downstream regardless of foundgenic
						if ($dbstrand eq '+') {
							$downstream{$name2}++;
							$distdown = $distl;
						} else {
							$upstream{$name2}++;
							$distup = $distl;
						}
					}
				#} elsif ($cdsstart == $cdsend+1) {				#non-coding RNA (could be microRNA, or could be due to lack of CDS annotation for mRNA such as NR_026730 or BC039000). Previously we already did cdsstart++ so here the cdsstart is more than cdsend
				#	if ($start >= $txstart and $start <= $txend or $end >= $txstart and $end <= $txend or $start <= $txstart and $end >= $txend) {
				#		$ncrna{$name2}++;
				#		$foundgenic++;
				#	}
				} else {							#query overlaps with coding region of gene
					
					#change in 2011jul24: handle ncRNA and protein coding gene together but with the ncRNA flag when printing out results in the future
					if ($cdsstart == $cdsend+1) {				#non-coding RNA (could be microRNA, or could be due to lack of CDS annotation for mRNA such as NR_026730 or BC039000). Previously we already did cdsstart++ so here the cdsstart is more than cdsend
						if ($start >= $txstart and $start <= $txend or $end >= $txstart and $end <= $txend or $start <= $txstart and $end >= $txend) {
							$ncrna{$name2}++;
							$foundgenic++;
						}
						
						#now treat this ncRNA as if it is a protein-coding gene
						($cdsstart, $cdsend) = ($txstart, $txend);
						$current_ncRNA++;		#current transcript is a noncoding transcript
					}
					
					
					
					my ($lenintron, $lenexon) = (0, 0);			#cumulative intron/exon length at a given exon
					my ($rcdsstart, $rcdsend, $rvarstart, $rvarend);			#start of coding and variant in reference mRNA sequence
					my @exonstart = @$exonstart;
					my @exonend = @$exonend;
					my $foundexonic;
					my $utr_canno;						#cDNA level annotation for UTR variants
					my ($refcut, $obscut);						#alternative obs by cutting the original obs (when a block substitution overlaps with both exon and intron, we should only use the exonic portion to calculate coding change
					if ($dbstrand eq '+') {					#forward strand, search from left to right (first exon to last exon)
						for my $k (0 .. @exonstart-1) {
							$k and $lenintron += ($exonstart[$k]-$exonend[$k-1]-1);		#calculate cumulative intron length
							$lenexon += ($exonend[$k]-$exonstart[$k]+1);
							if ($cdsstart >= $exonstart[$k]) {				#calculate CDS start accurately by considering intron length
								$rcdsstart = $cdsstart-$txstart-$lenintron+1;
								
								if ($cdsstart <= $exonend[$k]) {	#CDS start is within this exon
									$lenexon = ($exonend[$k]-$cdsstart+1);
								} else {				#CDS start is in previous exon
									#$lenexon += ($exonend[$k]-$exonstart[$k]+1);
								}
							}
							if ($cdsend >= $exonstart[$k]) {
								$rcdsend = $cdsend-$txstart-$lenintron+1;
							}
							
							if ($exonicsplicing) {			#when -exonicsplicing argument is set, exonic variants near exon/intron boundary is reported as "exonic;splicing" variant
								#splicing calculation (changed 2012may24)
								if (@exonstart != 1) {
									if ( $k == 0 and $start >= $exonend[$k]-$splicing_threshold+1 and $start <= $exonend[$k]+$splicing_threshold) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]+$splicing_threshold-1) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]+$splicing_threshold-1 or $start >= $exonend[$k]-$splicing_threshold+1 and $start <= $exonend[$k]+$splicing_threshold) {
											$splicing{$name2}++;		#when query start site is close to exon start or exon end
										}
									}
									
									if ($k == 0 and $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1 or $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
											$splicing{$name2}++;		#when query end site is close to exon start or exon end
										}
									}
									
									if ($k == 0 and $start <= $exonend[$k] and $end >= $exonend[$k]) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $start <= $exonstart[$k] and $end>=$exonstart[$k]) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($start <= $exonstart[$k] and $end>=$exonstart[$k] or $start <= $exonend[$k] and $end >= $exonend[$k]) {
											$splicing{$name2}++;		#when query encompass the exon/intron boundary
										}
									}
								}
							} else {
								#splicing calculation (changed 2013feb10)
								#splicing calculation (changed again 2013feb21 per Mitsuhiro Komura suggestion)
								if (@exonstart != 1) {
									if ( $k == 0 and $start >= $exonend[$k]+1 and $start <= $exonend[$k]+$splicing_threshold) {			#first exon
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]-1) {	#last exon
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {										#middle exon
										#if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]-1 or 
										#$start >= $exonend[$k]+1 and $start <= $exonend[$k]+$splicing_threshold or
										#$start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]-1 or
										#$end >= $exonend[$k]+1 and $end <= $exonend[$k]+$splicing_threshold
										if (hasOverlap ($start, $end, $exonstart[$k]-$splicing_threshold, $exonstart[$k]-1) or
										hasOverlap ($start, $end, $exonend[$k]+1, $exonend[$k]+$splicing_threshold)
										) {
											$splicing{$name2}++;		#when query start site is close to exon start or exon end
										}
									}
								}
							}
							
							
							
							#if name2 is already a splicing variant, but its detailed annotation (like c150-2A>G) is not available, and if this splicing leads to amino acid change (rather than UTR change)
							if ($splicing{$name2} and $start==$end and $start>=$cdsstart) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
									$lenexon -= ($exonend[$k]-$exonstart[$k]);		#formerly "$exonend[$k]-$exonstart[$k]+1"; changed this line 2011oct01 given German's comments. The coding portion for acceptor site should be added by one, but for donor site should be fine
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:c.$lenexon-" . ($exonstart[$k]-$start) . "$ref>$obs,";
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									#-------<---->-*--------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:c.$lenexon+" . ($start-$exonend[$k]) . "$ref>$obs,";
								}
							}
							#if name2 is already a splicing variant, but its detailed annotation (like c150-2A>G) is not available, and if it leads to UTR change
							if ($splicing{$name2} and $start==$end and $start<$cdsstart) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:UTR5,";
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									#-------<---->-*--------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:UTR5,";
								}
							} elsif ($splicing{$name2} and $start==$end and $start>$cdsend) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:UTR3,";
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									#-------<---->-*--------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:UTR3,";
								}
							}
							#20180404 IBM change: if name2 is already a splicing variant, but it is an indel, we will just use r.spl (original HGVS (http://www.hgvs.org/mutnomen/examplesRNA.html) or the updated definition (http://varnomen.hgvs.org/recommendations/RNA/variant/splicing/))
							if ($splicing{$name2} and $start!=$end) {
								if ($end >= $exonstart[$k]-$splicing_threshold and $end < $exonstart[$k]) {
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:r.spl,";
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									$splicing_anno{$name2} .= "$name:exon${\($k+1)}:r.spl,";
								}
							}
							
							if ($start < $exonstart[$k]) {
								if ($end >= $exonstart[$k]) {	#exonic 
									$rvarstart = $exonstart[$k]-$txstart-$lenintron+1;
									
									for my $m ($k .. @exonstart-1) {
										$m > $k and $lenintron += ($exonstart[$m]-$exonend[$m-1]-1);
										if ($end < $exonstart[$m]) {
											#query           --------
											#gene     <--**---******---****---->
											$rvarend = $exonend[$m-1]-$txstart-$lenintron+1 + ($exonstart[$m]-$exonend[$m-1]-1);
											last;
										} elsif ($end <= $exonend[$m]) {
											#query           -----------
											#gene     <--**---******---****---->
											$rvarend = $end-$txstart-$lenintron+1;
											last;
										}
									}
									if (not defined $rvarend) {
										$rvarend = $txend-$txstart-$lenintron+1;		#if this value is longer than transcript length, it suggest whole gene deletion
									}
									
									$refcut = substr ($ref, -($rvarend-$rvarstart+1));	#20170830
									$obscut = substr ($obs, -($rvarend-$rvarstart+1));	#20170830
									
									#here the trick begins to differentiate UTR versus coding exonic
									if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
										#query  ----
										#gene     <--*---*->
										#$utr5{$name2}++;		#positive strand for UTR5
										
										$utr_canno = "$name:c." . ($rvarstart-$rcdsstart) . '_' . ($rvarend-$rcdsstart) . "delins$obs";	#this is an indel
										$utr5{$name2} .= "$utr_canno,";
									} elsif ($start > $cdsend) {
										#query             ----
										#gene     <--*---*->
										#$utr3{$name2}++;		#positive strand for UTR3
										
										$utr_canno = "$name:c.*" . ($rvarstart-$rcdsend) . '_*' . ($rvarend-$rcdsend) . "delins$obs";
										$utr3{$name2} .= "$utr_canno,";
									} else {									
										$exonic{$name2}++;
										not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '+', $curlinecount, $k+1, $nextline, $refcut, $obscut];	#refseq CDS start, refseq variant start. obs is non-zero (obs is specified by user)
									}
									$foundgenic++;
									last;
								} elsif ($k and $start > $exonend[$k-1]) {	#intronic
									$intronic{$name2}++;
									$foundgenic++;
									last;
								}
							} elsif ($start <= $exonend[$k]) {	#exonic
								$rvarstart = $start-$txstart-$lenintron+1;
								
								for my $m ($k .. @exonstart-1) {
									$m > $k and $lenintron += ($exonstart[$m]-$exonend[$m-1]-1);
									if ($end < $exonstart[$m]) {
										#query              ------
										#gene     <--**---******---****---->
										$rvarend = $exonend[$m-1]-$txstart-$lenintron+1 + ($exonstart[$m]-$exonend[$m-1]-1);
										
										$refcut = substr ($ref, 0, $rvarend-$rvarstart+1);	#20170830 update obscut since the end is in intronic region
										$obscut = substr ($obs, 0, $rvarend-$rvarstart+1);	#20170830 update obscut since the end is in intronic region
										
										last;
									} elsif ($end <= $exonend[$m]) {
										#query           -----------
										#gene     <--**---******---****---->
										$rvarend = $end-$txstart-$lenintron+1;
										last;
									}
								}
								if (not defined $rvarend) {
									$rvarend = $txend-$txstart-$lenintron+1;		#if this value is longer than transcript length, it suggest whole gene deletion
								}
								
								#here is the trick begins to differentiate UTR versus coding exonic
								if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
									#query  ----
									#gene     <--*---*->
									#$utr5{$name2}++;		#positive strand for UTR5
									
									if ($start==$end and $ref ne '-' and $obs ne '-') {	#SNV
										$utr_canno = "$name:c." . ($rvarstart-$rcdsstart) . "$ref>$obs";
									} elsif ($start==$end and $ref eq '-' and $obs ne '-') {	#insertion
										$utr_canno = "$name:c." . ($rvarstart-$rcdsstart) . '_' . ($rvarstart-$rcdsstart+1) . "ins$obs";
									} elsif ($start==$end and $ref ne '-' and $obs eq '-') {	#single-base deletion
										$utr_canno = "$name:c." . ($rvarstart-$rcdsstart) . "del$obs";
									} elsif ($ref ne '-' and $obs eq '-') {		#multi-base deletion
										$utr_canno = "$name:c." . ($rvarstart-$rcdsstart) . '_' . ($rvarend-$rcdsstart) . "del$obs";
									} else {
										$utr_canno = "$name:c." . ($rvarstart-$rcdsstart) . '_' . ($rvarend-$rcdsstart) . "delins$obs";
									}
									$utr5{$name2} .= "$utr_canno,";
									
								} elsif ($start > $cdsend) {
									#query             ----
									#gene     <--*---*->
									#$utr3{$name2}++;		#positive strand for UTR3
									
									if ($start==$end and $ref ne '-' and $obs ne '-') {	#SNV
										$utr_canno = "$name:c.*" . ($rvarstart-$rcdsend) . "$ref>$obs";
									} elsif ($start==$end and $ref eq '-' and $obs ne '-') {	#insertion
										$utr_canno = "$name:c.*" . ($rvarstart-$rcdsend) . '_*' . ($rvarstart-$rcdsend+1) . "ins$obs";
									} elsif ($start==$end and $ref ne '-' and $obs eq '-') {	#single-base deletion
										$utr_canno = "$name:c.*" . ($rvarstart-$rcdsend) . "del$ref";
									} elsif ($ref ne '-' and $obs eq '-') {		#multi-base deletion
										$utr_canno = "$name:c.*" . ($rvarstart-$rcdsend) . '_*' . ($rvarend-$rcdsend) . "del$ref";
									} else {
										$utr_canno = "$name:c.*" . ($rvarstart-$rcdsend) . '_*' . ($rvarend-$rcdsend) . "delins$obs";
									}
									$utr3{$name2} .= "$utr_canno,";
									
								} else {
									$exonic{$name2}++;
									not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '+', $curlinecount, $k+1, $nextline, $refcut, $obscut];		#queryindex, refseq CDS start, refseq variant start
								}
								$foundgenic++;
								last;
							}
						}
					} elsif ($dbstrand eq '-') {		#process negative strand (in the future, this should be fused to the paragraph above for positive strands; for now, I keep them separate for easier debugging)
						
						for (my $k = @exonstart-1; $k>=0; $k--) {
							$k < @exonstart-1 and $lenintron += ($exonstart[$k+1]-$exonend[$k]-1);		#length of intron before current exon
							$lenexon += ($exonend[$k]-$exonstart[$k]+1);					#current cDNA position (for splicing calculation)
							if ($cdsend <= $exonend[$k]) {		#calculate CDS start accurately by considering intron length
								$rcdsstart = $txend-$cdsend-$lenintron+1;
								
								if ($cdsend >= $exonstart[$k]) {	#CDS start within this exon
									$lenexon = ($cdsend-$exonstart[$k]+1);		#reset lenexon since start is found
								} else {				#CDS start in prevous exon
									#$lenexon += ($exonend[$k]-$exonstart[$k]+1);
								}
							}
							if ($cdsstart <= $exonend[$k]) {
								#$rcdsend = $cdsend-$txstart-$lenintron+1;
								$rcdsend = $txend-$cdsstart-$lenintron+1;
							}

							if ($exonicsplicing) {			#when -exonicsplicing argument is set, exonic variants near exon/intron boundary is reported as "exonic;splicing" variant
								#splicing calculation (changed 2012may24)
								if (@exonstart != 1) {
									if ( $k == 0 and $start >= $exonend[$k]-$splicing_threshold+1 and $start <= $exonend[$k]+$splicing_threshold) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]+$splicing_threshold-1) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]+$splicing_threshold-1 or $start >= $exonend[$k]-$splicing_threshold+1 and $start <= $exonend[$k]+$splicing_threshold) {
											$splicing{$name2}++;		#when query start site is close to exon start or exon end
										}
									}
									
									if ($k == 0 and $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]+$splicing_threshold-1 or $end >= $exonend[$k]-$splicing_threshold+1 and $end <= $exonend[$k]+$splicing_threshold) {
											$splicing{$name2}++;		#when query end site is close to exon start or exon end
										}
									}
									
									if ($k == 0 and $start <= $exonend[$k] and $end >= $exonend[$k]) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $start <= $exonstart[$k] and $end>=$exonstart[$k]) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										if ($start <= $exonstart[$k] and $end>=$exonstart[$k] or $start <= $exonend[$k] and $end >= $exonend[$k]) {
											$splicing{$name2}++;		#when query encompass the exon/intron boundary
										}
									}
								}
							} else {
								#splicing calculation (changed 2013feb10)
								#splicing calculation (changed again 2013feb21 per Mitsuhiro Komura suggestion)
								if (@exonstart != 1) {
									if ( $k == 0 and $start >= $exonend[$k]+1 and $start <= $exonend[$k]+$splicing_threshold) {
										$splicing{$name2}++;
									} elsif ($k == @exonstart-1 and $end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]-1) {
										$splicing{$name2}++;
									} elsif ($k and $k < @exonstart-1) {
										#if ($end >= $exonstart[$k]-$splicing_threshold and $end <= $exonstart[$k]-1 or 
										#$start >= $exonend[$k]+1 and $start <= $exonend[$k]+$splicing_threshold or
										#$start >= $exonstart[$k]-$splicing_threshold and $start <= $exonstart[$k]-1 or
										#$end >= $exonend[$k]+1 and $end <= $exonend[$k]+$splicing_threshold
										if (hasOverlap ($start, $end, $exonstart[$k]-$splicing_threshold, $exonstart[$k]-1) or
										hasOverlap ($start, $end, $exonend[$k]+1, $exonend[$k]+$splicing_threshold)
										) {
											$splicing{$name2}++;		#when query start site is close to exon start or exon end
										}
									}
								}
							}
							
							#if name2 is already a splicing variant, but its detailed annotation (like c150-2A>G) is not available, and if this splicing leads to amino acid change (rather than UTR change)
							if ($splicing{$name2} and $start==$end and $start<=$cdsend) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-1-$k+1)}:c.$lenexon+" . ($exonstart[$k]-$start) . revcom($ref) . '>' . revcom ($obs) . ',';	#change ${\(@exonstart-$k+1)} to ${\(@exonstart-1-$k+1)} per Contini Elisa 20150604
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									#-------<---->-*--------<-->-------<------>----
									$lenexon -= ($exonend[$k]-$exonstart[$k]);	#formerly "$exonend[$k]-$exonstart[$k]+1"; changed this line 2011oct01 given German's comments. The coding portion for acceptor site should be added by one, but for donor site should be fine
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-1-$k+1)}:c.$lenexon-" . ($start-$exonend[$k]) . revcom($ref) . '>' . revcom($obs) . ',';	#change ${\(@exonstart-$k+1)} to ${\(@exonstart-1-$k+1)} per Contini Elisa 20150604
								}
							}
							#if name2 is already a splicing variant, but its detailed annotation (like c150-2A>G) is not available, and if it leads to UTR change
							if ($splicing{$name2} and $start==$end and $start>$cdsend) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-1-$k+1)}:UTR5,";	#change ${\(@exonstart-$k+1)} to ${\(@exonstart-1-$k+1)} per Contini Elisa 20150604
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									#-------<---->-*--------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-1-$k+1)}:UTR5,";	#change ${\(@exonstart-$k+1)} to ${\(@exonstart-1-$k+1)} per Contini Elisa 20150604
								}
							} elsif ($splicing{$name2} and $start==$end and $start<$cdsstart) {
								if ($start >= $exonstart[$k]-$splicing_threshold and $start < $exonstart[$k]) {
									#------*-<---->---------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-1-$k+1)}:UTR3,";	#change ${\(@exonstart-$k+1)} to ${\(@exonstart-1-$k+1)} per Contini Elisa 20150604
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									#-------<---->-*--------<-->-------<------>----
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-1-$k+1)}:UTR3,";	#change ${\(@exonstart-$k+1)} to ${\(@exonstart-1-$k+1)} per Contini Elisa 20150604
								}
							}
							#20180404 IBM change: if name2 is already a splicing variant, but it is an indel, we will just use r.spl (specified in both original HGVS (http://www.hgvs.org/mutnomen/examplesRNA.html) or the updated definition (http://varnomen.hgvs.org/recommendations/RNA/variant/splicing/))
							if ($splicing{$name2} and $start!=$end) {
								if ($end >= $exonstart[$k]-$splicing_threshold and $end < $exonstart[$k]) {
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-1-$k+1)}:r.spl,";
								} elsif ($start > $exonend[$k] and $start <= $exonend[$k]+$splicing_threshold) {
									$splicing_anno{$name2} .= "$name:exon${\(@exonstart-1-$k+1)}:r.spl,";
								}
								
							}
							
							#20170830: a major issue is that some mutations span intron/exon boundaries and we should only use the exonic portion as $obs, when calculating amino acid changes
							#an exmaple is chr9:135772997-135772998CC>TT: we should only use T as the $obs, not TT
							#however, the current code calculate exonic_variant_function using $obs extracted from $line, making this adjustment complicated
							#after careful thought, I think the best way forward (with mimimal change to code) is to add a $obscut field in the array to exonic_variant_function. when it is present, use this rather than $obs to calculate c. and p. changes
							if ($end > $exonend[$k]) {
								if ($start <= $exonend[$k]) {
									$rvarstart = $txend-$exonend[$k]-$lenintron+1;
									
									for (my $m = $k; $m >= 0; $m--) {
										$m < $k and $lenintron += ($exonstart[$m+1]-$exonend[$m]-1);
										if ($start > $exonend[$m]) {
											#query           --------
											#gene     <--**---******---****---->
											#$rvarend = $txend-$exonstart[$m]-$lenintron+1 - ($exonstart[$m+1]-$exonend[$m]-1);	#commented out 2011feb18
											$rvarend = $txend-$exonstart[$m+1]+1-$lenintron + ($exonstart[$m+1]-$exonend[$m]-1);	#fixed this 2011feb18
											last;		#finsih the cycle!!!!!!!!!!!!!!!!!!!
										} elsif ($start >= $exonstart[$m]) {		#start within exons
											#query               ----
											#gene     <--**---******---****---->
											$rvarend = $txend-$start-$lenintron+1;
											last;
										}
									}
									
									if (not defined $rvarend) {				#if rvarend is not found, then the whole tail of gene is covered
										$rvarend = $txend-$txstart-$lenintron+1;
									}
									
									$refcut = substr ($ref, 0, $rvarend-$rvarstart+1);	#20170830
									$obscut = substr ($obs, 0, $rvarend-$rvarstart+1);	#20170830
									
									#here is the trick begins to differentiate UTR versus coding exonic
									if ($end < $cdsstart) {					#usually disrupt/change 5' UTR region, unless the UTR per se is also separated by introns
										#query  ----
										#gene     <--*---*->
										#$utr3{$name2}++;		#negative strand for UTR5
										
										$utr_canno = "$name:c.*" . ($rvarstart-$rcdsend) . '_*' . ($rvarend-$rcdsend) . "delins" . revcom($obs);
										$utr3{$name2} .= "$utr_canno,";
										
									} elsif ($start > $cdsend) {
										#query             ----
										#gene     <--*---*->
										#$utr5{$name2}++;		#negative strand for UTR3
										
										$utr_canno = "$name:c." . ($rcdsstart-$rvarstart) . '_' . ($rcdsstart-$rvarend) . "delins" . revcom($obs);
										$utr5{$name2} .= "$utr_canno,";
										
									} else {
										$exonic{$name2}++;
										not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '-', $curlinecount, @exonstart-$k, $nextline, $refcut, $obscut];
									}
									$foundgenic++;
									last;
								} elsif ($k < @exonstart-1 and $end < $exonstart[$k+1]) {
									$intronic{$name2}++;
									$foundgenic++;
									last;
								}
							} elsif ($end >= $exonstart[$k]) {
								$rvarstart = $txend-$end-$lenintron+1;		#all the rvarstart, rvarend are with respect to the cDNA sequence (so rvarstart corresponds to end of variants)
								
								for (my $m = $k; $m >= 0; $m--) {
									$m < $k and $lenintron += ($exonstart[$m+1]-$exonend[$m]-1);
									if ($start > $exonend[$m]) {
										#query           ----
										#gene     <--**---******---****---->
										#$rvarend = $txend-$exonstart[$m]-$lenintron+1 - ($exonstart[$m+1]-$exonend[$m]-1);		#commented out 2011feb18 due to bug (10 42244567 42244600 CACCTTTGCTTGATATGATAATATAGTGCCAAGG - hetero)
										$rvarend = $txend-$exonstart[$m+1]+1 - $lenintron + ($exonstart[$m+1]-$exonend[$m]-1);		#fixed this 2011feb18
										
										$refcut = substr ($ref, -($rvarend-$rvarstart+1));	#20170830 update obscut since the start is in intronic region
										$obscut = substr ($obs, -($rvarend-$rvarstart+1));	#20170830 update obscut since the start is in intronic region
										
										last;			#finish the circle of counting exons!!!!!
									} elsif ($start >= $exonstart[$m]) {			#the start is right located within exon
										#query        -------
										#gene     <--**---******---****---->
										$rvarend = $txend-$start-$lenintron+1;
										last;						#finish the cycle
									}
								}
								if (not defined $rvarend) {					#if rvarend is not found, then the whole tail of gene is covered
									$rvarend = $txend-$txstart-$lenintron+1;
								}
								
								#here the trick begins to differentiate UTR versus coding exonic
								if ($end < $cdsstart) {			#usually disrupt/change 3' UTR region for - strand, unless the UTR per se is also separated by introns
									#query  ----
									#gene     <--*---*->
									#$utr3{$name2}++;		#negative strand for UTR5
									
									if ($start==$end and $ref ne '-' and $obs ne '-') {	#SNV
										#print STDERR "NOTICE: end=$end cdsstart=$cdsstart lenintron=$lenintron txstart=$txstart txend=$txend rvarstart=$rvarstart rvarend=$rvarend rcdsstart=$rcdsstart rcdsend=$rcdsend\n";
										$utr_canno = "$name:c.*" . ($rvarstart-$rcdsend) . revcom($ref). ">" . revcom($obs);
									} elsif ($start==$end and $ref eq '-' and $obs ne '-') {	#insertion
										$utr_canno = "$name:c.*" . ($rvarstart-$rcdsend-1) . '_*' . ($rvarstart-$rcdsend) . "ins" . revcom($obs);
									} elsif ($start==$end and $ref ne '-' and $obs eq '-') {	#single-base deletion
										$utr_canno = "$name:c.*" . ($rvarstart-$rcdsend) . "del" . revcom($ref);
									} elsif ($ref ne '-' and $obs eq '-') {		#multi-base deletion
										$utr_canno = "$name:c.*" . ($rvarend-$rcdsend) . '_*' . ($rvarstart-$rcdsend) . "del" . revcom($ref);
									} else {
										$utr_canno = "$name:c.*" . ($rvarend-$rcdsend) . '_*' . ($rvarstart-$rcdsend) . "delins" . revcom($obs);
									}
									$utr3{$name2} .= "$utr_canno,";
									
								} elsif ($start > $cdsend or $start==$cdsend&&$ref eq '-') {	#insertions at the first base in negative strand is not a real insertion since it is not translated
									#query             ----
									#gene     <--*---*->
									#$utr5{$name2}++;		#negative strand for UTR3
									
									if ($start==$end and $ref ne '-' and $obs ne '-') {	#SNV
										$utr_canno = "$name:c." . ($rvarstart-$rcdsstart) . revcom($ref) . ">". revcom($obs);
									} elsif ($start==$end and $ref eq '-' and $obs ne '-') {	#insertion
										$utr_canno = "$name:c." . ($rvarstart-$rcdsstart-1) . '_' . ($rvarstart-$rcdsstart) . "ins". revcom($obs);
									} elsif ($start==$end and $ref ne '-' and $obs eq '-') {	#single-base deletion
										$utr_canno = "$name:c." . ($rvarstart-$rcdsstart) . "del". revcom($ref);
									} elsif ($ref ne '-' and $obs eq '-') {		#multi-base deletion
										$utr_canno = "$name:c." . ($rvarend-$rcdsstart) . '_' . ($rvarstart-$rcdsstart) . "del". revcom($ref);
									} else {
										$utr_canno = "$name:c." . ($rvarend-$rcdsstart) . '_' . ($rvarstart-$rcdsstart) . "delins" . revcom($obs);
									}
									$utr5{$name2} .= "$utr_canno,";
									
								} else {
									$exonic{$name2}++;
									not $current_ncRNA and $obs and push @{$refseqvar{$name}}, [$rcdsstart, $rvarstart, $rvarend, '-', $curlinecount, @exonstart-$k, $nextline, $refcut, $obscut];
								}
								$foundgenic++;
								last;
							}
						}
					}
				}
			}
		}
		$foundgenic or $intergenic{''}++;		#changed $name2 to '' on 20110924
		$i =~ m/000000$/ and printerr "NOTICE: Finished analyzing $i query variants\n";

	
		my (@txname, %genename);
		my (%newsplicing);
		
		#process splicing annotation (change %splicing hash to %newsplicing, where gene name is replaced by gene name plus splicing annotation)
		if (%splicing) {
			if ($end-$start+1<=$indel_splicing_threshold) {		#make sure that long indel are not considered here
				for my $tempname (keys %splicing) {
					if ($splicing_anno{$tempname}) {
						$splicing_anno{$tempname} =~ s/,$//;	#remove the trailing comma
						$tempname .= "($splicing_anno{$tempname})";
					}
					
					$tempname =~ s/#\w+#\d+//g;		#20140711 (issue:multimap) to handle one transcript multiple mapping issue
					
					$newsplicing{$tempname}++;
				}
			} else {
				#%newsplicing = %splicing;			#20140711 (issue:multimap) commented out and replaced by below
				
				for my $tempname (keys %splicing) {
					$tempname =~ s/#\w+#\d+//g;		#20140711 (issue:multimap) to handle one transcript multiple mapping issue
					$newsplicing{$tempname}++;
				}
			}
		}
		if (%utr5) {
			for my $tempname (keys %utr5) {
				$utr5{$tempname} =~ s/,$//;
				$tempname .= "($utr5{$tempname})";
				
				$tempname =~ s/#\w+#\d+//g;		#20140711 (issue:multimap) to handle one transcript multiple mapping issue
				
				$newutr5{$tempname}++;
			}
		}
		if (%utr3) {
			for my $tempname (keys %utr3) {
				$utr3{$tempname} =~ s/,$//;
				$tempname .= "($utr3{$tempname})";
				
				$tempname =~ s/#\w+#\d+//g;		#20140711 (issue:multimap) to handle one transcript multiple mapping issue
				
				$newutr3{$tempname}++;
			}
		}
		
		if ($separate) {		#separately print out each effect on one line
			if (%exonic or %splicing or %intronic or %utr5 or %utr3 or %ncrna or %upstream or %downstream) {
				#if (%ncrna) {
				#	%exonic and print OUT "ncRNA_exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#	%splicing and $end-$start+1<=$splicing_threshold and print OUT "ncRNA_splicing\t", join (',', sort keys %newsplicing), "\t", $nextline, "\n";
				#	%intronic and print OUT "ncRNA_intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
				#	
				#	for my $key (keys %ncrna) {
				#		delete $exonic{$key};
				#		delete $splicing{$key};
				#		delete $intronic{$key};
				#	}
				#}
				if (%exonic) {
					my (@coding, @noncoding);
					for my $key (keys %exonic) {
						if ($ncrna{$key}) {
							push @noncoding, $key;
						} else {
							push @coding, $key;
						}
					}
					@coding and print $VARFUN "exonic\t", join(",", sort @coding), "\t", $nextline, "\n";
					@noncoding and print $VARFUN "ncRNA_exonic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				}
				#changed per Oscar 20131109
				#if (%splicing and $end-$start+1<=$indel_splicing_threshold) {
				#	my (@coding, @noncoding);
				#	for my $key (keys %splicing) {
				#		if ($ncrna{$key}) {
				if (%splicing) {
					my (@coding, @noncoding);
					for my $key (keys %newsplicing) {						
						$key =~ m/^([^\(]+)/;
						if ($ncrna{$1}) {
							push @noncoding, $key;
						} else {
							push @coding, $key;
						}
					}
					@coding and print $VARFUN "splicing\t", join(",", sort @coding), "\t", $nextline, "\n";
					@noncoding and print $VARFUN "ncRNA_splicing\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				}
				if (%intronic) {
					my (@coding, @noncoding);
					for my $key (keys %intronic) {
						if ($ncrna{$key}) {
							push @noncoding, $key;
						} else {
							push @coding, $key;
						}
					}
					@coding and print $VARFUN "intronic\t", join(",", sort @coding), "\t", $nextline, "\n";
					@noncoding and print $VARFUN "ncRNA_intronic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				}
				
				#the following paragraph is commented out on 2011oct02
				#%exonic and print OUT "exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#%splicing and $end-$start+1<=$splicing_threshold and print OUT "splicing\t", join (',', sort keys %newsplicing), "\t", $nextline, "\n";
				#%intronic and print OUT "intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
				%utr5 and print $VARFUN "UTR5\t", join(",", sort keys %newutr5), "\t", $nextline, "\n";
				%utr3 and print $VARFUN "UTR3\t", join(",", sort keys %newutr3), "\t", $nextline, "\n";
				#if (%ncrna) {
				#	if (%exonic) {
				#		print OUT "ncRNA_exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#	}
				#	if (%splicing and $end-$start+1<=$splicing_threshold) {		#a big deletion spanning splicing site is not really a "splicing" mutation
				#		print OUT "ncRNA_splicing\t",join (",", sort keys %newsplicing), "\t", $nextline, "\n";
				#	}
				#	if (%utr5) {	#ncRNA should not have UTR5 or UTR3. If such an output exists, then there is a bug that should be reported and debugged!!!!!!!!!!!!!!!!
				#		print OUT "ncRNA_UTR5\t", join(",", sort keys %utr5), "\t", $nextline, "\n";
				#	}
				#	if (%utr3) {	#ncRNA should not have UTR5 or UTR3. If such an output exists, then there is a bug that should be reported and debugged!!!!!!!!!!!!!!!!
				#		print OUT "ncRNA_UTR3\t", join(",", sort keys %utr3), "\t", $nextline, "\n";
				#	}
				#	if (%intronic) {
				#		print OUT "ncRNA_intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
				#	}
				#}
				%upstream and print $VARFUN "upstream\t", join(",", sort keys %upstream), "(dist=$distup)\t", $nextline, "\n";
				%downstream and print $VARFUN "downstream\t", join(",", sort keys %downstream), "(dist=$distdown)\t", $nextline, "\n";

			} elsif (%intergenic) {
				$genel ||= "NONE";
				$gener ||= "NONE";
				$distl ||= "NONE";
				$distr ||= "NONE";
				print $VARFUN "intergenic\t", "$genel(dist=$distl),$gener(dist=$distr)", "\t", $nextline, "\n";
			} else {
				die "FATAL ERROR: please report bug to ANNOVAR author with your input file\n";
			}
		} else {			
			if (@precedence) {
				my $foundmatch;
				for my $i (0 .. @precedence-2) {
					$precedence[$i] eq 'exonic' and %exonic and $foundmatch++;
					$precedence[$i] eq 'splicing' and %splicing and $foundmatch++;
					$precedence[$i] eq 'intronic' and %intronic and $foundmatch++;
					$precedence[$i] eq 'utr5' and %utr5 and $foundmatch++;
					$precedence[$i] eq 'utr3' and %utr3 and $foundmatch++;
					$precedence[$i] eq 'ncrna' and %ncrna and $foundmatch++;
					$precedence[$i] eq 'upstream' and %upstream and $foundmatch++;
					$precedence[$i] eq 'downstream' and %downstream and $foundmatch++;
					$precedence[$i] eq 'intergenic' and %intergenic and $foundmatch++;
					if ($foundmatch) {
						for my $j ($i+1 .. @precedence-1) {
							$precedence[$j] eq 'exonic' and %exonic = ();
							$precedence[$j] eq 'splicing' and %splicing = ();
							$precedence[$j] eq 'intronic' and %intronic = ();
							$precedence[$j] eq 'utr5' and %utr5 = ();
							$precedence[$j] eq 'utr3' and %utr3 = ();
							$precedence[$j] eq 'ncrna' and %ncrna = ();
							$precedence[$j] eq 'upstream' and %upstream = ();
							$precedence[$j] eq 'downstream' and %downstream = ();
							$precedence[$j] eq 'intergenic' and %intergenic = ();
						}
						last;
					}
				}
			}
			
		
				
			if (%exonic) {
				my (@coding, @noncoding);
				for my $key (keys %exonic) {
					if ($ncrna{$key}) {
						push @noncoding, $key;
					} else {
						push @coding, $key;
					}
				}
				if (@coding and %splicing and $end-$start+1<=$indel_splicing_threshold) {		#a big deletion spanning splicing site is not really a "splicing" mutation
					print $VARFUN "exonic;splicing\t", join(",", sort @coding), ";", join (",", sort keys %newsplicing), "\t", $nextline, "\n";
				} elsif (@coding) {
					print $VARFUN "exonic\t", join(",", sort @coding), "\t", $nextline, "\n";
				} elsif (@noncoding and %splicing and $end-$start+1<=$indel_splicing_threshold) {	#handle situations while a variant overlaps with both noncoding exonic and splicing
					print $VARFUN "ncRNA_exonic;splicing\t", join(",", sort @noncoding), ";", join (",", sort keys %newsplicing), "\t", $nextline, "\n";
				} elsif (@noncoding) {
					print $VARFUN "ncRNA_exonic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				}
			} elsif (%splicing) {
				my (@coding, @noncoding);
				for my $key (keys %newsplicing) {
					$key =~ m/^([^\(]+)/;
					if ($ncrna{$1}) {
						push @noncoding, $key;
					} else {
						push @coding, $key;
					}
				}
				if (@coding) {
					print $VARFUN "splicing\t", join (',', sort @coding), "\t", $nextline, "\n";
				} elsif (@noncoding) {
					print $VARFUN "ncRNA_splicing\t", join (',', sort @noncoding), "\t", $nextline, "\n";
				}
			} elsif (%ncrna) {
				my (@coding, @noncoding, @utr5, @utr3);
				for my $key (keys %intronic) {
					if ($ncrna{$key}) {
						push @noncoding, $key;
					} else {
						push @coding, $key;
					}
				}
				
				for my $key (keys %utr5) {
					if ($ncrna{$key}) {
						push @utr5, $key;
					} 
				}
				for my $key (keys %utr3) {
					if ($ncrna{$key}) {
						push @utr3, $key;
					}
				}
				#print OUT "ncRNA\t", join(",", sort keys %ncrna), "\t", $nextline, "\n";
				
				#if (%exonic) {
				#	if (%splicing and $end-$start+1<=$splicing_threshold) {		#a big deletion spanning splicing site is not really a "splicing" mutation
				#		print OUT "ncRNA_exonic;ncRNA_splicing\t", join(",", sort keys %exonic), ";", join (",", sort keys %newsplicing), "\t", $nextline, "\n";
				#	} else {
				#		print OUT "ncRNA_exonic\t", join(",", sort keys %exonic), "\t", $nextline, "\n";
				#	}
				#} elsif (%splicing) {
				#	print OUT "ncRNA_splicing\t", join (',', sort keys %newsplicing), "\t", $nextline, "\n";
				#} elsif (%utr5 or %utr3) {		#ncRNA should not have UTR5 or UTR3. If such an output exists, then there is a bug that should be reported and debugged!!!!!!!!!!!!!!!!
				if (@utr5 or @utr3) {		#example (chr11:18629303C>T chr11:46879972A>G are situated in both UTR3 and intron of a non-coding antisense transcript)
					if (@utr5 and @utr3) {
						print $VARFUN "ncRNA_UTR5;ncRNA_UTR3\t", join(",", sort keys %utr5), ";", join(",", sort keys %utr3), "\t", $nextline, "\n";		#use ";" to separate UTR5 and UTR3 genes
					} elsif (@utr5) {
						print $VARFUN "ncRNA_UTR5\t", join(",", sort keys %utr5), "\t", $nextline, "\n";
					} else {
						print $VARFUN "ncRNA_UTR3\t", join(",", sort keys %utr3), "\t", $nextline, "\n";
					}
				} elsif (@noncoding) {
					print $VARFUN "ncRNA_intronic\t", join(",", sort @noncoding), "\t", $nextline, "\n";
				} else {
					die "FATAL ERROR: please report bug to ANNOVAR author with your input line <$nextline>\n";
				}
			} elsif (%utr5 or %utr3) {
				if (%utr5 and %utr3) {
					print $VARFUN "UTR5;UTR3\t", join(",", sort keys %newutr5), ";", join(",", sort keys %newutr3), "\t", $nextline, "\n";		#use ";" to separate UTR5 and UTR3 genes
				} elsif (%utr5) {
					print $VARFUN "UTR5\t", join(",", sort keys %newutr5), "\t", $nextline, "\n";
				} else {
					print $VARFUN "UTR3\t", join(",", sort keys %newutr3), "\t", $nextline, "\n";
				}
			} elsif (%intronic) {
				print $VARFUN "intronic\t", join(",", sort keys %intronic), "\t", $nextline, "\n";
			} elsif (%upstream or %downstream) {
				if (%upstream and %downstream) {
					print $VARFUN "upstream;downstream\t", join(",", sort keys %upstream), "(dist=$distup);", join(",", sort keys %downstream), "(dist=$distdown)\t", $nextline, "\n";
				} elsif (%upstream) {
					print $VARFUN "upstream\t", join(",", sort keys %upstream), "(dist=$distup)\t", $nextline, "\n";
				} else {
					print $VARFUN "downstream\t", join(",", sort keys %downstream), "(dist=$distdown)\t", $nextline, "\n";
				}
			} elsif (%intergenic) {
				$genel ||= "NONE";
				$gener ||= "NONE";
				$distl ||= "NONE";
				$distr ||= "NONE";
				print $VARFUN "intergenic\t", "$genel(dist=$distl),$gener(dist=$distr)", "\t", $nextline, "\n";
			} else {
				die "FATAL ERROR: please report bug to ANNOVAR author with your input line <$nextline>\n";
			}
		}
	}
	%refseqvar and annotateExonicVariantsThread (\%refseqvar, $geneidmap, $cdslen, $mrnalen, $EXVARFUN);
}

# Perform the exonic_variant_function annotation (that is, whether the varint is missense, nonsense, etc), based on $refseqvar produced by processNextQueryBatchByGeneThread()
sub annotateExonicVariantsThread {
	my ($refseqvar, $geneidmap, $cdslen, $mrnalen, $EXVARFUN) = @_;
	my $refseqhash;
	my $function = {};
	my %varinfo;					#variants information (same as input line)
	my %unmatch_wtnt_ref;			#count how many user input variant has wildtype nucleotide that does not match reference genome (use hash, because if we use count, we are merely counting variant-transcript combinations)
	my @unmatch_example;
	
	$refseqhash = readSeqFromFASTADB ($refseqvar);

	for my $seqid (keys %$refseqvar) {
		for my $i (0 .. @{$refseqvar->{$seqid}}-1) {
			my ($refcdsstart, $refvarstart, $refvarend, $refstrand, $index, $exonpos, $nextline, $refcut, $obscut) = @{$refseqvar->{$seqid}->[$i]};	#20170830: $obscut is included now
			my ($wtnt3, $wtnt3_after, @wtnt3, $varnt3, $wtaa, $wtaa_after, $varaa, $varpos);		#wtaa_after is the aa after the wtaa
			my ($chr, $start, $end, $ref, $obs);
			my $canno;

			my ($pre_pad, $post_pad, $wt_aa_pad, $var_aa_pad);  # Hold padded seq
			my $refcdsend = $cdslen->{$seqid} + $refcdsstart - 1;  # the end of the CDS

			my @nextline = split (/\s+/, $nextline);
			($chr, $start, $end, $ref, $obs) = @nextline[@avcolumn];
			
			if ($obscut) {
				$ref = $refcut;
				$obs = $obscut;		#20170830 update obs to obscut, for indels that covers exon/intron boundaries
			}
			
			($ref, $obs) = (uc $ref, uc $obs);
			$zerostart and $start++;
			$chr =~ s/^chr//;
			
			$varinfo{$index} = $nextline;			
			
			if (not $refseqhash->{$seqid}) {					#this refseq do not have FASTA sequence so cannot be interrogated
				$function->{$index}{unknown} = "UNKNOWN";
				next;
			}
						
			my $fs = (($refvarstart-$refcdsstart) % 3);
			my $end_fs = (($refvarend-$refcdsstart) % 3);   # Needed to complete codon following end of multibase ref seq.
			if ($refvarstart-$fs-1 > length($refseqhash->{$seqid})) {
				printerr "WARNING: Potential database annotation error seqid=$seqid, refvarstart=$refvarstart, fs=$fs, seqlength=", length($refseqhash->{$seqid}), " refcdsstart=$refcdsstart, with inputline=$nextline\n";
				next;
			}

			$wtnt3 = substr ($refseqhash->{$seqid}, $refvarstart-$fs-1, 3);
			if (length ($refseqhash->{$seqid}) >= $refvarstart-$fs+3) {	#going into UTR
				$wtnt3_after = substr ($refseqhash->{$seqid}, $refvarstart-$fs+2, 3);
			} else {
				$wtnt3_after = '';					#last amino acid in the sequence without UTR (extremely rare situation) (example: 17        53588444        53588444        -       T       414     hetero)
			}
			@wtnt3 = split (//, $wtnt3);
			if (@wtnt3 != 3 and $refvarstart-$fs-1>=0) {			#some times there are database annotation errors (example: chr17:3,141,674-3,141,683), so the last coding frame is not complete and as a result, the cDNA sequence is not complete
				$function->{$index}{unknown} = "UNKNOWN";
				next;
			}

            ##################
            # Read pre and post seq padding
            #   - The Pre and Post padding do not include the triplet bases of the affected codon. Ie,
            #   these sequences end before and start after the affected codon.  
            #   - When start==end, the other codon bases are applied in that code block, completing the codon. (Except
            #   for single nucleotide deletion, where the nt3_after is added, so must remove first 3 bases 
            #   of variant post_pad in the $seq_padding codeblock in this scenario.)
            #   - When start!=end, must apply codon bases before the refvarstart, and after refvarend.
            #   Since refvarend may now have a different fs than refvar start, be sure to use $end_fs.
            ##################
		my (@pad_gene_info, $pad_begin, $pad_end, $do_trim, $is_fs); # do_trim: set if we need to trim postpad (for variants), ie, if wtnt3_after is used. # is_fs: set if this variant annotation is a frameshift indel
		if ($seq_padding) {
			@pad_gene_info = ($refcdsstart, $refcdsend, $refseqhash->{$seqid});
			($pre_pad, $post_pad) = get_pad_seq($refvarstart, $refvarend, $cDNA_pad, $fs, $end_fs, @pad_gene_info);
			$pad_begin = $refvarstart - $fs - (length $pre_pad);
			$pad_end = $refvarend - $end_fs + 2 + (length $post_pad);
			$do_trim = 0; 
			$is_fs = 0;   
		}

			
			if ($refstrand eq '-') {					#change the observed nucleotide to the reverse strand
				$obs = revcom ($obs);
				$ref = revcom ($ref);
			}
			
			#if ($start == $end) {		#20170830 sometimes the start!=end for block substitution but the refvarstart==refvarend since a portion is in introns
			if ($refvarstart == $refvarend) {
				if ($ref eq '-') {					#insertion variant
					#the insertion coordinate system in ANNOVAR always uses "position after the current site"
					#in positive strand, this is okay
					#in negative strand, the "after current site" becomes "before current site" during transcription
					#therefore, appropriate handling is necessary to take this into account
					#for example, for a trinucleotide GCC with frameshift of 1 and insertion of CCT
					#in positive strand, it is G-CTT-CC
					#but if the transcript is in negative strand, the genomic sequence should be GC-CCT-C, and transcript is G-AGG-GC
					if ($refstrand eq '+') {
						if ($fs == 1) {
							$varnt3 = $wtnt3[0] . $wtnt3[1] . $obs . $wtnt3[2];
						} elsif ($fs == 2) {
							$varnt3 = $wtnt3[0] . $wtnt3[1] . $wtnt3[2] . $obs;
						} else {
							$varnt3 = $wtnt3[0] . $obs . $wtnt3[1] . $wtnt3[2];
						}
					} elsif ($refstrand eq '-') {
						if ($fs == 1) {
							$varnt3 = $wtnt3[0] . $obs . $wtnt3[1] . $wtnt3[2];
						} elsif ($fs == 2) {
							$varnt3 = $wtnt3[0] . $wtnt3[1] . $obs . $wtnt3[2];
						} else {
							$varnt3 = $obs . $wtnt3[0] . $wtnt3[1] . $wtnt3[2];
						}
					}
					($wtaa, $wtaa_after, $varaa, $varpos) = (translateDNA ($wtnt3, $chr), translateDNA ($wtnt3_after, $chr), translateDNA ($varnt3, $chr), int(($refvarstart-$refcdsstart)/3)+1);
					$wtaa_after and $wtaa_after eq '*' and $wtaa_after = 'X';		#wtaa_after could be undefined, if the current aa is the stop codon (X) (example: 17        53588444        53588444        -       T)

					if ($refstrand eq '+') {
						if ($obs eq substr ($refseqhash->{$seqid}, $refvarstart-1, 1)) {
							$canno = "c." . ($refvarstart-$refcdsstart+1) . "dup$obs";
						} elsif ($obs eq substr ($refseqhash->{$seqid}, $refvarstart, 1)) {
							$canno = "c." . ($refvarstart-$refcdsstart+2) . "dup$obs";
						} else {
							$canno = "c." . ($refvarstart-$refcdsstart+1) .  "_" . ($refvarstart-$refcdsstart+2) . "ins$obs";		#cDNA level annotation
						}
					} elsif ($refstrand eq '-') {
						if ($obs eq substr ($refseqhash->{$seqid}, $refvarstart-1-1, 1)) {
							$canno = "c." . ($refvarstart-$refcdsstart)  . "dup$obs";
						} elsif ($obs eq substr ($refseqhash->{$seqid}, $refvarstart-1, 1)) {
							$canno = "c." . ($refvarstart-$refcdsstart+1)  . "dup$obs";
						} else {
							$canno = "c." . ($refvarstart-$refcdsstart) .  "_" . ($refvarstart-$refcdsstart+1) . "ins$obs";			#changed 20130613 (since "after" in genomic position becomes "before" in cDNA
						}
					}
					
					if (length ($obs) % 3 == 0) {
						if ($wtaa eq '*') {			#mutation on stop codon
							if ($varaa =~ m/\*/) {
								$varaa =~ s/\*.*/X/;	#delete all aa after stop codon, but keep the aa before
								$function->{$index}{nfsins} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "delins$varaa,";		#stop codon is stil present
							} else {
								$function->{$index}{stoploss} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "delins$varaa,";	#stop codon is lost
							}
							$is_fs++;
						} else {
							if ($varaa =~ m/\*/) {
								$varaa =~ s/\*.*/X/;	#delete all aa after stop codon, but keep the aa before
								$function->{$index}{stopgain} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "delins$varaa,";
								$is_fs++;
							} else {
								$function->{$index}{nfsins} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "delins$varaa,";
							}
						}
					} else {
						if ($wtaa eq '*') {			#mutation on stop codon
							if ($varaa =~ m/\*/) {		#in reality, this cannot be differentiated from non-frameshift insertion, but we'll still call it frameshift
								$varaa =~ s/\*.*/X/;	#delete all aa after stop codon, but keep the aa before
								$function->{$index}{fsins} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "delins$varaa,";
							} else {
								$function->{$index}{stoploss} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "delins$varaa,";
							}
						} else {
							if ($varaa =~ m/\*/) {
								$varaa =~ s/\*.*/X/;	#delete all aa after stop codon, but keep the aa before
								$function->{$index}{stopgain} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "_$wtaa_after" . ($varpos+1) . "delins$varaa,";
							} else {
								$function->{$index}{fsins} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "fs,";
							}
						}
						$is_fs++;
					}
				} elsif ($obs eq '-') {					#single nucleotide deletion
					$do_trim = 3;   # Trim first 3 nt of post_pad for variant, as wtnt3_after is being added here.
					my $deletent;
					if ($fs == 1) {
						$deletent = $wtnt3[1];
						$varnt3 = $wtnt3[0].$wtnt3[2].$wtnt3_after;
					} elsif ($fs == 2) {
						$deletent = $wtnt3[2];
						$varnt3 = $wtnt3[0].$wtnt3[1].$wtnt3_after;
					} else {
						$deletent = $wtnt3[0];
						$varnt3 = $wtnt3[1].$wtnt3[2].$wtnt3_after;
					}
					($wtaa, $varaa, $varpos) = (translateDNA ($wtnt3, $chr), translateDNA ($varnt3, $chr),  int(($refvarstart-$refcdsstart)/3)+1);
					
					$canno = "c." . ($refvarstart-$refcdsstart+1) . "del$deletent";
					if ($wtaa eq '*') {				#mutation on stop codon
						if ($varaa =~ m/\*/) {			#stop codon is still stop codon
							$function->{$index}{nfsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "X,";	#changed fsdel to nfsdel on 2011feb19
						} else {				#stop codon is lost
							$function->{$index}{stoploss} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos" . "$varaa,";
						}
					} else {
						if ($varaa =~ m/\*/) {			#new stop codon created
							$function->{$index}{stopgain} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "X,";
						} else {
							$function->{$index}{fsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos" . "fs,";
						}
					}
					$is_fs++;
				} elsif (length ($obs) > 1) {				#block substitution (since start==end, this changed from 1nt to several nt)
		                    	if ($fs == 1) {
		                        	$varnt3 = $wtnt3[0] . $obs . $wtnt3[2];
		                    	}
		                    	elsif ($fs == 2) {
		                        	$varnt3 = $wtnt3[0] . $wtnt3[1] . $obs;
		                    	}
		                    	else {
		                        	$varnt3 = $obs . $wtnt3[1] . $wtnt3[2];
		                    	}
		                    	
		                    	
		                    	#the input "chr14   62162525        62162525        G       GT" needs to be handled correctly c.3delinsGT
					#$canno = "c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "delins$obs";
					$canno = "c." . ($refvarstart-$refcdsstart+1) . "delins$obs";	
					
					if (($refvarend-$refvarstart+1-length($obs)) % 3 == 0) {
						#$function->{$index}{nfssub} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "delins$obs,";
						$function->{$index}{nfssub} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:" . $canno . ",";
					} else {
						#$function->{$index}{fssub} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "delins$obs,";
						$function->{$index}{fssub} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:" . $canno . ",";
						$is_fs++;
					}
				} else {						#single nucleotide substitution variant
					if ($fs == 1) {
						$varnt3 = $wtnt3[0] . $obs . $wtnt3[2];
						$canno = "c.$wtnt3[1]" . ($refvarstart-$refcdsstart+1) . $obs;
						$hgvs and $canno = "c." . ($refvarstart-$refcdsstart+1) . $wtnt3[1] . ">" . $obs;
						if ($ref and $wtnt3[1] ne $ref) {
							$unmatch_wtnt_ref{$index}++;
							@unmatch_example or push @unmatch_example, $nextline;
							$verbose and printerr "WARNING: ALLELE MISMATCH: strand=$refstrand user-specified-allele=$ref annovar-inferred-allele=$wtnt3[1] inputline=<$nextline>\n";
						}
					} elsif ($fs == 2) {
						$varnt3 = $wtnt3[0] . $wtnt3[1]. $obs;
						$canno = "c.$wtnt3[2]" . ($refvarstart-$refcdsstart+1) . $obs;
						$hgvs and $canno = "c." . ($refvarstart-$refcdsstart+1) . $wtnt3[2] . ">" . $obs;
						if ($ref and $wtnt3[2] ne $ref) {
							$unmatch_wtnt_ref{$index}++;
							@unmatch_example or push @unmatch_example, $nextline;
							$verbose and printerr "WARNING: ALLELE MISMATCH: strand=$refstrand user-specified-allele=$ref annovar-inferred-allele=$wtnt3[2] inputline=<$nextline>\n";
						}
					} else {
						$varnt3 = $obs . $wtnt3[1] . $wtnt3[2];
						$canno = "c.$wtnt3[0]" . ($refvarstart-$refcdsstart+1) . $obs;
						$hgvs and $canno = "c." . ($refvarstart-$refcdsstart+1) . $wtnt3[0] . ">" . $obs;
						if ($ref and $wtnt3[0] ne $ref) {
							$unmatch_wtnt_ref{$index}++;
							@unmatch_example or push @unmatch_example, $nextline;
							$verbose and printerr "WARNING: ALLELE MISMATCH: strand=$refstrand user-specified-allele=$ref annovar-inferred-allele=$wtnt3[0] inputline=<$nextline>\n";
						}
					}
					($wtaa, $varaa, $varpos) = (translateDNA ($wtnt3, $chr), translateDNA ($varnt3, $chr), int(($refvarstart-$refcdsstart)/3)+1);
					#print STDERR "seqid=$seqid wtaa=$wtaa varaa=$varaa varpos=$varpos nssnvfun=", $function->{$index}{nssnv}, "\n";
					if ($wtaa eq $varaa) {
						$wtaa eq '*' and ($wtaa, $varaa) = qw/X X/;		#change * to X in the output
						$function->{$index}{ssnv} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos$varaa,";
						#print STDERR "Found a synonymous NSSNV!!!!!\n";
					} elsif ($varaa eq '*') {
						$function->{$index}{stopgain} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa${varpos}X,";
						$is_fs++;  # mark this as $is_fs to trim at stopgain
					} elsif ($wtaa eq '*') {
						$function->{$index}{stoploss} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.X$varpos$varaa,";
						$is_fs++;  # mark this as $is_fs to extend beyond stoploss (although will not go past cdsstop)
					} else {
						if ($aamatrix) {
							$function->{$index}{nssnv} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos$varaa:AAMatrix=".$aamatrix->{$wtaa.$varaa}.",";
						} else {
							$function->{$index}{nssnv} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa$varpos$varaa,";
						}
						#print STDERR "Found a NSSNV!!!!!\n";
					}

				}

		                # padded
		                if ($seq_padding) {
		                    my $var_post_pad = $post_pad;
		                    if (! $wtnt3 ) {
		                        die "wtn3 not defined: $seqid\n";
		                    }
		                    if (! $varnt3 ) {
		                        die "varnt3 not defined: $seqid\n";
		                    }
		                    if ($do_trim) {
		                        if ($do_trim < (length $post_pad)) {
		                            $var_post_pad = substr($post_pad, $do_trim);
		                        }
		                        else {
		                            $var_post_pad = "";
		                            warn "Line $seqid - \$post_pad shorter than $do_trim\n";
		                        }
		                    }
		                    $canno ||= "";
		                    $wtaa ||= "";
		                    $varaa ||= "";
		                    $varpos ||= "";
		
		                    my $wt_pad = $pre_pad . $wtnt3 . $post_pad;
		                    my $var_pad = $pre_pad . $varnt3 . $var_post_pad;
		
		                    # Extend or trim aa sequence to first stop codon
		                    if ($is_fs) {
		                        if ($is_fs > 1) {
		                            warn "Line $seqid: \$is_fs set multiple times, which means multiple conditions satified."
		                                . " This should not happen!";
		                        }
		                        my $temp_post = ( get_pad_seq($refvarstart, $pad_end, -1, 0, 2, @pad_gene_info) )[1];
		                        my $temp_aa = translateDNA($var_pad . $temp_post, $chr);
		                        my $stop_pos = undef;
		                        if ($temp_aa && $temp_aa =~ /\*/g) {
		                            # get the position (in translation of $temp_post) of stop
		                            $stop_pos = pos($temp_aa) - (length translateDNA($var_pad, $chr));
		                        }
		                        #print STDERR "stop_pos: $stop_pos temp_aa: $temp_aa\n";
		                        if (defined $stop_pos && $stop_pos >= 0) {
		                            $var_pad .= ( get_pad_seq($refvarstart, $pad_end, $stop_pos*3, 0, 2, @pad_gene_info) )[1];
		                        }
		                        elsif (defined $stop_pos && $stop_pos < 0) {
		                            $var_pad = substr($var_pad, 0, $stop_pos*3);
		                        }
		                        elsif ($temp_aa) {
		                            $var_pad .= $temp_post;
		                        }
		                    }
		
		                    ($wt_aa_pad, $var_aa_pad) = ( translateDNA( $wt_pad, $chr ),
		                                                  translateDNA( $var_pad, $chr)
		                                                );
		                    print $pad_fh "$geneidmap->{$seqid}\t$seqid\t$refstrand\t$canno\t$wt_pad\t$var_pad\tp.$wtaa$varpos$varaa\t$wt_aa_pad\t$var_aa_pad\t$varinfo{$index}\n";
		                }

			} elsif ($obs eq '-') {				#deletion variant involving several nucleotides
				($wtaa, $varpos) = (translateDNA ($wtnt3, $chr), int(($refvarstart-$refcdsstart)/3)+1);		#wildtype amino acid, position of amino acid
				my ($varposend, $canno, $panno);		#the position of the last amino acid in the deletion
				if ($refvarstart<=$refcdsstart and $firstcodondel) {	#since the first amino acid is deleted, the whole gene is considered deleted
					$function->{$index}{fsdel} .= "$geneidmap->{$seqid}:$seqid:wholegene,";	#it is exonic variant, so the varend has to hit the first exon
				} elsif ($refvarstart<=$refcdsstart and not $firstcodondel) {	#the first amino acid is deleted but we still report the predicted functional consequence
					if ($refvarend >= $cdslen->{$seqid}+$refcdsstart) {	#3' portion of the gene is deleted
						$varposend = int ($cdslen->{$seqid}/3);
						$canno = "c." . ($refvarstart-$refcdsstart+1) . "_" . ($cdslen->{$seqid}+$refcdsstart-1) . "del";
						$is_fs++;
						$panno = '';
					} else {
						$varposend = int (($refvarend-$refcdsstart)/3) + 1;
						if ($refvarend-$refcdsstart > 0) {
							$canno = "c." . 1 . "_" . ($refvarend-$refcdsstart+1) . "del";	#added 20120618
							$panno = ":p.$wtaa${varpos}fs";
						} else {
							$canno = "c." . 1 . "_" . $cdslen->{$seqid} . "del";			#added 20191010 (the maximum cdot is the length of the transcript
							$panno = '';
						}
						($refvarend-$refvarstart+1) % 3 == 0 or $is_fs++;
					}
					#$function->{$index}{fsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.${varpos}_${varposend}del,";	#per user comments, this is changed to 'fs' below 20140711 (issue:blockfs)
					$function->{$index}{fsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno$panno,";
				} elsif ($refvarend >= $cdslen->{$seqid}+$refcdsstart) {	#3' portion of the gene is deleted
					$varposend = int ($cdslen->{$seqid}/3);		#cdslen should be multiples of 3, but just in case of database mis-annotation
					#$canno = "c." . ($refvarstart-$refcdsstart+1) . "_" . ($cdslen->{$seqid}+$refcdsstart-1) . "del";	#commented 20170601 due to error on 8       8887543 8887545 AAC     - (hg19)
					#this is a dilemma: the deletion is more than the coding sequence ends, so I have to pretend that the coding sequence is more than what it shoudl be in the c. notation
					$canno = "c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "del";
					#$function->{$index}{fsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.${varpos}_${varposend}del,";	#per user comments, this is changed to 'fs' below 20140711 (issue:blockfs)
					$function->{$index}{fsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa${varpos}fs,";
					$is_fs++;
				} elsif (($refvarend-$refvarstart+1) % 3 == 0) {	#regular nonframeshift deletion
					$varposend = int (($refvarend-$refcdsstart)/3) + 1;
					$canno = "c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "del";
					$function->{$index}{nfsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.${varpos}_${varposend}del,";
				} else {	#regular frameshift deletion
					$varposend = int (($refvarend-$refcdsstart)/3) + 1;
					$canno = "c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "del";
					#$function->{$index}{fsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.${varpos}_${varposend}del,";	#per user comments, this is changed to 'fs' below 20140711 (issue:blockfs)
					$function->{$index}{fsdel} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno:p.$wtaa${varpos}fs,";
					$is_fs++;
				}

		                # padded
		                if ($seq_padding) {
		                    my $real_ref = ($ref eq '-' ? '' : substr($refseqhash->{$seqid}, $refvarstart-1, $refvarend-$refvarstart+1));
		                    my $real_obs = ($obs eq '-' ? '' : $obs);
		                    $pre_pad .= substr($refseqhash->{$seqid}, $refvarstart-$fs-1, $fs);
		                    $post_pad = substr($refseqhash->{$seqid}, $refvarend, 2-$end_fs) . $post_pad;
		
		                    $canno ||= "";
		                    
		                    my $wt_pad = $pre_pad . $real_ref . $post_pad;
		                    my $var_pad = $pre_pad . $real_obs . $post_pad;
		                    
		                    # Extend or trim aa sequence to first stop codon
		                    if ($is_fs) {
		                        if ($is_fs > 1) {
		                            warn "Line $seqid: \$is_fs set multiple times, which means multiple conditions satified."
		                                . " This should not happen!";
		                        }
		                        my $temp_post = ( get_pad_seq($refvarstart, $pad_end, -1, 0, 2, @pad_gene_info) )[1];
		                        my $temp_aa = translateDNA($var_pad . $temp_post, $chr);
		                        my $stop_pos = undef;
		                        if ($temp_aa && $temp_aa =~ /\*/g) {
		                            # get the position (in translation of $temp_post) of stop
		                            $stop_pos = pos($temp_aa) - (length translateDNA($var_pad, $chr));
		                        }
		                        #print STDERR "stop_pos: $stop_pos temp_aa: $temp_aa\n";
		                        if (defined $stop_pos && $stop_pos >= 0) {
		                            $var_pad .= ( get_pad_seq($refvarstart, $pad_end, $stop_pos*3, 0, 2, @pad_gene_info) )[1];
		                        }
		                        elsif (defined $stop_pos && $stop_pos < 0) {
		                            $var_pad = substr($var_pad, 0, $stop_pos*3);
		                        }
		                        elsif ($temp_aa) {
		                            $var_pad .= $temp_post;
		                        }
		                    }
		                    
		                    ($wt_aa_pad, $var_aa_pad) = ( translateDNA( $wt_pad, $chr ),
		                                                  translateDNA( $var_pad, $chr)
		                                                );
		
					print $pad_fh "$geneidmap->{$seqid}\t$seqid\t$refstrand\t$canno\t$wt_pad\t$var_pad\tp.${varpos}_${varposend}del\t$wt_aa_pad\t$var_aa_pad\t$varinfo{$index}\n";
				}

			} else {							#block substitution event


				#20191010: the block below is changed (1) add delins in the string (2) change the start number to be 1 less if cdot starts before 1
				if ($refvarstart-$refcdsstart < 0) {	#the start is before the cds, because cdot change from -1 to 1 without 0, we have to treat it differently here
					$canno = "c." . ($refvarstart-$refcdsstart) . "_" . ($refvarend-$refcdsstart+1) . "delins$obs";	#for example, 37_16-88923279-CCGCCATG-CA__NM_000512 should be c.-1_6T rather than c.0_6T
				} else {
					$canno = "c." . ($refvarstart-$refcdsstart+1) . "_" . ($refvarend-$refcdsstart+1) . "delins$obs";
				}            

				if (($refvarend-$refvarstart+1-length($obs)) % 3 == 0) {
					$function->{$index}{nfssub} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno,";
				} else {
					$function->{$index}{fssub} .= "$geneidmap->{$seqid}:$seqid:exon$exonpos:$canno,";
					$is_fs++;
				}

		                # padded
		                if ($seq_padding) {
		                    my $real_ref = ($ref eq '-' ? '' : substr($refseqhash->{$seqid}, $refvarstart-1, $refvarend-$refvarstart+1));
		                    my $real_obs = ($obs eq '-' ? '' : $obs);
		                    $pre_pad .= substr($refseqhash->{$seqid}, $refvarstart-$fs-1, $fs);
		                    $post_pad = substr($refseqhash->{$seqid}, $refvarend, 2-$end_fs) . $post_pad;
		                    
		                    $canno ||= "";
		                    
		                    my $wt_pad = $pre_pad . $real_ref . $post_pad;
		                    my $var_pad = $pre_pad . $real_obs . $post_pad;
		
		                    # Extend or trim aa sequence to first stop codon
		                    if ($is_fs) {
		                        if ($is_fs > 1) {
		                            warn "Line $seqid: \$is_fs set multiple times, which means multiple conditions satified."
		                                . " This should not happen!";
		                        }
		                        my $temp_post = ( get_pad_seq($refvarstart, $pad_end, -1, 0, 2, @pad_gene_info) )[1];
		                        my $temp_aa = translateDNA($var_pad . $temp_post, $chr);
		                        my $stop_pos = undef;
		                        if ($temp_aa && $temp_aa =~ /\*/g) {
		                            # get the position (in translation of $temp_post) of stop
		                            $stop_pos = pos($temp_aa) - (length translateDNA($var_pad, $chr));
		                        }
		                        #print STDERR "stop_pos: $stop_pos temp_aa: $temp_aa\n";
		                        if (defined $stop_pos && $stop_pos >= 0) {
		                            $var_pad .= ( get_pad_seq($refvarstart, $pad_end, $stop_pos*3, 0, 2, @pad_gene_info) )[1];
		                        }
		                        elsif (defined $stop_pos && $stop_pos < 0) {
		                            $var_pad = substr($var_pad, 0, $stop_pos*3);
		                        }
		                        elsif ($temp_aa) {
		                            $var_pad .= $temp_post;
		                        }
		                    }
		
		                    ($wt_aa_pad, $var_aa_pad) = ( translateDNA( $wt_pad, $chr ),
		                                                  translateDNA( $var_pad, $chr)
		                                                );
		
		                    print $pad_fh "$geneidmap->{$seqid}\t$seqid\t$refstrand\t$canno\t$wt_pad\t$var_pad\t\t$wt_aa_pad\t$var_aa_pad\t$varinfo{$index}\n";
		                }

			}

		}
	}
	
	for my $index (sort {$a<=>$b} keys %$function) {
		#my $lineindex = $index + ($batchcount-1)*$batchsize;		#this may be the second batch if input is over 5 million lines
		my $lineindex = $index;
		if ($separate) {			#print out each type of exonic mutations separately (one effect in one line), rather than printing out only the most important function
			my $sortout;
			if ($function->{$index}{fsins}) {
				$sortout = cleanExonicAnnotation ($function->{$index}{fsins});
				print $EXVARFUN "line$lineindex\t", "frameshift insertion\t$sortout\t", $varinfo{$index}, "\n";
			}
			if ($function->{$index}{fsdel}) {
				$sortout = cleanExonicAnnotation ($function->{$index}{fsdel});
				print $EXVARFUN "line$lineindex\t", "frameshift deletion\t$sortout\t", $varinfo{$index}, "\n";
			}
			if ($function->{$index}{fssub}) {
				$sortout = cleanExonicAnnotation ($function->{$index}{fssub});
				print $EXVARFUN "line$lineindex\t", "frameshift substitution\t$sortout\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{stopgain}) {
				$sortout = cleanExonicAnnotation ($function->{$index}{stopgain});
				print $EXVARFUN "line$lineindex\t", "stopgain\t$sortout\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{stoploss}) {
				$sortout = cleanExonicAnnotation ($function->{$index}{stoploss});
				print $EXVARFUN "line$lineindex\t", "stoploss\t$sortout\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{nfsins}) {
				$sortout = cleanExonicAnnotation ($function->{$index}{nfsins});
				print $EXVARFUN "line$lineindex\t", "nonframeshift insertion\t$sortout\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{nfsdel}) {
				$sortout = cleanExonicAnnotation ($function->{$index}{nfsdel});
				print $EXVARFUN "line$lineindex\t", "nonframeshift deletion\t$sortout\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{nfssub}) {
				$sortout = cleanExonicAnnotation ($function->{$index}{nfssub});
				print $EXVARFUN "line$lineindex\t", "nonframeshift substitution\t$sortout\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{nssnv}) {
				$sortout = cleanExonicAnnotation ($function->{$index}{nssnv});
				print $EXVARFUN "line$lineindex\t", "nonsynonymous SNV\t$sortout\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{ssnv}) {
				$sortout = cleanExonicAnnotation ($function->{$index}{ssnv});
				print $EXVARFUN "line$lineindex\t", "synonymous SNV\t$sortout\t", $varinfo{$index}, "\n";
			} 
			if ($function->{$index}{unknown}) {
				$sortout = cleanExonicAnnotation ($function->{$index}{unknown});
				print $EXVARFUN "line$lineindex\t", "unknown\t$sortout\t", $varinfo{$index}, "\n";
			}
		} else {				#print out only the most important functional changes (for example, chr3:9931279-9931279 G->A can be both non-synonymous and synonymous mutations based on UCSC gene model)
			print $EXVARFUN "line$lineindex\t";
			my $sortout;
			if ($sortout = $function->{$index}{fsins}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				$sortout = cleanExonicAnnotation ($sortout);
				print $EXVARFUN "frameshift insertion\t$sortout\t";
			} elsif ($sortout = $function->{$index}{fsdel}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				$sortout = cleanExonicAnnotation ($sortout);
				print $EXVARFUN "frameshift deletion\t$sortout\t";
			} elsif ($sortout = $function->{$index}{fssub}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				$sortout = cleanExonicAnnotation ($sortout);
				print $EXVARFUN "frameshift substitution\t$sortout\t";
			} elsif ($sortout = $function->{$index}{stopgain}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				$sortout = cleanExonicAnnotation ($sortout);
				print $EXVARFUN "stopgain\t$sortout\t";
			} elsif ($sortout = $function->{$index}{stoploss}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				$sortout = cleanExonicAnnotation ($sortout);
				print $EXVARFUN "stoploss\t$sortout\t";
			} elsif ($sortout = $function->{$index}{nfsins}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				$sortout = cleanExonicAnnotation ($sortout);
				print $EXVARFUN "nonframeshift insertion\t$sortout\t";
			} elsif ($sortout = $function->{$index}{nfsdel}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				$sortout = cleanExonicAnnotation ($sortout);
				print $EXVARFUN "nonframeshift deletion\t$sortout\t";
			} elsif ($sortout = $function->{$index}{nfssub}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				$sortout = cleanExonicAnnotation ($sortout);
				print $EXVARFUN "nonframeshift substitution\t$sortout\t";
			} elsif ($sortout = $function->{$index}{nssnv}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				$sortout = cleanExonicAnnotation ($sortout);
				print $EXVARFUN "nonsynonymous SNV\t$sortout\t";
			} elsif ($sortout = $function->{$index}{ssnv}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				$sortout = cleanExonicAnnotation ($sortout);
				print $EXVARFUN "synonymous SNV\t$sortout\t";
			} elsif ($sortout = $function->{$index}{unknown}) {
				$exonsort and $sortout = sortExonicAnnotation ($sortout);
				$sortout = cleanExonicAnnotation ($sortout);
				print $EXVARFUN "unknown\t$sortout\t";
			}
			print $EXVARFUN $varinfo{$index}, "\n";
		}
	}
	
	if (%unmatch_wtnt_ref) {
		printerr "\n";
		printerr "---------------------------------------------------------------------------------------\n";
		printerr "WARNING: ${\(scalar keys %unmatch_wtnt_ref)} exonic SNPs have WRONG reference alleles specified in your input file!\n";
		printerr "WARNING: An example input line is <$unmatch_example[0]>\n";
		printerr "WARNING: ANNOVAR can still annotate exonic_variant_function for the mutation correctly!\n";
		printerr "WARNING: you may have used wrong -buildver, or specified incorrect reference allele, or used outdated mRNA FASTA file!\n";
		printerr "---------------------------------------------------------------------------------------\n\n";
	}
}

# Clean the exonic annotation string when it has extra comment that should not be printed out
sub cleanExonicAnnotation {
	my ($anno) = @_;
	$anno =~ s/#[\w\.\-]+#\d+//g;	#to handle situations where the chromosome name has dot or - sign
	return $anno;
}

# Sort the exonc annotation string so that transcripts are printed out in a way that lower exon number is printed out first
sub sortExonicAnnotation {
	my ($anno) = @_;
	my @anno1 = split (/,/, $anno);
	my @anno2;
	#example: ATG16L1:NM_198890:exon5:c.A409G:p.T137A,ATG16L1:NM_017974:exon8:c.A841G:p.T281A
	#example: DNAJA4:NM_001130183:wholegene,
	#example: UNKNOWN
	for my $i (0 .. @anno1-1) {
		my @temp = split (/:/, $anno1[$i]);
		if (@temp<3) {
			$temp[0] eq 'UNKNOWN' or print STDERR "Potential Bug: please report the annotation string (<$anno>) to ANNOVAR developers\n";
		} elsif (@temp==3) {
			$temp[2] eq 'wholegene' or print STDERR "Potential Bug: please report the annotation string (<$anno>) to ANNOVAR developers\n";
		} elsif ($temp[2] =~ s/^exon//) {		#some are wholegene, not exon
			push @anno2, [$anno1[$i], @temp];
		} elsif ($temp[2]=~ s/wholegene.*/0/) {
			push @anno2, [$anno1[$i], @temp];
		}
	}
	if (@anno2 and @anno1==@anno2) {
		@anno2 = sort {$a->[3] <=> $b->[3] or $a->[2] cmp $b->[2]} @anno2;		#first sort by exon number, then by transcript name
		my @anno3 = map {$_->[0]} @anno2;
		return join (',', @anno3);
	} else {
		return $anno;
	}
}

# Control filter-based annotation by calling filterNextBatchThread, with support for multi-threading if $start_line, $end_line and $cur_thread is specified
sub filterQueryThread {
	my ($filtered_file, $dropped_file, $invalid_file, $start_line, $end_line, $cur_thread) = @_;
	
	my ($FIL, $DROPPED, $INVALID);
	open ($FIL, ">$filtered_file") or die "Error: cannot write to output file $outfile.$filtered_file: $!\n"; 
	open ($DROPPED, ">$dropped_file") or die "Error: cannot write to output file $dropped_file: $!\n";
	open ($INVALID, ">$invalid_file") or die "Error: cannot write to output file $invalid_file: $!\n";
			
	open (QUERY, $queryfile) or die "Error: cannot read from query file $queryfile: $!\n";
	
	my (%variant, $filedone, $batchdone);
	my ($linecount, $batchlinecount, $invalid, $invalidcount) = (0, 0);
	my ($chr, $start, $end, $ref, $obs, $info);
	
	
	if ($thread and $start_line) {
		for (1 .. $start_line-1) {
			$_ = <QUERY>;
		}
	}
	while (1) {
		$_ = <QUERY>;
		
		if ($thread) {
			if ($. > $end_line) {
				$_ = undef;		#do not read any line any more
			}
		}
		
		if (not defined $_) {
			$filedone++;
		} else {
			s/[\r\n]+$//;
			
			if (m/^#/ and $comment) {				#comment line start with #, do not include this is $linecount
				print $FIL "$_\n";
				print $DROPPED "#comment\t#comment\t$_\n";
				next;
			}
			
			$linecount++;
			$batchlinecount++;
			if ($batchlinecount == $batchsize) {
				$batchdone++;
			}
			
			if ($memfree or $memtotal) {		#if these arguments are specified
				if ($linecount =~ m/00000$/) {						#about 40Mb memory per 10k lines for a typical input dataset
					my ($availmem, $allmem) = currentAvailMemory();
					$verbose and printerr "NOTICE: Current available system memory is $availmem kb (this program uses $allmem bytes memory), after reading $linecount query\n";
					if ($availmem and $availmem <= $memfree+50_000) {		#some subsequent steps may take ~50Mb memory, so here we try to allocate some more memory
						$batchdone++;
					}
					if ($memtotal and $allmem >= $memtotal-50_000) {	#when --memtotal is specified, ensure that program use less memory
						$batchdone++;
					}
				}
			}
	
			($invalid, $chr, $start, $end, $ref, $obs) = detectInvalidInput ($_);
						
			if ($invalid) {
				print $INVALID $_, "\n";	#invalid record found
				$invalidcount++;
				next;
			}
			
			if ($start == $end and $ref eq '-') {	#insertion
				$obs = "0$obs";
			} elsif ($obs eq '-') {			#deletion
				$obs = $end-$start+1;
			} elsif ($end>$start or $start==$end and length($obs)>1) {	#block substitution	#fixed the bug here 2011feb19
				$obs = ($end-$start+1) . $obs;
			}
			
			if (exists $variant{$chr, $start, $obs}) {
				$variant{$chr, $start, $obs} .= "\n$_";
			} else {
				$variant{$chr, $start, $obs} = "$ref\n$_";	#store the information for the variant in newline-separated strings, with the first field being ref allele
			}
		}
		
		if ($filedone or $batchdone) {
			%variant and printerr "NOTICE: Processing next batch with ${\(scalar keys %variant)} unique variants in $batchlinecount input lines\n";
			%variant and filterNextBatchThread (\%variant, $FIL, $DROPPED, $INVALID);		#filter the next batch of variants (print output to FIL, DROPPED), then clear %variant hash for next batch
			%variant = ();
			$batchlinecount = 0;				#reset the line count for this batch
			$batchdone = 0;
		}
		if ($filedone) {
			close ($FIL);
			close ($DROPPED);
			close ($INVALID);
			last;
		}
	}
}

# Perform filter-based annotation on a batch of variants
sub filterNextBatchThread {
	my ($variant, $FIL, $DROPPED, $INVALID) = @_;
	my $dbfile;
	
	if ($dbtype1 eq 'generic') {
		$dbfile = File::Spec->catfile ($dbloc, $genericdbfile);
	} elsif ($dbtype1 eq 'vcf') {
		$dbfile = File::Spec->catfile ($dbloc, $vcfdbfile);
	} else {
		$dbfile = File::Spec->catfile ($dbloc, "${buildver}_$dbtype1.txt");
	}
	
	-f $dbfile or die "Error: the database file $dbfile is not present\n";
	
	my (@record, $chr, $start, $end, $ref, $obs, $score, @otherinfo, $qual, $fil, $info);		#@otherinfo: other information about the variant in addition to score (by default, only score for the variant will be printed out)
	my ($rsid, $strand, $ucscallele, $twoallele, $class, $af, $attribute);
	my $count_invalid_dbline;

	my ($BIN, $DBSIZE, $NUMBIN) = (0, 0);
	my %index = ();			#holds all bin index, used by DBM only and not by IDX (20140725)
	my $bb = {};			#a subset of %index, which corresponds to the input variants
	my $flag_idx_search = 0;	#indicate if index-based search algorithm is used (faster speed for a small number of input variants)
	
	if ($dbm) {				#dbm is no longer used in ANNOVAR but we provide it here for backward compatibility, read https://en.wikipedia.org/wiki/DBM_(computing)
		if (-f "$dbfile.dir" and -f "$dbfile.pag") {
			if (dbmopen (%index, $dbfile, 0)) {		#read only mode
				$BIN = $index{BIN} || undef;
				$DBSIZE = $index{DBSIZE} || undef;
				$NUMBIN = $index{NUMBIN} || 'unknown';
				if (defined $BIN and $BIN =~ m/^\d+$/ and $BIN > 0 and defined $DBSIZE and $DBSIZE == -s $dbfile) {
					printerr "NOTICE: DBM information: BIN=$BIN DBSIZE=$DBSIZE\n";
					$flag_idx_search++;	#okay to proceed with DBM search
				
					my $count_search = 0;
					foreach my $k ( keys %$variant ) {
						my ($chrom, $pos) = split ($;, $k);
						
						my $bin = $pos - ($pos % $BIN);
						defined $index{"$chrom\t$bin"} or next;		#this BIN is not in index database, skipping it
				
						$bb->{"$chrom\t$bin"} or $count_search++;	#if this bin has not been added, add it to candidate
						$bb->{"$chrom\t$bin"} = $index{"$chrom\t$bin"};
						
					}
					printerr "NOTICE: Total of ", scalar (keys %$bb), " bins will be scanned\n";
					printerr "NOTICE: DBM index loaded. Total number of bins is $NUMBIN and the number of bins to be scanned is $count_search\n";
				} else {
					printerr "WARNING: DBM index $dbfile.dir and $dbfile.pag contain invalid information. Switching to conventional filter-based search.\n";
				}
				
			} else {
				printerr "WARNING: Failed to open DBM index $dbfile.dir and $dbfile.pag. Switching to conventional filter-based search.\n";
			}
			
			
		}
	}

	if (not $flag_idx_search and -f "$dbfile.idx") {		#DBM failed or not used, and idx file exists
		
		if (open(IDX, "$dbfile.idx")) {
			my $line = <IDX>;
			if (not $line =~ m/BIN\t(\d+)\t(\d+)/) {
				printerr "WARNING: Malformed database index file $dbfile.idx.\n";
			} elsif ($2 != -s $dbfile) {			#file version is different, do not use index file in this case
				printerr "WARNING: Your index file $dbfile.idx is out of date and will not be used. ANNOVAR can still generate correct results without index file.\n";
			} else {
				($BIN, $DBSIZE) = ($1, $2);

				foreach my $k ( keys %$variant ) {
					my ($chrom, $pos) = split ($;, $k);
					
					my $bin = $pos - ($pos % $BIN);
					
					$bb->{"$chrom\t$bin"} = 0;	#flag this bin to be searched later
				}				
				
				my ($count_total, $count_search) = qw/0 0/;
				while ( $line = <IDX> ) {
					$line =~ s/[\r\n]+$//;
					my ( $chrom, $bin, $offset0, $offset1 ) = split (/\t/, $line);
					$chrom =~ s/^chr//;		#delete the chr in snp135, etc 
					defined $offset1 or next;	#invalid input line in the index file
					
					if (defined $bb->{"$chrom\t$bin"}) {
						$bb->{"$chrom\t$bin"} or $count_search++;
						$bb->{"$chrom\t$bin"} = "$offset0,$offset1";
					}
					$count_total++;
				}
				if ($count_total and $count_search/$count_total < $indexfilter_threshold) {
					$flag_idx_search++;
					printerr "NOTICE: Database index loaded. Total number of bins is $count_total and the number of bins to be scanned is $count_search\n";
				}
			}
			close (IDX);
		} else {
			printerr "WARNING: Failed to read from filter database index $dbfile.idx: $!. Switching to text scanning-based filter operation\n";
		}
	}

	if (not $flag_idx_search) {
		$bb = {1, join (',', 0, -s "$dbfile")};
		%index = (1, join (',', 0, -s "$dbfile"));
	}

	open (DB, $dbfile) or die "Error: cannot read from input database file $dbfile: $!\n";
	printerr "NOTICE: Scanning filter database $dbfile...";
	
	foreach my $b (sort keys %$bb) {
		
		$bb->{$b} or next;					#DB does not have this bin to search against (value is zero)
		
		my ($chunk_min, $chunk_max) = split (/,/, $bb->{$b});
		
		defined $chunk_max or die "Error: b=$b index=$index{$b} bb=$bb->{$b}\n";

		seek(DB, $chunk_min, 0);				#place file pointer to the chunk_min
		my $chunk_here = $chunk_min;

		while (<DB>) {
			my $line_length = length($_);			#calculate line length of the DB

			$chunk_here += $line_length;
	
			my (@obs2, @score2, @start2, @ref2);				#for 1000G2010 data set in VCF format, some tri-allelic SNPs are present; in the future, some quad-allelic SNPs may be also present in VCF files
			s/[\r\n]+$//;
			m/\S/ or next;					#skip empty lines in the database file (sometimes this occurs)
			m/^#/ and next;					#skip the comment line
			if ($dbtype eq 'avsift') {
				@record = split (/\t/, $_);
				@record == 8 or die "Error: invalid record found in DB file $dbfile (8 tab-delimited fields expected): <$_>\n";
				($chr, $start, $end, $ref, $obs, $score, @otherinfo) = @record;
				$chr =~ s/^chr//i;		#added 20130513
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				
				if (defined $score and defined $sift_threshold) {
					if ($reverse) {
						$score > $sift_threshold and next;
					} else {
						$score < $sift_threshold and next;
					}
				}
			} elsif ($dbtype =~ m/^ljb\w+_/) {
				@record = split (/\t/, $_);
				@record >= 5 or die "Error: invalid record found in DB file $dbfile (at least 5 tab-delimited fields expected): <$_>\n";
				($chr, $start, $end, $ref, $obs, $score, @otherinfo) = @record;
				$chr =~ s/^chr//i;		#added 20130513
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				if (defined $score and defined $score_threshold and $score=~/\d{1,}\.{0,1}\d{0,}/) {
					if ($reverse) {
						$score > $score_threshold and next;
					} else {
						$score < $score_threshold and next;
					}
				}
			} elsif ($dbtype =~ m/^snp\d+/) {
				@record = split (/\t/, $_, -1);		#-1 is required before some dbSNP records have many empty tab fields in the end
				@record == 18 or @record == 26 or die "Error: invalid record found in dbSNP database file $dbfile (18 or 26 fields expected but found ${\(scalar @record)}): <$_>\n" . join("\n",@record);
				$record[1] =~ s/^chr// or die "Error: invalid record found in DB file (2nd field should start with 'chr'): <$_>\n";
				($chr, $start, $end, $rsid, $strand, $ucscallele, $twoallele, $class) = @record[1,2,3,4,6,8,9,11];
				$start++;			#UCSC use zero-start system
				$chr =~ s/^chr//i;		#added 20130513
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				unless ($class eq 'single' or $class eq 'deletion' or $class eq 'in-del' or $class eq 'insertion') {	#enum('unknown','single','in-del','het','microsatellite','named','mixed','mnp','insertion','deletion')
					next;
				}
		
				my @allele = split (/\//, $twoallele);
				
				#before Jan 2011, only di-allelic SNPs are handled in ANNOVAR
				#@allele == 2 or next;		#many entries have no allele information (for example, rs71010435)
				#in Jan 2011 version, I decided to handle tri-allelic and quad-allelic SNP as well
				
				@allele >= 2 or next;		#Jan 2011 modification
				if ($strand eq '-') {					#handle reverse strand annotation (the vast majority of records in dbSNP should be already in + strand)
					for my $i (0 .. @allele-1) {
						$allele[$i] = revcom ($allele[$i]);
					}
					#$ucscallele = revcom ($ucscallele);		#added Jan 24, 2011 (per Kevin Ha) removed Feb 10, 2011 (per Eric Stawiski)
					#note that some SNPs (e.g., rs28434453) may have multiple location in diferent chromosome or strand; I may want to handle this by a special flag in the future
					#585        chr1    13301   13302   rs28434453      0       -       C       C       C/T     genomic single etc...
					#1367    chr15   102517867       102517868       rs28434453      0       +       G       G       C/T     genomic single etc...
				}
				
				#in-del is usually annotated below, so they require special treatment
				#587     chr1    384538  384539  rs3971283       0       +       T       T       -/ATT   genomic in-del  unknown 0       0       unknown exact   3
				if ($class eq 'in-del') {					#indel are usually annotated as -/xxx, where xxx is the alternative allele
					$obs = length ($ucscallele) . $allele[1];		#prefix a number before the alleles, indicating block substitution
					if (@allele > 2) {
						for my $i (2 .. @allele-1) {
							push @obs2, length ($ucscallele) . $allele[$i];
							push @start2, $start;
							push @score2, $rsid;
						}
					}
				} elsif ($class eq 'insertion') {
					$start--;
					$obs = "0$allele[1]";
					if (@allele > 2) {
						for my $i (2 .. @allele-1) {
							push @obs2, "0$allele[$i]";
							push @start2, $start;
							push @score2, $rsid;
						}
					}
				} elsif ($class eq 'deletion') {
					$obs = length ($ucscallele);
				} else {
					for my $i (0 .. @allele-1) {
						if ($ucscallele eq $allele[$i]) {
							@obs2 = @allele;
							splice (@obs2, $i, 1);
							for my $j (0 .. @obs2-1) {
								push @score2, $rsid;
								push @start2, $start;
							}
						}
					}
					if (@obs2) {
						$obs = shift @obs2;
						$start = shift @start2;
						$score = shift @score2;
					} else {
						$verbose and printerr ("Database error: wildtype base $ucscallele is not part of the allele description in <$_>\n");
						next;
					}
				}
				$score = $rsid;
			} elsif ($dbtype =~ m/^1000g_(\w+)/ or $dbtype =~ m/^1000g2010_(\w+)/ or $dbtype =~ m/^1000g201\d\w\w\w_(\w+)/ or $dbtype1 =~ m/[A-Z][A-Z][A-Z]\.sites.\d{4}_\d{2}/) {	#20191010: since too many users asked about FAQ #4, I decided to just allow this keyword (previously invalid) in command line
				@record = split (/\t/, $_);
				@record == 5 or @record == 6 or die "Error: invalid record found in 1000G database file $dbfile (5 or 6 fields expected): <$_>\n";
				($chr, $start, $ref, $obs, $af) = @record;			#there is no "END" in 1000G input file
				$chr =~ s/^chr//i;		#added 20130513
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				if ($maf_threshold) {
					if ($af > 0.5) {					#the frequency is the non-reference allele frequency, which could exceed 0.5
						if ($reverse) {
							1-$af <= $maf_threshold or next;
						} else {
							1-$af >= $maf_threshold or next;
						}
					} else {
						if ($reverse) {
							$af <= $maf_threshold or next;
						} else {
							$af >= $maf_threshold or next;
						}
					}
				}
				$score = $af;
				if (defined $score and defined $score_threshold and $score=~/\d{1,}\.{0,1}\d{0,}/) {
					if ($reverse) {
						$score > $score_threshold and next;
					} else {
						$score < $score_threshold and next;
					}
				}
			} elsif ($dbtype eq 'vcf') {			#vcf file is adopted by 1000 Genomes Project; it can describe both SNPs and indels, and it may contain both summary level statistics and individual level genotype calls
				($chr, $start, $rsid, $ref, $obs, $qual, $fil, $info) = split (/\t/, $_);
				$chr =~ s/^chr//i;			#added 20130513
				if ($chromosome) {
					$valichr{$chr} or next;
				}
	
				my ($ac, $an);

				if ($info =~ m/AF=([^;]+)/ or $info =~ m/AF1=([^;]+)/) {		#AC=18;AF=9.967e-03;AN=1806;DB;DP=3721
					($score, @score2) = split (/,/, $1);	#this is ideal case scenario; in reality, 1000G does not supply AF for alternative alleles
					if ($obs =~ m/(\w+),(\w+)/) {		#1000G November; this format is not really valid because it does not handle tri-allelic SNP
						($obs, @obs2) = split (/,/, $obs);
						if (@obs2 and not @score2) {
							@score2 = map {$score} @obs2;
						}
					}
				} elsif ($info =~ m/AC=([^;]+);AN=(\d+)/) {
					my ($alleles, $count) = ($1, $2);
					if ($alleles =~ m/^(\d+),(.+)/) {
						$score = sprintf ("%.3f", $1/$count);
						@score2 = split (/,/, $2);
						@score2 = map {sprintf("%.3f", $_/$count)} @score2;
						($obs, @obs2) = split (/,/, $obs);				#the obs is composed of two alleles
					} else {
						$af = sprintf ("%.3f", $alleles/$count);
						$score = $af;
						#this is an invalid record in 1000GJuly: 1       2266231 rs11589451      C       T,A     .       PASS    AA=c;AC=20;AN=120;DP=237
						if ($obs =~ m/(\w+),/) {
							$count_invalid_dbline++;
							$verbose and printerr "WARNING: Invalid input line found in $dbfile (more than one alleles are observed, but only one is annotated with allelic counts): <$_>\n";
							$obs = $1;		#2011jun: instead of using "next", I decided to still process this type of variants
						}
					}
				} else {
					$score = 'NA';
					if ($obs =~ m/^(\w+),/) {
						($obs, @obs2) = split (/,/, $obs);
						@score2 = map {'NA'} @obs2;
					}
				}
				
				if ($infoasscore) {                             #when -infoasscore is set, print out the information field as the score for VCF file
					$score = $info;
					@score2 = map {$info} @obs2;
				} elsif ($idasscore) {
					$score = $rsid;
					@score2 = map {$rsid} @obs2;
				}
				
				
				for my $i (0 .. @obs2-1) {		#if there are tri-allelic variants or quad-allelic variants or more alleles
					($start2[$i], $ref2[$i], $obs2[$i]) = reformatStartRefObs ($start, $ref, $obs2[$i]);
					defined $start2[$i] or die "Error: unable to parse VCF input line: <$_>\n";
				}

				($start, $ref, $obs) = reformatStartRefObs ($start, $ref, $obs);	#after handling obs2, then change start/ref/obs
				defined $start or die "Error: unable to parse VCF input line: <$_>\n";
			} else {
				#$dbtype eq 'generic' or print STDERR "NOTICE: the --dbtype $dbtype is assumed to be in generic ANNOVAR database format\n";
				($chr, $start, $end, $ref, $obs, $score, @otherinfo) = split (/\t/, $_);
				defined $obs or die "Error: invalid database entry (generic ANNOVAR database should contain at least 4 tab-delimited fields per line): <$_>\n";
				($ref, $obs) = (uc $ref, uc $obs);		#make sure to use upper case, as query is always in upper case
				defined $obs or die "Error: the generic database file must contains at least five tab-delimited fields per line (but observed line: $_)\n";
				defined $score or $score = "NA";
				$chr =~ s/^chr//i;		#changed 20120712; when genericdb contains "chr", all variants will be filtered.
				if ($chromosome) {
					$valichr{$chr} or next;
				}
				defined $obs or die "Error: invalid record found in DB file $dbfile (at least 5 fields expected for 'generic' dbtype): <$_>\n";
				
				if ($start == $end and $ref eq '-') {	#insertion
					$obs = "0$obs";
				} elsif ($obs eq '-') {			#deletion
					$obs = $end-$start+1;
				} elsif ($start != $end or $start==$end and length($obs)>1) {		#block substitution fixed 20130430
					$obs = ($end-$start+1) . $obs;
				}
				if (defined $score and defined $score_threshold and $score=~/\d{1,}\.{0,1}\d{0,}/) {
					if ($reverse) {
						$score > $score_threshold and next;
					} else {
						$score < $score_threshold and next;
					}
				}
			}
			if ($variant->{$chr, $start, $obs}) {
				my ($ref, @info) = split (/\n/, $variant->{$chr, $start, $obs});	#most likely, only one piece of information
				for my $i (0 .. @info-1) {
					if ($otherinfo) {
						if ($infosep) {
							map {s/,/\\x23/g} @otherinfo;
							print $DROPPED join ("\t", $dbtype, join ('#', $score, @otherinfo)), "\t", $info[$i], "\n";
						} else {
							map {s/,/\\x2c/g} @otherinfo;
							print $DROPPED join ("\t", $dbtype, join (",", $score, @otherinfo)), "\t", $info[$i], "\n";
						}
					} else {
						print $DROPPED join ("\t", $dbtype, $score), "\t", $info[$i], "\n";
					}
				}
				delete $variant->{$chr, $start, $obs};
			}
			if (@obs2) {					#this block is to handle multi-allelic variants, basically it copies the above paragraph
				for my $j (0 .. @obs2-1) {
					if ($variant->{$chr, $start2[$j], $obs2[$j]}) {
						my ($ref, @info) = split (/\n/, $variant->{$chr, $start2[$j], $obs2[$j]});	#most likely, only one piece of information
						for my $i (0 .. @info-1) {
							if (@otherinfo) {
								if ($infosep) {
									map {s/,/\\x23/g} @otherinfo;
									print $DROPPED join ("\t", $dbtype, join ('#', $score2[$j], @otherinfo)), "\t", $info[$i], "\n";
								} else {
									map {s/,/\\x2c/g} @otherinfo;
									print $DROPPED join ("\t", $dbtype, join (",", $score2[$j], @otherinfo)), "\t", $info[$i], "\n";
								}
							} else {
								print $DROPPED join ("\t", $dbtype, $score2[$j]), "\t", $info[$i], "\n";
							}
						}
						delete $variant->{$chr, $start2[$j], $obs2[$j]};
					}
				}
			}
	
			if ( $chunk_here > $chunk_max ) {
				last;
			}
		}
	}

	for my $key (keys %$variant) {
		my ($chr, $start, $obs) = split ($;, $key);		#hash key separator
		my ($ref, @info) = split (/\n/, $variant->{$key});

		for my $i (0 .. @info-1) {
			print $FIL $info[$i], "\n";
		}
	}
	printerr "Done\n";
	$count_invalid_dbline and printerr "WARNING: $count_invalid_dbline lines in dbfile $dbfile were ignored due to invalid formats\n";
}

# Reformat the $start, $ref, $obs
sub reformatStartRefObs {
	my ($start, $ref, $obs) = @_;
	if (length ($ref) == 1 and length ($obs) == 1) {#single base substitution
		1;					#the obs and obs2 is already handled
	} elsif ($obs =~ m/^\-((\w)(\w*))$/) {		#deletion (1000G March)
		$2 eq $ref or $ref eq 'N' or die "Error: mismatch of deleted allele and reference allele: <$_>\n";
		$obs = length ($1);
	} elsif ($obs =~ m/^\+(\w+)$/) {		#insertion (1000G March)
		$obs = "0$1";
	} elsif ($ref =~ m/^[ACGTN]+$/ and $obs =~ m/^[ACGTN]+$/) {
		if (length ($ref) > length ($obs)) {			#deletion or block substitution with shorter size
			my $head = substr ($ref, 0, length ($obs));
			if ($head eq $obs) {
				$start += length ($obs);
				$obs = length ($ref) - length ($obs);
				$ref = substr ($ref, length ($obs));	#added 20130513
			} else {
				$obs = length ($ref) . $obs;
			}
		} else {					#insertion or block substitution with equal size
			my $head = substr ($obs, 0, length ($ref));
			if ($head eq $ref) {
				$start += (length ($ref)-1);
				$obs = '0' . substr ($obs, length ($ref));
				$ref = '-';			#added 20130513
			} else {
				$obs = length ($ref) . $obs;
			}
		}
	} else {
		return undef;
	}
	return ($start, $ref, $obs);
}

# Perform region-based annotation with support for multi-threading, when $start_line, $end_line and $cur_thread are specified. Unlike gene or filter-based annotation which needs to process by batches, the region-based annotation is performed on each input query on the fly
sub annotateQueryByRegionThread {
	my ($regionout_file, $invalid_file, $start_line, $end_line, $cur_thread) = @_;
	my ($REGIONOUT, $INVALID, $QUERY);
	open ($REGIONOUT, ">$regionout_file") or die "Error: cannot write to output file $outfile.${buildver}_$dbtype1: $!\n";
	open ($INVALID, ">$invalid_file") or die "Error: cannot write to output file $outfile.invalid_input: $!\n";

	open ($QUERY, $queryfile) or die "Error: cannot read from --queryfile ($queryfile): $!\n";

	my $regiondb = {};
	if ($dbtype eq 'gff3') {
		($regiondb) = readGFF3RegionAnnotation ();
	} elsif ($dbtype eq 'bed') {
		($regiondb) = readBedRegionAnnotation ();
	} elsif ($dbtype eq 'generic') {
		($regiondb) = readGenericRegionAnnotation ();
	} else {
		($regiondb) = readUCSCRegionAnnotation ();
	}
	
	my ($chr, $start, $end, $ref, $obs);
	my ($linecount, $invalidcount, $invalid) = qw/0 0/;
	
	if ($thread and $start_line) {
		for (1 .. $start_line-1) {
			$_ = <$QUERY>;
		}
	}
	while (<$QUERY>) {
		if ($thread) {
			if ($. > $end_line) {				#reach to the last line for this thread
				last;
			}
		}
		
		s/[\r\n]+$//;
		if (m/^#/ and $comment) {				#comment line start with #, do not include this is $linecount
			print $REGIONOUT "#comment\t#comment\t$_\n";
			next;
		}
		
		($invalid, $chr, $start, $end, $ref, $obs) = detectInvalidInput ($_);
		
		if ($invalid) {
			print $INVALID $_, "\n";			#invalid record found
			$invalidcount++;
			next;
		}
		$linecount++;

		my $bin1 = int ($start/$genomebinsize);		#start bin
		my $bin2 = int ($end/$genomebinsize);		#end bin (usually same as start bin, unless the query is really big that spans multiple megabases)
		my ($foundhit, $score, $name);
		for my $bin ($bin1 .. $bin2) {
			for my $nextgene (@{$regiondb->{$chr, $bin}}) {
				my ($txstart, $txend, $txscore, $txname) = @$nextgene;
				
				if ($end < $txstart) {
					#db:            <------------------------->
					#query: <--->
					last;						#if genomic region is too far away from end, end the search of the bins
				} elsif ($end <= $txend) {				#query contained completely within db region
					if ($start >= $txstart) {
						#db:      <-------------------------->
						#query:       <------------------>
					} else {					#query overlap but upstream of db region
						#db:       <------------------------->
						#query: <---------------------->
						if ($minqueryfrac) {
							if (($end-$txstart+1)/($end-$start+1) < $minqueryfrac) {
								next;
							}
						}
					}
					$foundhit++;
					$score ||= $txscore; $name ||= $txname;
					if (defined $txscore) {
						if ($score < $txscore) {
							$score = $txscore;
							$name=$txname;
						}
						if ($score == $txscore and defined $name and $name ne $txname) {
							$name .= "$/$txname";
						}
					} else {
						$name .= "$/$txname";
					}
					if ($dbtype1 eq 'cytoBand') {			#a new chromosome band is encountered
						$name ne $txname and $name .= "$/$txname";
					}
				} elsif ($start <= $txend) {
					if ($start >= $txstart) {			#query overlap but downstream of db region
						#db:      <------------------------>
						#query:        <----------------------->
						if ($minqueryfrac) {
							if (($txend-$start+1)/($end-$start+1) < $minqueryfrac) {
								next;
							}
						}
					} else {
						#db region completely contained within query
						#db:      <------------------------->
						#query: <------------------------------>
						if ($minqueryfrac) {
							if (($txend-$txstart+1)/($end-$start+1) < $minqueryfrac) {
								next;
							}
						}
					}
					$foundhit++;
					$score ||= $txscore; $name ||= $txname;
					if (defined $txscore) {
						if ($score < $txscore) {
							$score = $txscore;
							$name=$txname;
						}
						if ($score == $txscore and defined $name and $name ne $txname) {
							$name .= "$/$txname";
						}
					} else {
						$name .= "$/$txname";
					}
					if ($dbtype1 eq 'cytoBand') {			#a new chromosome band is encountered
						$name ne $txname and $name .= "$/$txname";
					}
				} else {
					#query            ---
					#gene  <-*----*->
				}
			}
		}
		$linecount =~ m/000000$/ and printerr "NOTICE: Finished processing $linecount variants in queryfile\n";
		if ($foundhit) {
			$name ||= '';
			my @name = split (/$\//, $name);
			my %name = map {$_, 1} @name;
			@name = keys %name; 
				
			if ($gff3attr) {
				$name = join (";;", @name);		#when multiple attributes are available, separate them by double semicolon ;;
				print $REGIONOUT "$dbtype\t$name\t$_\n";
			} else {
				if ($dbtype1 eq 'cytoBand') {
					map {s/^chr//} @name;
					if (@name >= 2) {
						$name[$#name] =~ s/^\d+//;
						$name = $name[0] . '-' . $name[$#name];
					} else {
						$name = $name[0];
					}
					print $REGIONOUT "$dbtype\t$name\t$_\n";
				} else {
					$name = join (",", @name);
					print $REGIONOUT "$dbtype\t", $score?"Score=$score;":"", $name?"Name=$name":"", "\t", $_, "\n";
				}
			}
		}
	}
	close ($QUERY);
	close ($REGIONOUT);
	close ($INVALID);
	
	printerr "NOTICE: Finished region-based annotation on $linecount genetic variants\n";
}

# Read region annotation database stored in GFF3 files	
sub readGFF3RegionAnnotation {
	my ($dbfile);
	my ($regioncount, $dbcount) = (0, 0);
	my (@record, %regiondb);
	
	$dbfile = File::Spec->catfile ($dbloc, $gff3dbfile);
	-f $dbfile or die "Error: required database $dbfile does not exists.\n";
	
	open (DB, $dbfile) or die "Error: cannot read from database file $dbfile: $!\n";
	printerr "NOTICE: Reading annotation database $dbfile ... ";
	$_ = <DB>;
	$_ =~ m/^##gff-version\s+3/ or die "Error: invalid header line found in the GFF3 database $dbfile (expect to see '##gff-version 3'): <$_>\n";
	while (<DB>) {
		m/^#/ and next;			#skip comments line
		m/^##FASTA/ and last;		#reached the FASTA sequence section of GFF3 file
		$dbcount++;
		s/[\r\n]+$//;			#deleting the newline characters
		@record = split (/\t/, $_);
		@record == 9 or die "Error: invalid records found in the GFF3 database $dbfile (9 fields expected): <$_>\n";
		my ($chr, $start, $end, $score, $attribute) = @record[0,3,4,5,8]; 
		$chr=~s/^chr//;			#sometimes the chr prefix is present and should be removed (query usually does not contain this chr prefix)
		my $name;
		
		$score eq '.' and $score = undef;	#score in GFF3 is either a floating point number or undefined as "."
		if (defined $score_threshold and defined $score) {
			if ($reverse) {
				$score > $score_threshold and next;
			} else {
				$score < $score_threshold and next;			#if --score_threshold is set, the low scoring segment will be skipped
			}
		}
		
		if ($gff3attr) {		#use the entire 9th column (attributes column) as the name in the output
			$name = $attribute;
		} else {			#otherwise, use the ID field (which is unique for every record in the GFF3 file) as the name in the output
			my @feature = split (/;/, $attribute);
			for my $i (0 .. @feature-1) {
				if ($feature[$i] =~ m/^\s*ID=(.+?)\s*$/) {
					$name = $1;
					last;
				}
			}
			defined $name or $name = 'NA';
		}
		
		#GFF3 does not strictly require existence of ID field, so the following line is commented out on 20130208
		#defined $name or die "Error: invalid record in GFF3 database $dbfile (ID field not found): <$_>\n";
		
		my ($bin1, $bin2) = (int($start/$genomebinsize), int($end/$genomebinsize));
		for my $nextbin ($bin1 .. $bin2) {
			push @{$regiondb{$chr, $nextbin}}, [$start, $end, $score, $name];
		}
		$regioncount++;
		if ($verbose and $dbcount =~ m/000000$/) {
			my ($availmem, $allmem) = currentAvailMemory();
			printerr "NOTICE: Current system available memory is $availmem kb (this ANNOVAR program used $allmem kb)\n";
		}
	}
	close (DB);
	for my $key (keys %regiondb) {						#pre-sort gene DB by txstart to faciliate future use
		@{$regiondb{$key}} = sort {$a->[0] <=> $b->[0]} @{$regiondb{$key}};
	}
	printerr "Done with $regioncount regions from $dbcount GFF3 records\n";
	return (\%regiondb);
}

# Read region annotation database stored in BED files
sub readBedRegionAnnotation {
	my ($dbfile);
	my ($regioncount, $dbcount) = (0, 0);
	my (@record, %regiondb);
	my ($chr, $start, $end);
	
	$dbfile = File::Spec->catfile ($dbloc, $bedfile);

	-f $dbfile or die "Error: required bedfile $dbfile does not exists.\n";
	
	open (DB, $dbfile) or die "Error: cannot read from database file $dbfile: $!\n";
	printerr "NOTICE: Reading annotation database $dbfile ... ";

	while (<DB>) {
		$dbcount++;
		s/[\r\n]+$//;							#deleting the newline characters
		m/^browser\s/ and next;
		m/^track\s/ and next;
		@record = split (/\t/, $_);
		
		($chr, $start, $end) = @record;
		defined $end or die "Error: invalid record found in BED file (at least 3 tab-delimited records expected): <$_>\n";
		$start =~ m/^\d+$/ or die "Error: invalid record found in BED file (second column must be a positive integer): <$_>\n";
		$end =~ m/^\d+$/ or die "Error: invalid record found in BED file (third column must be a positive integer): <$_>\n";

		$chr =~ s/^chr//;
		$start++;										#due to the zero-opening coordinate system in UCSC
		
		my $score = '';
		for my $i (0 .. @colsWanted-1) {
			defined $record[$colsWanted[$i]-1] and $score ||= ',' . $record[$colsWanted[$i]-1];
		}
		$score =~ s/^,//;
		$score ||= 'NA';

		my ($bin1, $bin2) = (int($start/$genomebinsize), int($end/$genomebinsize));
		for my $nextbin ($bin1 .. $bin2) {
			push @{$regiondb{$chr, $nextbin}}, [$start, $end, 0, $score];
		}
		$regioncount++;
		if ($verbose and $dbcount =~ m/000000$/) {
			my ($availmem, $allmem) = currentAvailMemory();
			printerr "NOTICE: Current system available memory is $availmem kb (this ANNOVAR program used $allmem kb)\n";
		}
	}
	close (DB);
	
	for my $key (keys %regiondb) {						#pre-sort gene DB by txstart to faciliate future use
		@{$regiondb{$key}} = sort {$a->[0] <=> $b->[0]} @{$regiondb{$key}};
	}
	printerr "Done with $regioncount regions\n";
	return (\%regiondb);
}

# Read region annotation database stored in ANNOVAR generic files (five column format)
sub readGenericRegionAnnotation {
	my ($dbfile);
	my ($regioncount, $dbcount) = (0, 0);
	my (@record, %regiondb);
	my ($chr, $start, $end);
	
	$dbfile = File::Spec->catfile ($dbloc, $genericdbfile);

	-f $dbfile or die "Error: required genericdbfile $dbfile does not exists.\n";
	
	open (DB, $dbfile) or die "Error: cannot read from database file $dbfile: $!\n";
	printerr "NOTICE: Reading annotation database $dbfile ... ";

	while (<DB>) {
		$dbcount++;
		s/[\r\n]+$//;							#deleting the newline characters

		@record = split (/\t/, $_);
		
		($chr, $start, $end) = @record;
		defined $end or die "Error: invalid record found in BED file (at least 3 tab-delimited records expected): <$_>\n";
		$start =~ m/^\d+$/ or die "Error: invalid record found in BED file (second column must be a positive integer): <$_>\n";
		$end =~ m/^\d+$/ or die "Error: invalid record found in BED file (third column must be a positive integer): <$_>\n";

		$chr =~ s/^chr//;
		#$start++;										#due to the zero-opening coordinate system in UCSC (comment out this line since annovar use 1-based coordinate)
		
		my $score = '';
		for my $i (0 .. @colsWanted-1) {
			defined $record[$colsWanted[$i]-1] and $score ||= ',' . $record[$colsWanted[$i]-1];
		}
		$score =~ s/^,//;
		$score ||= 'NA';

		my ($bin1, $bin2) = (int($start/$genomebinsize), int($end/$genomebinsize));
		for my $nextbin ($bin1 .. $bin2) {
			push @{$regiondb{$chr, $nextbin}}, [$start, $end, 0, $score];
		}
		$regioncount++;
		if ($verbose and $dbcount =~ m/000000$/) {
			my ($availmem, $allmem) = currentAvailMemory();
			printerr "NOTICE: Current system available memory is $availmem kb (this ANNOVAR program used $allmem kb)\n";
		}
	}
	close (DB);
	
	for my $key (keys %regiondb) {						#pre-sort gene DB by txstart to faciliate future use
		@{$regiondb{$key}} = sort {$a->[0] <=> $b->[0]} @{$regiondb{$key}};
	}
	printerr "Done with $regioncount regions\n";
	return (\%regiondb);
}

# Read region annotation database stored in standard UCSC annotation database files
sub readUCSCRegionAnnotation {
	my ($dbfile);
	my ($regioncount, $dbcount) = (0, 0);
	my (@record, %regiondb);
	my ($chr, $start, $end, $score, $normscore, $name);
	my ($expectedLength, @positionCols, @scoreCols, @colsToOutput);
	
	if ($dbtype1 =~ m/^mce(\d+way)$/) {
		$dbfile = File::Spec->catfile ($dbloc, "${buildver}_phastConsElements$1.txt");
	} else {
		$dbfile = File::Spec->catfile ($dbloc, "${buildver}_$dbtype1.txt");
	}
	-f $dbfile or die "Error: required database $dbfile does not exists. Please use 'annotate_variation.pl -buildver $buildver -downdb $dbtype $dbloc' to download annotation database.\n";

	#################$$$
	### The following SWITCH structure is modified Jan 2011 to faciliate future expansion
	### $expectedLength is the number of cols expected in each line
	### @postionCols => location of ($chr,$start,$end) columns
	### @scoreCols => location of ($score, $normscore) columns leave empty is set not present (then set to zero below) ; WARNING must be empty or of length 2
	### @colsToOutPut => location of ($name) columns to put into $name concatinated with ":" below

	if ($dbtype1 =~ m/^phastConsElements\d+way/) {
		$expectedLength=6;
		@positionCols=(1,2,3);
		@scoreCols=(4,5);		#normalized score
		@colsToOutput=(4);		#lod=xxx is the Name output
	} elsif ($dbtype1 eq 'evofold') {
		$expectedLength=10;
		@positionCols=(1,2,3);
		@scoreCols=(5,5);
		@colsToOutput=(4);
	} elsif ($dbtype1 eq 'tfbsConsSites') {
		$expectedLength=8;
		@positionCols=(1,2,3);
		@scoreCols=(7,5);
		@colsToOutput=(4);
	} elsif ($dbtype1 eq 'wgRna') {
		$expectedLength=10;
		@positionCols=(1,2,3);
		@scoreCols=(5,5);
		@colsToOutput=(4);
	} elsif ($dbtype1 eq 'targetScanS') {
		$expectedLength=7;
		@positionCols=(1,2,3);
		@scoreCols=(5,5);
		@colsToOutput=(4);
	} elsif ($dbtype1 eq 'genomicSuperDups') {
		$expectedLength=30;
		@positionCols=(1,2,3);
		@scoreCols=(27,27);
		@colsToOutput=(4);
	} elsif ($dbtype1 eq 'omimGene') {
		$expectedLength=5;
		@positionCols=(1,2,3);
		@scoreCols=();
		@colsToOutput=(4);
	} elsif ($dbtype1 eq 'gwasCatalog') {
		$expectedLength=23;
		@positionCols=(1,2,3);
		@scoreCols=();
		@colsToOutput=(10);
	} elsif ($dbtype1 eq 'dgv') {
		$expectedLength=16;
		@positionCols=(1,2,3);
		@scoreCols=();
		@colsToOutput=(4); 
	} elsif ($dbtype1 eq 'rmsk') {
		$expectedLength=17;
		@positionCols=(5,6,7);
		@scoreCols=();
		@colsToOutput=(10); 
	} elsif ($dbtype1 eq 'cytoBand') {		#special handling required
		$expectedLength=5;
		@positionCols=(0,1,2);
		@scoreCols=();
		@colsToOutput=(0,3);
	} elsif ($dbtype1 =~ m/^chr\w+_chainSelf$/) {	#example: chr1_selfChain
		$expectedLength=13;
		@positionCols=(2,4,5);
		@scoreCols=(12,12);
		@colsToOutput=(11);
	} elsif ($dbtype1 =~ m/^chr\w+_chain\w+$/) {	#example: chr1_chainPanTro2
		$expectedLength=12;
		@positionCols=(2,4,5);
		@scoreCols=();
		@colsToOutput=(11);
	} elsif ($dbtype1 =~ m/^snp\d+/) {	
		#$expectedLength=18;			#dbSNP132 in hg19 has now 26 fields!
		$expectedLength='';
		@positionCols=(1,2,3);
		@scoreCols=();
		@colsToOutput=(4);
	} elsif ($dbtype1 eq 'gerp++elem') {
		$expectedLength = 5;			#previously this was 7, changed to 5 on 20120222
		@positionCols=(0,1,2);
		@scoreCols=();
		@colsToOutput=(3);
	} elsif ($dbtype1 eq 'stsMap') {
		$expectedLength='';
		@positionCols=(0,1,2);
		@scoreCols=();
		@colsToOutput=(3);
	} elsif ($dbtype1 =~ m/^wgEncode/) {
		$expectedLength='';
		@positionCols=(1,2,3);
		@scoreCols=(5);
		@colsToOutput=(4);
	} else {
		#other UCSC format if file is not defined above
		$expectedLength='';
		@positionCols=(1,2,3);
		@scoreCols=();
		@colsToOutput=(4);
	}
	
	if ($scorecolumn) {
		@scoreCols = ($scorecolumn, $scorecolumn);
	}
	if ($poscolumn) {
		@positionCols=split(/,/, $poscolumn);
	}
	
	open (DB, $dbfile) or die "Error: cannot read from database file $dbfile: $!\n";
	printerr "NOTICE: Reading annotation database $dbfile ... ";
	
	if ($expectedLength eq '') {		# if DB is unknown "generic format" use first line to get $expectedLength : file rewound afterwards
		my $line = <DB>;
		@record = split (/\t/, $line);
		$expectedLength=@record;
		seek (DB, 0, 0);
	};
	
	########$$ Check to see if user has defined columns to output (intergers or all allowed)
	if (defined $colsWanted) { 
		if ($colsWanted[0] eq 'all') {
			@colsToOutput= 0 .. ($expectedLength-1);
		} elsif ($colsWanted[0] eq 'none') {
			@colsToOutput = ();
		} else{
			@colsToOutput = @colsWanted;
		}
	};

	########$$ check that the columns requested exist in the current DB
	for my $i (0 .. @colsToOutput-1) {
		if ($colsToOutput[$i] > $expectedLength) {
			die "Error: The DB file $dbfile has only $expectedLength columns but output column $colsToOutput[$i] is requested by --colsWanted!\n";
		}
	}

	while (<DB>) {
		$dbcount++;
		s/[\r\n]+$//;							#deleting the newline characters
		@record = split (/\t/, $_, -1);					#-1 is required so that trailing empty strings are not discarded (this is common in dbSNP files)
		
		@record == $expectedLength or die "Error: invalid record in dbfile $dbfile ($expectedLength fields expected): <$_>\n";
		($chr, $start, $end) = @record[@positionCols];
		if (@colsToOutput) {		#I think there should always be a Name in the output column
			$name = join (':', @record[@colsToOutput]);
		}
		
		if(@scoreCols){
			($score, $normscore)=(@record[@scoreCols])
		} else{
			($score, $normscore) = qw/0 0/;
		}
		
		#########$$ Unusual exceptions for phastCons
		if ($dbtype1 =~ m/^phastConsElements\d+way/) {
			$score =~ s/^lod=// or die "Error: invalid lod score designation (no 'lod=' found) in dbfile $dbfile: <$_>\n";
		} ##lod= in the score for conservation tracks
		
		#########$$ Unusual exceptions for cytoBand
		if ($dbtype1 eq 'cytoBand' and not defined $colsWanted) {	#the name for chromosome band is concatenated as single word
			$name =~ s/://;
		}
		
		if (defined $score_threshold) {
			if ($reverse) {
				$score > $score_threshold and next;		#typically, people won't really do this
			} else {
				$score < $score_threshold and next;		#if --score_threshold is set, the low scoring segment will be skipped
			}
		}
		if (defined $normscore_threshold) {
			if ($reverse) {
				$normscore > $normscore_threshold and next;
			} else {
				$normscore < $normscore_threshold and next;		#if --normscore_threshold is set, the low scoring segment will be skipped
			}
		}
		
		$chr =~ s/^chr//;
		#$start++;										#due to the zero-opening coordinate system in UCSC
		$start =~ m/^\d+$/ or die "Error: invalid record found in region annotation database: <$_>\n";	#address warning messages when people use wrong region annotation databases 20150423
		$end =~ m/^\d+$/ or die "Error: invalid record found in region annotation database: <$_>\n";	
		$start == $end or $start++;								#changed on 20120517, because for indels, there is no need to add start by 1
		
		my ($bin1, $bin2) = (int($start/$genomebinsize), int($end/$genomebinsize));
		for my $nextbin ($bin1 .. $bin2) {
			if ($rawscore) {								#print out rawscore, rather than normalized score (default)
				$normscore = $score;
			}
			if (defined $name) {
				push @{$regiondb{$chr, $nextbin}}, [$start, $end, $normscore, $name];
			} else {			#name is not requested in the output
				push @{$regiondb{$chr, $nextbin}}, [$start, $end, $normscore];
			}
		}
		$regioncount++;
		if ($verbose and $dbcount =~ m/000000$/) {
			my ($availmem, $allmem) = currentAvailMemory();
			printerr "NOTICE: Current system available memory is $availmem kb (this ANNOVAR program used $allmem kb)\n";
		}
	}
	close (DB);
	
	for my $key (keys %regiondb) {						#pre-sort gene DB by txstart to faciliate future use
		@{$regiondb{$key}} = sort {$a->[0] <=> $b->[0]} @{$regiondb{$key}};
	}
	printerr "Done with $regioncount regions";
	if (defined $score_threshold or $normscore_threshold) {
		printerr " (that passed --score_threhsold or --normscore_threshold from a total of $dbcount regions)\n";
	} else {
		printerr "\n";
	}
	return (\%regiondb);
}

# Translate a DNA sequence to peptides, depending on whether autosomal or mitochondria is supplied. Currently it does not support any non-standard codon (but users can modify source code in the beginning of this file)
sub translateDNA {
	my ($seq, $chr1) = @_;
	my $protein = '';
	$seq = uc $seq;
	#length ($seq) % 3 == 0 or printerr "WARNING: length of DNA sequence to be translated is not multiples of 3: <length=${\(length $seq)}>\n";
	while ($seq =~ m/(...)/g) {
		defined $codon1{$1} or printerr "WARNING: invalid triplets found in DNA sequence to be translated: <$1>\n";
		if (uc $chr1 eq 'M' or uc $chr1 eq 'MT') {		#process chrM variants correctly
			$protein .= $codon1m{$1};
		} else {
			$protein .= $codon1{$1};
		}
	}
	return $protein;
}

# Translate a RNA sequence to peptides
sub translateRNA {
	my ($seq) = @_;
	my $protein = '';
	$seq = uc $seq;
	#length ($seq) % 3 == 0 or printerr "WARNING: length of RNA sequence to be translated is not multiples of 3: <length=${\(length $seq)}>\n";
	while ($seq =~ m/(...)/g) {
		defined $codonr1{$1} or printerr "WARNING: invalid triplets found in RNA sequence to be translated: <$1>\n";
		$protein .= $codonr1{$1};
	}
	return $protein;
}

# Calculate the reverse complement of a DNA sequence
sub revcom {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr/acgtACGT/tgcaTGCA/;
	return ($seq);
}

# Read a subset of sequence from a FASTA file, if the ID of the sequence is included in $refseqvar
sub readSeqFromFASTADB {
	my ($refseqvar) = @_;
	my (%seqhash);
	my $seqdbfile;
	
	#the four statements below should be condensed in the future (they are identical)
	$seqdbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1" . "Mrna.fa");

	my ($seqid, $curseq) = ('', '');
	my ($orfwarning, %badorf);		#warn that ORF is not complete; store the list of transcripts without complete ORF

	-f $seqdbfile or die "Error: FASTA sequence file $seqdbfile does not exist. Please use 'annotate_variation.pl --downdb $dbtype $dbloc' download the database.\n";
	open (SEQ, $seqdbfile) or die "Error: cannot read from seqdbfile $seqdbfile: $!\n";
	printerr "NOTICE: Reading FASTA sequences from $seqdbfile ... ";
	while (<SEQ>) {			#if there are multiple sequences for the same identifier, only the first one will be used, but the first one may not be a complete ORF for various reasons
		if (m/^>(\S+)/) {
			if ($refseqvar->{$seqid}) {
				#not defined $seqhash{$seqid} and $seqhash{$seqid} = $curseq;		#finish reading the sequence for seqid and save it (unless the sequence is already read from the file)
				if (defined $seqhash{$seqid}) {				#this transcript is already read before
					if (not $orfwarning and $badorf{$seqid}) {	#the previously identified transcript has incomplete ORF
						$seqhash{$seqid} = $curseq;
						$badorf{$seqid} = 0;
					}
				} else {
					$seqhash{$seqid} = $curseq;
				}
			}
			$seqid = $1;
			
			#20140711 (issue:multimap) to solve the one transcript multiple mapping issue I added below
			if (m/leftmost exon at (chr)?([\w\.]+):(\d+)/) {		#tomato chromosome names are SL2.40ch00
				$seqid .= "#$2#$3";
			}
			
			$curseq = '';
			if (m/does not have correct ORF annotation/) {
				$orfwarning++;
				$badorf{$seqid}++;
			} else {
				undef $orfwarning;
			}
		} else {
			if ($refseqvar->{$seqid}) {
				s/[\r\n]+$//;
				$curseq .= uc $_;					#only use upper case characters
			}
		}
	}
	if ($refseqvar->{$seqid}) {
		#not defined $seqhash{$seqid} and $seqhash{$seqid} = $curseq;
		if (defined $seqhash{$seqid}) {
			if (not $orfwarning and $badorf{$seqid}) {
				$seqhash{$seqid} = $curseq;
				$badorf{$seqid} = 0;
			}
		} else {
			$seqhash{$seqid} = $curseq;
		}
	}
	close (SEQ);
	printerr "Done with ", scalar keys %seqhash, " sequences\n";
	if (keys %seqhash < keys %$refseqvar) {
		my (@seqnotfound, @seqnotfound_example);
		for $seqid (keys %$refseqvar) {
			exists $seqhash{$seqid} or push @seqnotfound, $seqid;
		}
		printerr "WARNING: A total of ${\(scalar @seqnotfound)} sequences cannot be found in $seqdbfile\n";
		@seqnotfound_example = splice (@seqnotfound, 0, 3);
		printerr " (example: @seqnotfound_example)\n";
	}
	
	if (%badorf) {			#although these transcript sequences have been read into memory, it is best to delete them (Nov 2011 release)
		printerr "WARNING: A total of ${\(scalar keys %badorf)} sequences will be ignored due to lack of correct ORF annotation\n";
		for my $key (keys %badorf) {
			delete $seqhash{$key};
		}
	}
	
	return (\%seqhash);
}

# Read UCSC knownGene cross-reference file. This file provides gene name for UCSC knownGene transcript identifiers.
sub readKgXref {
	my ($inputfile) = @_;
	my (%gene_xref);
	open (XREF, $inputfile) or die "Error: cannot read from kgxref file $inputfile: $!\n";
	while (<XREF>) {
		m/^#/ and next;
		s/[\r\n]+$//;
		my @record = split (/\t/, $_, -1);
		@record >= 8 or die "Error: invalid record found in knownGene cross-reference file (>=8 fields expected but found ${\(scalar @record)}): <$_>\n";
		#some genes were given names that are prefixed with "Em:" which should be removed due to the presence of ":" in exonic variant annotation
		#Em:AC006547.7 Em:AC005003.4 Em:U62317.15 Em:AC008101.5 Em:AC004997.11 Em:U51561.2
		$record[4] =~ s/^Em:/Em./;
		$record[4] =~ s/\s//g;				#for example, gene name for uc001btm.2 is "BAI 2" with a space in the name (this gene is now obselete though)
		if ($gene_xref{$record[0]}) {			#BC003168 occur twice in kgxref file (OSBPL10, BC003168)
			if ($gene_xref{$record[0]} =~ m/^(BC|AK)\d+$/) {
				$gene_xref{$record[0]} = $record[4];
			}
		} else {
			$gene_xref{$record[0]} = $record[4];
		}
	}
	close (XREF);
	return (\%gene_xref);
}

# Read gene definition from UCSC gene annotation database file with UCSC-specific format
sub readUCSCGeneAnnotation {			#read RefGene annotation database from the UCSC Genome Browser, convert 0-based coordinates to 1-based coordinates
	my ($dbloc) = @_;
	my ($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes);
	my (%genedb, %geneidmap, %name2count, %cdslen, %mrnalen);
	my ($genecount, $ncgenecount) = (0, 0);
	
	my $dbfile;
	my $kgxref;
	my %iscoding;		#this gene name is a coding gene (if it has coding and noncoding transcripts, ignore all noncoding transcripts)
	
	if ($dbtype1 eq 'refGene') {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");
	} elsif ($dbtype1 eq 'knownGene') {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");
		my $kgxreffile = File::Spec->catfile($dbloc, $buildver . "_kgXref.txt");
		-f $kgxreffile or die "Error: the knownGene cross-reference file $kgxreffile does not exist. Please use 'annotate_variation.pl --downdb knownGene $dbloc' to download the database.\n";
		$kgxref = readKgXref ($kgxreffile);
	} elsif ($dbtype1 eq 'ensGene') {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");
	} else {
		$dbfile = File::Spec->catfile($dbloc, $buildver . "_$dbtype1.txt");		#added 2011feb18
		#die "FATAL ERROR: the dbype $dbtype1 is not supported in the readUCSCGeneAnnotation() subroutine.\n";		#commented 2011feb18
	}
	-f $dbfile or die "Error: The gene annotation database $dbfile does not exist. Please use 'annotate_variation.pl --downdb $dbtype $dbloc -build $buildver' to download the database.\n";

	open (GENEDB, $dbfile) or die "Error: cannot read from gene annotaion database $dbfile: $!\n";
	printerr "NOTICE: Reading gene annotation from $dbfile ... ";
	while (<GENEDB>) {
		s/[\r\n]+$//;							#deleting the newline characters
		my @record = split (/\t/, $_);

		if ($dbtype1 eq 'refGene') {
			@record == 15 or @record == 16 or die "Error: invalid record in $dbfile (expecting 15 or 16 tab-delimited fields in refGene file): <$_>\n";
			@record == 15 and unshift @record, 0;		#some refGene files have only 15 records (such as those generated by gff3ToGenePred)
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes) = @record[1..15];		#human hg18, mouse
		} elsif ($dbtype1 eq 'knownGene') {
			@record >= 11 or die "Error: invalid record in $dbfile (>=11 fields expected in knownGene file): <$_>\n";	#mm8=11, hg18=hg19=12
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend) = @record[0..9];
			$name2 = $kgxref->{$name} || $name;
		} elsif ($dbtype1 eq 'ensGene') {
			@record == 16 or die "Error: invalid record in $dbfile (expecting 16 fields in ensGene file): <$_>\n";
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes) = @record[1..15];
		} else {
			@record >= 11 or die "Error: invalid record in $dbfile (>=11 fields expected in $dbtype1 gene definition file): <$_>\n";
			($name, $chr, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, $exoncount, $exonstart, $exonend, $id, $name2, $cdsstartstat, $cdsendstat, $exonframes) = @record[1..15];
			$name2 or $name2=$name;		#changed on 20130821 (deleted "defined $name2" so that ccdsGene can be handled, by relacing gene name by name2 field)
			#die "FATAL ERROR: the --dbtype $dbtype is not supported in readUCSCGeneAnnotation() subroutine.\n";		#commented 2011feb18
		}



		#handle situations where the same transcript is mapped to several chromosomes or regions (for example, NM_019105 is mapped to chr6, chr6_cox_hap1, chr6_qbl_hap2; NM_002538 is mapped to chr5 positive and negative strand and also in chr5_h2_hap1)
		if ($chr =~ m/hap\d+$/) {
			next;			#this is a temporary solution on 2011feb19, to ignore alternative haplotype chromosomes
		}
	
		#$chr =~ s/^chr// or die "Error: invalid record found in $dbfile (chrom field not found): <$_>\n";						#UCSC always prefix "chr" to the chromosome identifier, so this is a good check to make sure that the file is the correct file
		$chr =~ s/^chr//;			#some genomes like zebrafish does not start with chr in their chromosome names.


		#for transcript-based gene definition, it is possible for the same transcript to map to multiple locations, therefore the "name" is not unique. One example is given below:
		#BROWSER | SIZE IDENTITY CHROMOSOME  STRAND    START     END              QUERY      START  END  TOTAL
		#browser |  1119  100.0%         22     +  24322314  24326106             NM_000854     1  1119  1136
		#browser |  1119   99.9%         22     -  24299601  24303393             NM_000854     1  1119  1136
		#in previous versions of ANNOVAR, this issue is solved in a random manner and some users complain about this behavior, but there is no real good solution
		#20140711 (issue:multimap): I implemented the functionality to use alternative transcript mappings in calling gene function
		#the key change is simply to use transcript+chr+xstart as the new "identifier" of the transcript so it is guaranteed to be unique
		$name = join ("#", $name, $chr, $txstart);
		
		$dbstrand eq '+' or $dbstrand eq '-' or die "Error: invalid dbstrand information found in $dbfile (dbstrand has to be + or -): <$_>\n";		#dbstrand is important to know and cannot be optional
		my @exonstart = split (/,/, $exonstart); 			#remove trailing comma
		my @exonend = split (/,/, $exonend);				#remove trailing comma
		$exoncount == @exonstart or die "Error: invalid record found in $dbfile (exoncount discordance): <$exoncount vs ${\(scalar @exonstart)}>\n";
		@exonstart == @exonend or die "Error: invalid record found in $dbfile (exonstart and exonend count discordance): <${\(scalar @exonstart)} vs ${\(scalar @exonend)}>\n";
		$txstart++; $cdsstart++; map {$_++} @exonstart;			#convert 0-based coordinate to 1-based coordinate

		#LOGIC here:
		#first calcluate mRNA length, and if the transcript maps to multiple locations with discordant mRNA length, only consider the leftmost chromosome and leftmost coordinate (because the FASTA file is sorted in this manner)

		my $cdslength = 0;
		my $mrnalength = 0;
		for my $i (0 .. @exonstart-1) {
			$mrnalength += $exonend[$i]-$exonstart[$i]+1;
		}
		for my $i (0 .. @exonstart-1) {					#this calculation is valid regardless of strand
			#$mrnalength += $exonend[$i]-$exonstart[$i]+1;
			if ($cdsstart >= $exonstart[$i] and $cdsstart <= $exonend[$i]) {
				if ($cdsend <= $exonend[$i]) {
					$cdslength = $cdsend-$cdsstart+1;
					last;
				} else {
					$cdslength += $exonend[$i]-$cdsstart+1;
					next;
				}
			}
			if ($cdslength and $cdsend < $exonstart[$i]) {
				die "FATAL ERROR: impossible scenario for $name in $dbfile (cdsend is less than exon start)";
			} elsif ($cdslength and $cdsend <= $exonend[$i]) {
				$cdslength += $cdsend-$exonstart[$i]+1;
				last;
			} elsif ($cdslength and $cdsend > $exonend[$i]) {
				$cdslength += $exonend[$i]-$exonstart[$i]+1;
			}
		}
		
		
		
		if ($cdsstart != $cdsend+1) {		#coding gene
			if (defined $mrnalen{$name} and $mrnalen{$name} != $mrnalength) {
				$verbose and printerr "WARNING: $name occurs more than once in $dbfile with different mRNA length. The first occurences with identical mRNA length will be uesd in analysis.\n";
				next;
			}
			
			if (defined $cdslen{$name} and $cdslen{$name} != $cdslength) {
				$verbose and printerr "WARNING: $name occurs more than once in $dbfile with different CDS length. The first occurences with identical CDS length will be uesd in analysis.\n";
				next;
			}
			
			$iscoding{$name2}++;		#name2 is a coding gene, and if there is a noncoding transcript, ignore such transcripts in future analysis
		} else {		#noncoding gene
			1;
		}
		
		$cdslen{$name} = $cdslength;
		$mrnalen{$name} = $mrnalength;
				
		my ($bin1, $bin2) = (int(($txstart - $neargene)/$genomebinsize), int(($txend + $neargene)/$genomebinsize));
		for my $nextbin ($bin1 .. $bin2) {
			push @{$genedb{$chr, $nextbin}}, [$name, $dbstrand, $txstart, $txend, $cdsstart, $cdsend, [@exonstart], [@exonend], $name2];
		}
		$geneidmap{$name} = $name2;
		$genecount++;
		$name2count{$name2}++;
		$cdsstart == $cdsend+1 and $ncgenecount++;			#non-coding gene has the same start and end site
	} 
	close (GENEDB);
	
	my %badgene;
	for my $key (keys %genedb) {
		my @newgenedb;
		for my $geneinfo (@{$genedb{$key}}) {
			if (not $cdslen{$geneinfo->[0]} and $iscoding{$geneinfo->[8]}) {
				$badgene{$geneinfo->[0]}++;
				$verbose and printerr "WARNING: $geneinfo->[0] will be ignored in analysis, because it is a non-coding transcript but the associated gene has another coding transcript\n";
			} else {
				push @newgenedb, $geneinfo;
			}
		}
		@{$genedb{$key}} = @newgenedb;
	}
	
	for my $key (keys %genedb) {						#pre-sort gene DB by txstart to faciliate future use
		@{$genedb{$key}} = sort {$a->[2] <=> $b->[2]} @{$genedb{$key}};
	}
	printerr "Done with $genecount transcripts (including $ncgenecount without coding sequence annotation) for ", scalar (keys %name2count), " unique genes\n";
	$verbose and %badgene and printerr "NOTICE: ", scalar (keys %badgene), " noncoding transcripts will be ignored, because their associated genes have annotated coding transcript\n";
	return (\%genedb, \%geneidmap, \%cdslen, \%mrnalen);
}


##################
# Get the cDNA sequence before the indicated start, and after the indicated end. $cDNA_pad gives the length
#   of the padding sequence.  $fs and $fs_end are used to indicate codon position of the refvarstart and 
#   refvarend.
#   !!!! NOTE !!!!!
#       - The padded sequence omits other bases in the affected codon! Therefore, these must be filled in
#         separately! To avoid this "trimming" behavior, set $fs to 0, $fs_end to 2.
#       - set @_[2] (get_pad_seq.$cDNA_pad) to < 0 to pad from beginning, and until end of gene.
##################
sub get_pad_seq {

    my ($refvarstart, $refvarend, $cDNA_pad, $fs, $end_fs, $refcdsstart, $refcdsend, $refseq) = @_;
    my ($pre_pad, $post_pad);

    # pre_pad
    if ($cDNA_pad < 0 || $refvarstart-$fs-$cDNA_pad < $refcdsstart) {
        $pre_pad = substr ($refseq, $refcdsstart-1, $refvarstart-$fs-$refcdsstart);
    }
    else {
        $pre_pad = substr ($refseq, $refvarstart-$fs-1-$cDNA_pad, $cDNA_pad);
    }

    # post_pad
    if ($cDNA_pad < 0 || $refvarend-$end_fs+3 + $cDNA_pad > $refcdsend) {
        $post_pad = substr ($refseq, $refvarend-$end_fs+2, $refcdsend - ($refvarend-$end_fs+2));
    }
    else {
        $post_pad = substr ($refseq, $refvarend-$end_fs+2, $cDNA_pad);
    }

    return ($pre_pad, $post_pad);

}

# Download user-specified databases from UCSC or ANNOVAR website or other websites specified by the user
sub downloadDB {
	my ($cwd, $msg, $sc);
	
	$cwd = Cwd::cwd();		#save current directory information to go back later
	
	-w $dbloc or die "Error: the directory $dbloc is not writable by the current user\n";
	chdir ($dbloc) or die "Error: the directory $dbloc cannot be accessed\n";
	
	my (@urlin, @filein, @fileout, %fail);		#the fail hash contains index of files that fail to be downloaded
	my $count_success;
	my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
	my $mrna_reminder;	#remind users to generate mRNA file
	if ($dbtype1 =~ m/^refGene/) {		#change refGene to m/^refGene/ so users can download refGeneWithVer as well
		if ($buildver =~ m/^hg\d+$/ and not $webfrom) {
			print STDERR "NOTICE: data retrieval from ANNOVAR repository by default\n";
			$webfrom = 'annovar';
		}
		
		if ($buildver =~ m/^hg\d+$/ and $webfrom eq 'annovar') {
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_$dbtype1.txt.gz";
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_${dbtype1}Mrna.fa.gz";
			$dbtype1 eq 'refGene' and push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_${dbtype1}Version.txt.gz";		#by default, add the version file to download list
		} else {
			push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/refGene.txt.gz";
			$mrna_reminder++;
		}
	} elsif ($dbtype1 =~ m/^knownGene/) {
		if ($buildver =~ m/^hg\d+$/ and not $webfrom) {
			print STDERR "NOTICE: data retrieval from ANNOVAR repository by default\n";
			$webfrom = 'annovar';
		}
		
		if ($buildver =~ m/^hg\d+$/ and $webfrom eq 'annovar') {
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_knownGene.txt.gz";
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_kgXref.txt.gz";
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_knownGeneMrna.fa.gz";
		} else {
			push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/knownGene.txt.gz";
			push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/kgXref.txt.gz";
			$mrna_reminder++;
		}		
	} elsif ($dbtype1 =~ m/^ensGene/) {
		if ($buildver =~ m/^hg\d+$/ and not $webfrom) {
			print STDERR "NOTICE: data retrieval from ANNOVAR repository by default\n";
			$webfrom = 'annovar';
		}
		
		if ($buildver =~ m/^hg\d+$/ and $webfrom eq 'annovar') {
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_ensGene.txt.gz";
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_ensGeneMrna.fa.gz";
		} else {
			push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/ensGene.txt.gz";
			$mrna_reminder++;
		}
	} elsif ($dbtype1 eq 'seq') {
		push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/bigZips/chromFa.zip";		#example: hg18, hg19
		push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/bigZips/chromFa.tar.gz";	#example: panTro2
		push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/bigZips/$buildver.chromFa.tar.gz";	#UCSC changed directory location for hg38, so this line is added July 2016
		push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/bigZips/$buildver.fa.gz";	#example: bosTau4
	} elsif ($dbtype1 eq '1000g' or $dbtype1 eq '1000g2010' or $dbtype1 eq '1000g2010jul') {
		$buildver eq 'hg18' or die "Error: currently the --dbtype of '$dbtype1' only support --buildver of 'hg18'\n";
		push @urlin, "http://www.openbioinformatics.org/annovar/download/hg18_$dbtype1.zip";
	} elsif ($dbtype1 =~ m/^1000g(\d{4})(\w{3})$/) {	#1000g2010nov, 1000g2011may, 1000g2012feb, 1000g2012apr, etc
		push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_$dbtype1.zip";
	} elsif ($dbtype1 eq 'null') {
		1;
	} else {
		$webfrom ||= 'ucsc';
		if ($webfrom eq 'annovar') {
			push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_$dbtype1.txt.gz";
			$dbtype1 ne 'avdblist' and push @urlin, "http://www.openbioinformatics.org/annovar/download/${buildver}_$dbtype1.txt.idx.gz";
		} elsif ($webfrom eq 'ucsc') {
			push @urlin, "http://hgdownload.cse.ucsc.edu/goldenPath/$buildver/database/$dbtype1.txt.gz";
		} else {
			push @urlin, "$webfrom/$dbtype1.txt.gz";
		}
	}
	
	
	@filein = @urlin;
	map {s/.+\///} @filein;
	@fileout = @filein;
	map {s/\.gz$//; s/\.zip$//} @fileout;
	
	if ($wget) {
		$msg = qx/wget --help 2>&1/ || '';		#collect the output of the system command
	} else {
		$msg = '';					#when --nowget is specified, do not use wget to retrieve files from Internet
	}
	if ($msg =~ m/Usage/) {
		checkProgramUpdate ("wget");
		for my $i (0 .. @urlin-1) {
			printerr "NOTICE: Downloading annotation database $urlin[$i] ... ";
			if ($verbose) {
				$sc = "wget -t 1 -T 30 -O $filein[$i] $urlin[$i]";
			} else {
				$sc = "wget -t 1 -T 30 -q -O $filein[$i] $urlin[$i]";
			}
			if (system ($sc)) {	#time-out is 30 seconds, with 1 retry attempt
				printerr "Failed\n";
				$verbose and print "WARNING: unable to execute system command: <$sc>\n";
				unlink ($filein[$i]);		#delete the temporary files generated by wget
				$fail{$i}++;
			} else {
				printerr "OK\n";
				$count_success++;
			}
		}
	} else {
		eval {
			require Net::FTP;
			require LWP::UserAgent;
		};
		if ($@) {
			printerr "WARNING: cannot retrieve remote files automatically (by 'wget' command or by standard Net::FTP/LWP::UserAgent Perl module).\n";
			printerr "Please manually download the following file, uncompress the files to $dbloc directory, then add a ${buildver}_ prefix to the file names.\n";
			printerr join ("\n", @urlin), "\n";
			exit (100);
		}
		
		checkProgramUpdate ("lwp");
		my ($http, $ftp);
		for my $i (0 .. @urlin-1) {
			printerr "NOTICE: Downloading annotation database $urlin[$i] ... ";
			if ($urlin[$i] =~ m/^http/) {
				$http = LWP::UserAgent->new (timeout=>10, show_progress=>$verbose);
				$http->env_proxy;
				
				my $response = $http->get ($urlin[$i], ':content_file'=>$filein[$i]);
				if ($response->is_success) {
					printerr "Done\n";
					$count_success++;
				} else {
					printerr "Failed\n";
					$verbose and printerr "WARNING: cannot retrieve remote files ($urlin[$i]) via LWP::UserAgent Perl module: ", $response->status_line, "\n";
					$fail{$i}++;
				}
			} elsif ($urlin[$i] =~ m#^ftp://([^\\\/]+)#) {		#for hgdownload.cse.ucsc.edu, ftp-trace.ncbi.nih.gov, ftp.ensembl.org, etc
				my $urlroot = $1;
				if ($ftp = Net::FTP->new($urlroot, Timeout=>10, Debug=>$verbose)) {
					$ftp->login("anonymous", 'anonymous@');		#these are the typical username and password for almost all FTP sites that offer public access
					$ftp->binary();
					my $url = $urlin[$i];
					$url =~ s#ftp://[\w\.\-]+/##;		#remove the URL root
					if (not $ftp->get($url)) {
						printerr "Failed\n";
						$verbose and printerr "WARNING: cannot retrieve remote file ($url) in FTP server $urlroot\n";
						$fail{$i}++;
					} else {
						printerr "Done\n";
						$count_success++;
					}
				} else {
					printerr "Failed\n";
					$verbose and printerr "WARNING: cannot retrieve remote file ($urlin[$i]) via Net::FTP Perl module\n";
					$fail{$i}++;
				}
				
			} else {
				die "Error: The URL $urlin[$i] uses an unsupported protocol. Download cannot continue\n";
			}
		}
	}
	
	$count_success and printerr "NOTICE: Uncompressing downloaded files\n";
	for my $i (0 .. @filein-1) {
		$fail{$i} and next;
		if ($filein[$i] =~ m/\.zip$/) {
			$msg = qx/unzip --help 2>&1/ || '';		#collect the output of the system command
			if ($msg =~ m/Usage/i) {
				if ($verbose) {
					system ("unzip -o $filein[$i]");
				} else {
					system ("unzip -o -q $filein[$i]");
				}
			} else {
				printerr "ERROR: unzip is not installed in your system.\nPlease manually uncompress the files (@filein) at the $dbloc directory", $dbtype1 eq 'seq'?".\n":", and rename them by adding ${buildver}_ prefix to the file names.\n";
				exit (101);
			}
			unlink ($filein[$i]);						#delete the ZIP file
		} elsif ($filein[$i] =~ m/\.tar\.gz$/) {		#panTro2 FASTA sequence is stored as tar.gz rather than zip
			$msg = qx/tar --help 2>&1/ || '';			#collect the output of the system command
			if ($msg =~ m/Usage/i or $msg =~ m/Options/i) {	#BSD-derived version of `tar` on Mac OS does not list "Usage" information (it prints "Option" instead)
				system ("tar -x -z -f $filein[$i]");
			} else {
				printerr "ERROR: tar/gunzip is not installed in your system.\nPlease manually uncompress the files (@filein) at the $dbloc directory", $dbtype1 eq 'seq'?".\n":", and rename them by adding ${buildver}_ prefix to the file names.\n";
				exit (102);
			}
		} elsif ($filein[$i] =~ m/\.gz$/) {
			$msg = qx/gunzip --help 2>&1/ || '';		#collect the output of the system command
			if ($msg =~ m/Usage/i) {
				system ("gunzip -f $filein[$i]");
			} else {
				printerr "ERROR: gunzip is not installed in your system.\nPlease manually uncompress the files (@filein) at the $dbloc directory", $dbtype1 eq 'seq'?".\n":", and rename them by adding ${buildver}_ prefix to the file names.\n";
				exit (103);
			}
		}
	}

	for my $i (0 .. @fileout-1) {
		$fail{$i} and next;							#skip the file that failed to be downloaded
		my $fileout = $fileout[$i];
		$dbtype1 eq 'seq' and next;					#the zip file contains dozens of FASTA files so cannot rename them automatically
		if (not $fileout =~ m/^${buildver}_/) {		#if the buildver is not the prefix of the files
			rename ($fileout, "${buildver}_$fileout") or die "Error: cannot rename $fileout to ${buildver}_$fileout\n";
			$fileout = "${buildver}_$fileout";
		}
		if (not $fileout =~ m/\.txt$/ and not $fileout =~ m/\.fa$/ and not $fileout =~ m/\.idx$/) {
			rename ($fileout, "$fileout.txt");
		}
	}
	
	$count_success and printerr "NOTICE: Finished downloading annotation files for $buildver build version, with files saved at the '$dbloc' directory\n";
	$cwd and chdir ($cwd);
	if (%fail) {
		my @failindex = keys %fail;
		if ($dbtype1 eq 'seq' and @failindex == 2) {	#not really a fail, because for seq, ANNOVAR attempts on tar.gz and zip file
			1;
		} else {
			printerr "WARNING: Some files cannot be downloaded, including ", join (', ', @urlin[@failindex]), "\n";
		}
		

	}
	
	if ($mrna_reminder) {
		printerr "---------------------------ADDITIONAL PROCEDURE---------------------------\n";
		printerr "--------------------------------------------------------------------------\n";
		printerr "NOTICE: the FASTA file for the genome is not available to download but can be generated by the ANNOVAR software.\n";
		printerr "PLEASE RUN THE FOLLOWING TWO COMMANDS CONSECUTIVELY TO GENERATE THE FASTA FILES (you may need to change -seqdir to -seqfile for some genomes):\n\n";
		printerr "\tannotate_variation.pl --buildver $buildver --downdb seq $dbloc/${buildver}_seq\n";
		printerr "\tretrieve_seq_from_fasta.pl $dbloc/${buildver}_$dbtype1.txt -seqdir $dbloc/${buildver}_seq -format $dbtype1 -outfile $dbloc/${buildver}_${dbtype1}Mrna.fa\n";
		printerr "--------------------------------------------------------------------------\n";
		printerr "--------------------------------------------------------------------------\n";
	}
	
	if (not $count_success and not $dbtype1 eq 'null') {
		exit (1);	#exit to show that nothing was downloaded
	}
}

# Check currently available memory in the system. This function is ugly written and I may consider replace it in the future, but at least it works now
sub currentAvailMemory {
	my ($availmem, $allmem) = (0, 0);
	if ($^O eq "MSWin32") {		#no easy solution to get available memory from Windows.
		($availmem, $allmem) = (0, 0);
	} elsif ($^O eq 'linux' or $^O eq 'aix' or $^O eq 'solaris') {
		if (open (TOP, "top -b -n 1 2>&1 |")) {
			my $index;
			while (<TOP>) {
				if (m/^Mem:.+\s(\d+)k free/) {
					$availmem = $1;
				}
				s/^\s+//;
				my @field = split (/\s+/, $_);
				@field >= 10 or next;			#make sure that the PID lines are reached
				if ($field[0] eq 'PID') {
					for my $i (0 .. @field-1) {
						$field[$i] eq 'RES' and $index = $i;
					}
				}
				if ($field[0] eq $$) {
					defined $index or die "Error: invalid output from top command: the line with PID and RES is not found\n";
					$allmem = $field[$index];
					if ($allmem =~ m/^([\d\.]+)(\w)$/) {
						if ($2 eq 'g') {
							$allmem = $1 * 1_000_000;
						} elsif ($2 eq 'm') {
							$allmem = $1 * 1_000;
						} elsif ($2 eq 'k') {
							$allmem = $1;
						} else {
							printerr "WARNING: unrecognizable output from top command: <$_>\n";
						}
					}
					last;
				}
			}
		}
	} else {
		($availmem, $allmem) = (0, 0);
	}
	return ($availmem, $allmem);
}

# Print to both STDERR and a LOG file
sub printerr {
	print STDERR @_;
	print LOG @_;
}

# Check whether update for the ANNOVAR software is available and print out what features are included in the update. Whenever user specifies -downdb, this subroutine is automatically executed
sub checkProgramUpdate {
	my ($method) = @_;
	my $sc;
	my ($curdate, $webdate, $webdate1);
	my (@webcontent);
	$method eq 'wget' or $method eq 'lwp' or die "Error: update checking method can be only 'wget' or 'lwp'";
	printerr "NOTICE: Web-based checking to see whether ANNOVAR new version is available ... ";
	$DATE =~ m/Date: (\d+)\-(\d+)-(\d+)/ or printerr "Failed\n" and return; 
	$curdate = $1.$2.$3;
	if ($method eq 'wget') {
		$sc = "wget -t 1 -T 30 -q -O .annovar_date http://www.openbioinformatics.org/annovar/download/annovar_date";
		if (system ($sc)) {
			printerr "Failed\n";
			return;
		} else {
			if (not open (AVDATE, ".annovar_date")) {
				printerr "Cannot access version information\n";
			} else {
				printerr "Done\n";
				@webcontent = <AVDATE>;		#$LAST_CHANGED_DATE =	'$Date: 2020-06-07 23:56:37 -0400 (Sun,  7 Jun 2020) $';
				close (AVDATE);
				unlink (".annovar_date");
			}
		}
	} elsif ($method eq 'lwp') {
		my $http = LWP::UserAgent->new (timeout=>10);
		$http->env_proxy;
		my $response = $http->get("http://www.openbioinformatics.org/annovar/download/annovar_date");
		if ($response->is_success) {
			printerr "Done\n";
			$_ = $response->decoded_content;
			@webcontent = split (/\n/, $_);
		} else {
			printerr "Failed\n";
			return;
		}
	}
	
	$webdate = $webcontent[0];
	$webdate =~ s/[\r\n]+$//;
	$webdate1 = $webdate;
	$webdate1 =~ s/\-//g;			#remove the - sign in webdate
	if ($curdate < $webdate1) {
		printerr "----------------------------UPDATE AVAILABLE------------------------------\n";
		printerr "--------------------------------------------------------------------------\n";
		printerr "WARNING: A new version of ANNOVAR (dated $webdate) is available!\n";
		printerr "         Download from http://www.openbioinformatics.org/annovar/\n";
	
		if (@webcontent >= 2) {
			printerr "Changes made in the $webdate version:\n";
			for my $i (1 .. @webcontent-1) {
				if ($webcontent[$i] =~ m/^(\d{4})\-(\d{2})\-(\d{2})[\r\n]+$/) {
					$webdate = "$1-$2-$3";
					$webdate1 = "$1$2$3";
					if ($curdate >= $webdate1) {	#the current version is more recent than this date
						last;
					} else {
						printerr "Changes made in the $webdate version:\n";
					}
				} else {
					printerr "         * $webcontent[$i]";
				}
			}
		}
		printerr "--------------------------------------------------------------------------\n";
		printerr "--------------------------------------------------------------------------\n";
	}
}

sub hasOverlap {	#deteremine whether two regions have overlaps or not
	my ($s1, $e1, $s2, $e2) = @_;
	if ($s1 > $s2) {
		($s1, $e1, $s2, $e2) = ($s2, $e2, $s1, $e1);
	}
	if ($e1 >= $s2) {	#ending of region 1 is more than start of region 2
		return 1;
	} else {
		return 0;
	}
}

=head1 SYNOPSIS

 annotate_variation.pl [arguments] <query-file|table-name> <database-location>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
        
        Arguments to download databases or perform annotations
            --downdb			download annotation database
            --geneanno			annotate variants by gene-based annotation (infer functional consequence on genes)
            --regionanno		annotate variants by region-based annotation (find overlapped regions in database)
            --filter			annotate variants by filter-based annotation (find identical variants in database)
        
        Arguments to control input and output
            --outfile <file>		output file prefix
            --webfrom <string>		specify the source of database (ucsc or annovar or URL) (downdb operation)
            --dbtype <string>		specify database type
            --buildver <string>		specify genome build version (default: hg18 for human)
            --time			print out local time during program run
            --comment			print out comment line (those starting with #) in output files 
            --exonsort			sort the exon number in output line (gene-based annotation)
            --transcript_function	use transcript name rather than gene name (gene-based annotation)
            --hgvs			use HGVS format for exonic annotation (c.122C>T rather than c.C122T) (gene-based annotation)
            --separate			separately print out all functions of a variant in several lines (gene-based annotation)
            --seq_padding		create a new file with cDNA sequence padded by this much either side (gene-based annotation)
            --(no)firstcodondel		treat first codon deletion as wholegene deletion (default: ON) (gene-based annotation)
            --aamatrix <file>		specify an amino acid substitution matrix file (gene-based annotation)
            --colsWanted <string>	specify which columns to output by comma-delimited numbers (region-based annotation)
            --scorecolumn <int>		the column with scores in DB file (region-based annotation)
            --poscolumn <string>	the comma-delimited column with position information in DB file (region-based annotation)
            --gff3dbfile <file>		specify a DB file in GFF3 format (region-based annotation)
            --gff3attribute		output all fields in GFF3 attribute (default: ID and score only)
            --bedfile <file>		specify a DB file in BED format file (region-based annotation)
            --genericdbfile <file>	specify a DB file in generic format (filter-based annotation)
            --vcfdbfile <file>		specify a DB file in VCF format (filter-based annotation)
            --otherinfo			print out additional columns in database file (filter-based annotation)
            --infoasscore		use INFO field in VCF file as score in output (filter-based annotation)
            --idasscore			use ID field in VCF file as score in output (filter-based annotation)
            --infosep			use # rather than , to separate fields when -otherinfo is used

        
        Arguments to fine-tune the annotation procedure
            --batchsize <int>		batch size for processing variants per batch (default: 5m)
            --genomebinsize <int>	bin size to speed up search (default: 100k for -geneanno, 10k for -regionanno)
            --expandbin <int>		check nearby bin to find neighboring genes (default: 2m/genomebinsize)
            --neargene <int>		distance threshold to define upstream/downstream of a gene
            --exonicsplicing		report exonic variants near exon/intron boundary as 'exonic;splicing' variants
            --score_threshold <float>	minimum score of DB regions to use in annotation
            --normscore_threshold <float> minimum normalized score of DB regions to use in annotation
            --reverse			reverse directionality to compare to score_threshold
            --rawscore			output includes the raw score (not normalized score) in UCSC Browser Track
            --minqueryfrac <float>	minimum percentage of query overlap to define match to DB (default: 0)
            --splicing_threshold <int>	distance between splicing variants and exon/intron boundary (default: 2)
            --indel_splicing_threshold <int>	if set, use this value for allowed indel size for splicing variants (default: --splicing_threshold)
            --maf_threshold <float>	filter 1000G variants with MAF above this threshold (default: 0)
            --sift_threshold <float>	SIFT threshold for deleterious prediction for -dbtype avsift (default: 0.05)
            --precedence <string>	comma-delimited to specify precedence of variant function (default: exonic>intronic...)
            --indexfilter_threshold <float>	controls whether filter-based annotation use index if this fraction of bins need to be scanned (default: 0.9)
            --thread <int>		use multiple threads for filter-based annotation
            --maxgenethread <int>	max number of threads for gene-based annotation (default: 6)
            --mingenelinecount <int>	min line counts to enable threaded gene-based annotation (default: 1000000)
       
       Arguments to control memory usage
            --memfree <int>		ensure minimum amount of free system memory (default: 0)
            --memtotal <int>		limit total amount of memory used by ANNOVAR (default: 0, unlimited, in the order of kb)
            --chromosome <string>	examine these specific chromosomes in database file
            

 Function: annotate a list of genetic variants against genome annotation 
 databases stored at local disk.
 
 Example: #download annotation databases from ANNOVAR or UCSC and save to humandb/ directory
 	  annotate_variation.pl -downdb -webfrom annovar refGene humandb/
 	  annotate_variation.pl -downdb -buildver mm9 refGene mousedb/
 	  annotate_variation.pl -downdb -buildver hg19 -webfrom annovar esp6500siv2_all humandb/
 	
 	  #gene-based annotation of variants in the varlist file (by default --geneanno is ON)
 	  annotate_variation.pl -geneanno -buildver hg19 ex1.avinput humandb/
 	  
 	  #region-based annotate variants
 	  annotate_variation.pl -regionanno -buildver hg19 -dbtype cytoBand ex1.avinput humandb/
 	  annotate_variation.pl -regionanno -buildver hg19 -dbtype gff3 -gff3dbfile tfbs.gff3 ex1.avinput humandb/
 	  
 	  #filter rare or unreported variants (in 1000G/dbSNP) or predicted deleterious variants
 	  annotate_variation.pl -filter -dbtype 1000g2015aug_all -maf 0.01 ex1.avinput humandb/
 	  annotate_variation.pl -filter -buildver hg19 -dbtype snp138 ex1.avinput humandb/
 	  annotate_variation.pl -filter -dbtype dbnsfp35a -buildver hg38 ex1.avinput humandb/
 	  annotate_variation.pl -filter -dbtype gnomad211_exome -buildver hg19 ex1.avinput humandb/
 
 Version: $Date: 2020-06-07 23:56:37 -0400 (Sun,  7 Jun 2020) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--downdb>

download annotation databases from UCSC Genome Browser, Ensembl, 1000 Genomes 
Project, ANNOVAR website or other resources. The annotation databases are 
required for functional annotation of genetic variants.

=item B<--geneanno>

perform gene-based annotation. For each variant, examine whether it hit exon, 
intron, intergenic region, or close to a transcript, or hit a non-coding RNA 
gene, or is located in a untranslated region (see *.variant_function output 
file). In addition, for an exonic variant, determine whether it causes splicing 
change, non-synonymous amino acid change, synonymous amino acid change or 
frameshift changes (see *.exonic_variant_function output file).

=item B<--regionanno>

perform region-based annotation. For each variant, examine whether its genomic 
region (one or multiple base pairs) overlaps with a specific genomic region, 
such as the most conserved elements, the predicted transcription factor binding 
sites, the specific cytogeneic bands, the evolutionarily conserved RNA secondary 
structures.

=item B<--filter>

perform filter-based annotation. For each variants, filter it against a 
variation database, such as the 1000 Genomes Project database, to identify whether it 
has been reporte in the database. Exact match of nucleotide position and nucleotide 
composition are required.

=item B<--outfile>

specify the output file prefix. Several output files will be generated using 
this prefix and different suffixes. A directory name can also be specified as 
part of the argument, so that the output files can be written to a different 
directory than the current directory.

=item B<--webfrom>

specify the source of database (ucsc or annovar or URL) in the downdb operation. 
By default, files from UCSC Genome Browser annotation database will be 
downloaded.

=item B<--dbtype>

specify the database type to be used in gene-based, region-based or filter-based 
annotations. For gene-based annotation, by default refGene annotations from the 
UCSC Genome Browser will be used for annotating variants. However, users can 
switch to use Ensembl annotations, or use the UCSC Gene annotations, or the 
GENCODE Gene annotations, or other types of gene annotations. For region-based 
annotations, users can select any UCSC annotation databases (by providing the 
database name), or alternatively select a Generic Feature Format version 3 
(GFF3) formatted file for annotation (by providing 'gff3' as the --dbtype and 
providing the --gff3dbfile argument), or select a BED file (by providing '--
dbtype bed' and --bedfile arguments). For filter-based annotations, users can 
select a dbSNP file, a 1000G file, a generic format file (with simple columns 
including chr, start, end, reference, observed, score), a VCF format file (which is 
a widely used format for variants exchange), or many other types of formats.

=item B<--buildver>

genome build version to use. By default, the hg18 build for human genome is 
used. The build version will be used by ANNOVAR to identify corresponding database files 
automatically, for example, when gene-based annotation is used for hg18 build, 
ANNOVAR will search for the hg18_refGene.txt file, but if the hg19 is used as --
buildver, ANNOVAR will examine hg19_refGene.txt instead.

=item B<--time>

print out the local time during execution of the program

=item B<--comment>

specify that the program should include comment lines in the output files. 
Comment lines are defined as any line starting with #. By default, these lines 
are not recognized as valid ANNOVAR input and are therefore written to the 
INVALID_INPUT file. This argument can be very useful to keep columns headers in 
the output file, if the input file use comment line to flag the column headers 
(usually the first line in the input file).

=item B<--exonsort>

sort the exon number in output line in the exonic_variant_function file during 
gene-based annotation. If a mutation affects multiple transcripts, the ones with 
the smaller exon number will be printed before the transcript with larger exon 
number in the output.

=item B<--transcript_function>

use transcript name rather than gene name in output, for gene-based annotation

=item B<--hgvs>

use HGVS format for exonic annotation (c.122C>T rather than c.C122T) for gene-based annotation

=item B<--separate>

for gene-based annotation, separate the effects of each variant, so that each 
effect (intronic, exonic, splicing) is printed in one output line. By default, 
all effects are printed in the same line, in the comma-separated form of 
'UTR3,UTR5' or 'exonic,splicing'.

=item B<--seq_padding>

create a new file with cDNA sequence padded by this much either side (gene-based annotation)

=item B<--firstcodondel>

if the first codon of a gene is deleted, then the whole gene will be treated as 
deleted in gene-based annotation. By default, this option is ON.

=item B<--aamatrixfile>

specify an amino acid substitution matrix, so that the scores are printed in the 
exonic_variant_function file in gene-based annotation. The matrix file is tab-
delimited, and an example is included in the ANNOVAR package.

=item B<--colsWanted>

specify which columns are desired in the output for -regionanno. By default, 
ANNOVAR inteligently selects the columns based on the DB type. However, users 
can use a list of comma-delimited numbers, or use 'all', or use 'none', to 
request custom output columns.

=item B<--scorecolumn>

specify the the column with desired output scores in UCSC database file (for 
region-based annotation). The default usually works okay.

=item B<--poscolumn>

the comma-delimited column with position information in DB file (region-based 
annotation). The default usually works okay.

=item B<--gff3dbfile>

specify the GFF3-formatted database file used in the region-based annotation. 
Please consult http://www.sequenceontology.org/resources/gff3.html for detailed 
description on this file format. Note that GFF3 is generally not compatible with 
previous versions of GFF.

=item B<--gff3attribute>

output should contain all fields in GFF3 file attribute column (the 9th column). 
By default, only the ID in the attribute and the scores for the GFF3 file will 
be printed.

=item B<--bedfile>

specify a DB file in BED format file in region-based annotation. Please consult 
http://genome.ucsc.edu/FAQ/FAQformat.html#format1 for detailed descriptions on 
this format.

=item B<--genericdbfile>

specify the generic format database file used in the filter-based annotation.

=item B<--vcfdbfile>

specify the database file in VCF format in the filter-based annotation. VCF has 
been a popular format for summarizing SNP and indel calls in a population of 
samples, and has been adopted by 1000 Genomes Project in their most recent data 
release.

=item B<--otherinfo>

print out additional columns in database file in filter-based annotation. This 
argument is useful when the annotation database contains more than one 
annotation columns, so that all columns will be printed out and separated by 
comma (by default).

=item B<--idasscore>

when annotating against a VCF file, treat the ID field in VCF file as the 
score to be printed in the output, in filter-based annotation. By default the 
score is the allele frequency inferred from VCF file.

=item B<--infoasscore>

when annotating against a VCF file, treat the INFO field in VCF file as the 
score to be printed in the output, in filter-based annotation. By default the 
score is allele frequency inferred from VCF file.

=item B<--infosep>

use '#' rather than ',' to separate multiple fields when -otherinfo is used in 
annotation. This argument is useful when the annotation string itself contains 
comma, to help users clearly separate different annotation fields.

=item B<--batchsize>

this argument specifies the batch size for processing variants by gene-based 
annotation. Normally 5 million variants (usually one human genome will have 
about 3-5 million variants depending on ethnicity) are annotated as a batch, to 
reduce the amounts of memory. The users can adjust the parameters: larger values 
make the program slightly faster, at the expense of slightly larger memory 
requirements. In a 64bit computer, the default settings usually take 1GB memory 
for gene-based annotation for human genome for a typical query file, but this 
depends on the complexity of the query (note that the query has a few required 
fields, but may have many optional fields and those fields need to be read and 
kept in memory).

=item B<--genomebinsize>

the bin size of genome to speed up search. By default 100kb is used for gene-
based annotation, so that variant annotation focused on specific bins only 
(based on the start-end site of a given variant), rather than searching the 
entire chromosomes for each variant. By default 10kb is used for region-based 
annotation. The filter-based annotations look for variants directly so no bin is 
used.

=item B<--expandbin>

expand bin to both sides to find neighboring genes/regions. For gene-based 
annotation, ANNOVAR tries to find nearby genes for any intergenic variant, with 
a maximum number of nearby bins to search. By default, ANNOVAR will 
automatically set this argument to search 2 megabases to the left and right of 
the variant in genome.

=item B<--neargene>

the distance threshold to define whether a variant is in the upstream or 
downstream region of a gene. By default 1 kilobase from the start or end site of 
a transcript is defined as upstream or downstream, respectively. This is useful, 
for example, when one wants to identify variants that are located in the 
promoter regions of genes across the genome.

=item B<--exonicsplicing>

report exonic variants near exon/intron boundary as 'exonic;splicing' variants. 
These variants are technically exonic variants, but there are some literature 
reports that some of them may also affect splicing so a keyword is preserved 
specifically for them.

=item B<--score_threshold>

the minimum score to consider when examining region-based annotations on UCSC 
Genome Browser tables. Some tables do not have such scores and this argument 
will not be effective.

=item B<--normscore_threshold>

the minimum normalized score to consider when examining region-based annotations 
on UCSC Genome Browser tables. The normalized score is calculated by UCSC, 
ranging from 0 to 1000, to make visualization easier. Some tables do not have 
such scores and this argument will not be effective.

=item B<--reverse>

reverse the criteria for --score_threshold and --normscore_threshold. So the 
minimum score becomes maximum score for a result to be printed.

=item B<--rawscore>

for region-based annotation, print out raw scores from UCSC Genome Browser 
tables, rather than normalized scores. By default, normalized scores are printed 
in the output files. Normalized scores are compiled by UCSC Genome Browser for 
each track, and they usually range from 0 to 1000, but there are some 
exceptions.

=item B<--minqueryfrac>

The minimum fraction of overlap between a query and a database record to decide 
on their match. By default, any overlap is regarded as a match, but this may not 
work best when query consist of large copy number variants.

=item B<--splicing_threshold>

distance between splicing variants and exon/intron boundary, to claim that a 
variant is a splicing variant. By default, 2bp is used. ANNOVAR is relatively 
more stringent than some other software to claim variant as regulating splicing. 
In addition, if a variant is an exonic variant, it will not be reported as 
splicing variant even if it is within 2bp to an exon/intron boundary.

=item B<--indel_splicing_threshold>

If set, max size of indel allowed to be called a splicing variant (if boundary within
--splicing_threshold bases of an intron/exon junction.) If not set, this is equal to
the --splicing_threshold, as per original behavior.

=item B<--maf_threshold>

the minor allele frequency (MAF) threshold to be used in the filter-based 
annotation for the 1000 Genomes Project databases. By default, any variant 
annotated in the 1000G will be used in filtering.

=item B<--sift_threshold>

the default SIFT threshold for deleterious prediction for -dbtype avsift 
(default: 0.05). This argument is obselete, since the recommended database for 
SIFT annotation is LJB database now, rather than avsift database.

=item B<--thread>

specify the number of threads to use in filter-based annotation. The Perl and
all components in the system needs to support multi-threaded analysis to use
this feature. It is recommended when your database is stored at a SSD drive,
which results in nearly linear speed up of annotation for large genome files.

=item B<--maxgenethread>

specify the maximum number of threads for gene-based annotation (default: 6).
Generally speaking, too many threads for gene-based annotation will negatively
impacts the performance.

=item B<--mingenelinecount>

specify the minimum line counts to enable threaded gene-based annotation
(default: 1000000). For input files with less lines, the threaded annotation
will not be used, since it actually cost more time than non-threaded annotation.


=item B<--memfree>

the minimum amount of free system memory that ANNOVAR should ensure to have.

=item B<--memtotal>

the total amount of memory that ANNOVAR should use at most. By default, this 
value is zero, meaning that there is no limit on that. Decreasing this threshold 
reduce the memory requirement by ANNOVAR, but may increase the execution time.

=item B<--chromosome>

examine these specific chromosomes in database file. The argument takes comma-
delimited values, and the dash can be correctly recognized. For example, 5-10,X 
represent chromosome 5 through chromosome 10 plus chromosome X.

=back

=head1 DESCRIPTION

ANNOVAR is a software tool that can be used to functionally annotate a list of 
genetic variants, possibly generated from next-generation sequencing 
experiments. For example, given a whole-genome resequencing data set for a human 
with specific diseases, typically around 3 million SNPs and around half million 
insertions/deletions will be identified. Given this massive amounts of data (and 
candidate disease- causing variants), it is necessary to have a fast algorithm 
that scans the data and identify a prioritized subset of variants that are most 
likely functional for follow-up Sanger sequencing studies and functional assays.

Currently, these various types of functional annotations produced by ANNOVAR can 
be (1) gene-based annotations (the default behavior), such as exonic variants, 
intronic variants, intergenic variants, downstream variants, UTR variants, 
splicing site variants, stc. For exonic variants, ANNOVAR will try to predict 
whether each of the variants is non-synonymous SNV, synonymous SNV, 
frameshifting change, nonframeshifting change. (2) region-based annotation, to 
identify whether a given variant overlaps with a specific type of genomic 
region, for example, predicted transcription factor binding site or predicted 
microRNAs.(3) filter-based annotation, to filter a list of variants so that only 
those not observed in variation databases (such as 1000 Genomes Project and 
dbSNP) are printed out.

Detailed documentation for ANNOVAR should be viewed in ANNOVAR website 
(http://annovar.openbioinformatics.org/). Below is description on commonly 
encountered file formats when using ANNOVAR software.

=over 8

=item * B<variant file format>

A sample variant file contains one variant per line, with the fields being chr, 
start, end, reference allele, observed allele, other information. The other 
information can be anything (for example, it may contain sample identifiers for 
the corresponding variant.) An example is shown below:

	16      49303427        49303427        C       T       rs2066844       R702W (NOD2)
	16      49314041        49314041        G       C       rs2066845       G908R (NOD2)
	16      49321279        49321279        -       C       rs2066847       c.3016_3017insC (NOD2)
	16      49290897        49290897        C       T       rs9999999       intronic (NOD2)
	16      49288500        49288500        A       T       rs8888888       intergenic (NOD2)
	16      49288552        49288552        T       -       rs7777777       UTR5 (NOD2)
	18      56190256        56190256        C       T       rs2229616       V103I (MC4R)

=item * B<database file format: UCSC Genome Browser annotation database>

Most but not all of the gene annotation databases are directly downloaded from 
UCSC Genome Browser, so the file format is identical to what was used by the 
genome browser. The users can check Table Browser (for example, human hg18 table 
browser is at http://www.genome.ucsc.edu/cgi-bin/hgTables?org=Human&db=hg18) to 
see what fields are available in the annotation file. Note that even for the 
same species (such as humans), the file format might be different between 
different genome builds (such as between hg16, hg17 and hg18). ANNOVAR will try 
to be smart about guessing file format, based on the combination of the --
buildver argument and the number of columns in the input file. In general, the 
database file format should not be something that users need to worry about.

=item * B<database file format: GFF3 format for gene-based annotations)>

As of June 2010, ANNOVAR cannot perform gene-based annotations using GFF3 input 
files, and any annotations on GFF3 is region-based. I suggest that users 
download gff3ToGenePred tool from UCSC and convert GFF3-based gene annotation to 
UCSC format, so that ANNOVAR can perform gene-based annotation for your species 
of interests.

=item * B<database file format: GFF3 format for region-based annotations)>

Currently, region-based annotations can support the Generic Feature Format 
version 3 (GFF3) formatted files. The GFF3 has become the de facto golden 
standards for many model organism databases, such that many users may want to 
take a custom annotation database and run ANNOVAR on them, and it would be the 
most convenient if the custom file is made with GFF3 format. 

=item * B<database file format: generic format for filter-based annotations)>

The 'generic' format is designed for filter-based annotation that looks for 
exact variants. The format is almost identical to the ANNOVAR input format, with 
chr, start, end, reference allele, observed allele and scores (higher scores are 
regarded as better).

=item * B<database file format: VCF format for filter-based annotations)>

ANNOVAR can directly interrogate VCF files as database files. A VCF file 
may contain summary information for variants (for example, this variant has MAF 
of 5% in this population), or it may contain the actual variant calls for each 
individual in a specific population.

=item * B<sequence file format>

ANNOVAR can directly examine FASTA-formatted sequence files. For mRNA sequences, 
the name of the sequences are the mRNA identifier. For genomic sequences, the 
name of the sequences in the files are usually chr1, chr2, chr3, etc, so that 
ANNOVAR knows which sequence corresponds to which chromosome. Unfortunately, 
UCSC uses things like chr6_random to annotate un-assembled sequences, as opposed 
to using the actual contig identifiers. This causes some issues (depending on 
how reads alignment algorithms works), but in general should not be something 
that user need to worry about. If the users absolutely care about the exact 
contigs rather than chr*_random, then they will need to re-align the short reads 
at chr*_random to a different FASTA file that contains the contigs (such as the 
GRCh36/37/38), and then execute ANNOVAR on the newly identified variants.

=item * B<invalid input>

If the query file contains input lines with invalid format, ANNOVAR will skip 
such line and continue with the annotation on next lines. These invalid input 
lines will be written to a file with suffix invalid_input. Users should manually 
examine this file and identify sources of error.

=back

--------------------------------------------------------------------------------

ANNOVAR is free for academic, personal and non-profit use.

For questions or comments, please contact $Author: Kai Wang <kaichop@gmail.com> $. 

=cut
