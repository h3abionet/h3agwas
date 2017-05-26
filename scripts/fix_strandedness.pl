#!/usr/bin/perl
# fix_strandedness.pl Fixes the strand of map and ped files using illumina data
# and is released under the terms of the GNU GPL version 3, or any
# later version, at your option. See the file README and COPYING for
# more information.
# Copyright 2017 by Don Armstrong <don@donarmstrong.com>.


use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

=head1 NAME

fix_strandedness.pl - Fixes the strand of map and ped files using illumina data

=head1 SYNOPSIS

fix_strandedness.pl [options]

 Options:
   --debug, -d debugging level (Default 0)
   --help, -h display this help
   --man, -m display manual

=head1 OPTIONS

=over

=item B<--debug, -d>

Debug verbosity. (Default 0)

=item B<--help, -h>

Display brief usage information.

=item B<--man, -m>

Display this manual.

=back

=head1 EXAMPLES

fix_strandedness.pl

=cut

use Scalar::Util qw(looks_like_number);
use Bio::DB::Fasta;
use vars qw($DEBUG);
use Text::LevenshteinXS qw(distance);

my %options = (debug           => 0,
               help            => 0,
               man             => 0,
              );

GetOptions(\%options,
           'ped=s',
           'map=s',
           'vcf=s',
           'ref=s',
           'merge=s',
           'illumina_data|illumina-data=s',
           'ped_out|ped-out=s',
           'map_out|map-out=s',
           'debug|d+','help|h|?','man|m');

pod2usage() if $options{help};
pod2usage({verbose=>2}) if $options{man};

$DEBUG = $options{debug};

my @USAGE_ERRORS;
for (qw(merge vcf ped map ped_out map_out illumina_data ref)) {
    if (not defined $options{$_}) {
        my $a = $_;
        $a =~ s/_/-/g;
        push @USAGE_ERRORS,"You must provide the --".$a." option";
    }
}

pod2usage(join("\n",@USAGE_ERRORS)) if @USAGE_ERRORS;


my $map = open_compressed_file($options{map}) or
    die "Unable to open $options{map} for reading: $!";
my $ped = open_compressed_file($options{ped}) or
    die "Unable to open $options{ped} for reading: $!";

my $ped_out; my $map_out;
open($ped_out,'>:encoding(UTF-8)',$options{ped_out}) or
    die "Unable to open $options{ped_out} for writing: $!";
open($map_out,'>:encoding(UTF-8)',$options{map_out}) or
    die "Unable to open $options{map_out} for writing: $!";


# read VCF file
my $vcf = open_compressed_file($options{vcf}) or
    die "Unable to open $options{vcf} for reading: $!";
my $rsids = read_vcf_file($vcf);
close($vcf);

my $illumina = open_compressed_file($options{illumina_data}) or
    die "Unable to open $options{illumina_data} for reading: $!";
# read illumina SNPs in
my $snps = read_illumina($illumina);
close($illumina);

my $db = Bio::DB::Fasta->new($options{ref}) or
    die "Unable to open $options{ref} for reading: $!";
for my $snp (sort keys %{$snps}) {
    next if $snps->{$snp}{chr} eq "0"; # skip chr 0
    next if $snps->{$snp}{chr} eq "XY"; # this isn't a chr
    my $il_gen_seq = uc($snps->{$snp}{genom_seq});
    if (length($il_gen_seq) > 10) {
        $il_gen_seq =~ s/^.*(.{10})$/$1/;
    }
    my $seq;
    my $il_gen_seq_len = length($il_gen_seq);
    if ($snps->{$snp}{flip}) {
        $il_gen_seq = flip_strand($il_gen_seq);
        $seq = $db->seq($snps->{$snp}{chr},
                        $snps->{$snp}{pos}+1,
                        $snps->{$snp}{pos}+$il_gen_seq_len,
                       );
    } else {
        $seq = $db->seq($snps->{$snp}{chr},
                        $snps->{$snp}{pos}-$il_gen_seq_len,
                        $snps->{$snp}{pos}-1
                       );
    }
    $snps->{$snp}{ambiguous} = 0;
    if ($seq ne $il_gen_seq) {
        # OK, this looks like a case where illumina has the wrong
        # genomic sequence
        my $rev_il_gen_seq;
        my $rev_seq;
        if ($snps->{$snp}{flip}) {
            ## we flipped the strand previously; flip it back
            $rev_il_gen_seq = flip_strand($il_gen_seq);
            $rev_seq = $db->seq($snps->{$snp}{chr},
                                $snps->{$snp}{pos}-$il_gen_seq_len,
                                $snps->{$snp}{pos}-1
                               );
        } else {
            ## we haven't flipped it; flip it now.
            $rev_il_gen_seq = flip_strand($il_gen_seq);
            $rev_seq = $db->seq($snps->{$snp}{chr},
                                $snps->{$snp}{pos}+1,
                                $snps->{$snp}{pos}+$il_gen_seq_len,
                               );
        }
        if ($rev_seq eq $rev_il_gen_seq) {
            $snps->{$snp}{flip} = $snps->{$snp}{flip} ? 0 : 1;
            print STDERR "$snp was wrong; flipping now $snps->{$snp}{flip}\n";
        } else {
            my $rev_dist = distance($rev_seq,$rev_il_gen_seq);
            my $dist = distance($seq,$il_gen_seq);
            if ($rev_dist < $dist and $rev_dist < $il_gen_seq_len * 0.4) {
                $snps->{$snp}{flip} = $snps->{$snp}{flip} ? 0 : 1;
                print STDERR "$snp was wrong; flipping now $snps->{$snp}{flip}\n";
            }elsif ($dist < $rev_dist and $dist < $il_gen_seq_len * 0.4) {
                # do nothing; the forward is the best match
            } else {
                # use Data::Printer;
                # p $snps->{$snp};
                # print STDERR "distance (rev): ".distance($rev_seq,$rev_il_gen_seq)."\n";
                # print STDERR "distance: ".distance($seq,$il_gen_seq)."\n";
                print STDERR "Can't match $rev_seq ($seq) with $rev_il_gen_seq ($il_gen_seq) at $snp\n";
                # compare alleles to VCF, flip if possible, otherwise complain
                if (exists $rsids->{$snp}) {
                    my $ref=$rsids->{$snp}{ref};
                    my $alt=$rsids->{$snp}{alt};
                    my $alleles=$snps->{$snp}{alleles};
                    $alleles=~s/[\[\]\/]//g;
                    if ($alleles eq flip_strand($alleles)) {
                        print STDERR "rsid $snp is also ambiguous, should be removed\n";
                        $snps->{$snp}{ambiguous} = 1;
                    } else {
                        print STDERR "illumina alleles: $alleles ref: $ref alt: $alt ";
                        if ($alleles eq "$ref$alt" or
                            $alleles eq "$alt$ref"
                           ) {
                            print STDERR "flip 0\n";
                            $snps->{$snp}{flip} = 0;
                        } else {
                            print STDERR "flip 1\n";
                            $snps->{$snp}{flip} = 1;
                        }
                    }
                }

            }
        }
    }
}

# read map file
my $map_rsids = read_map_file($map);
close($map);

# read merge file
my $merge = open_compressed_file($options{merge}) or
    die "Unable to open $options{merge} for reading: $!";
my $merged_rsids =
    read_merge_file($merge,
                   {(map {($_->{rsid},1)}
                     @{$map_rsids}),
                     (map {($_,1)} keys %{$snps})});
close($merge);

for my $m_rsid (keys %{$merged_rsids}) {
    $snps->{$merged_rsids->{$m_rsid}} =
        $snps->{$m_rsid} unless exists
        $snps->{$merged_rsids->{$m_rsid}};
}

# write map file
write_map_file($map_out,$map_rsids);
close($map_out);

my %reported_rsids;
while (<$ped>) {
    chomp;
    # The fields in a PED file are
    # 
    # Family ID
    # Sample ID
    # Paternal ID
    # Maternal ID
    # Sex (1=male; 2=female; other=unknown)
    # Affection (0=unknown; 1=unaffected; 2=affected)
    # Genotypes (space or tab separated, 2 for each marker. 0=missing)
    my ($fam,$sam,$pat,$mat,$sex,$aff,@genotypes) = split /\s+/;
    print {$ped_out} join(" ",$fam,$sam,$pat,$mat,$sex,$aff);
    for my $i (0..(@genotypes/2-1)) {
        my ($a,$b) = ($genotypes[$i*2],
                      $genotypes[$i*2+1]);
        my $mr = $map_rsids->[$i];
        next if $mr->{duplicate};
        if (not exists $snps->{$mr->{rsid}}) {
            print STDERR "Don't know about $mr->{rsid}";
        }
        if (exists $snps->{$mr->{rsid}} and
            $snps->{$mr->{rsid}}{flip}) {
            $a = flip_strand($a);
            $b = flip_strand($b);
        }
        if ((defined $snps->{$mr->{rsid}}{ambiguous} and
             $snps->{$mr->{rsid}}{ambiguous}) or
            $snps->{$mr->{rsid}}{chr} eq "0" or
            $snps->{$mr->{rsid}}{chr} eq "XY"
           ) {
            $a = "0";
            $b = "0";
        }
        if (not $reported_rsids{$mr->{rsid}} and
            $a ne 0 and $b ne 0 and
            exists $rsids->{$mr->{rsid}} and
            not $rsids->{$mr->{rsid}}{multiallelic} and
            not (($rsids->{$mr->{rsid}}{ref} eq $a or
                  $rsids->{$mr->{rsid}}{alt} eq $a) and
                 ($rsids->{$mr->{rsid}}{ref} eq $b or
                  $rsids->{$mr->{rsid}}{alt} eq $b))
           ) {
            $reported_rsids{$mr->{rsid}} = 1;
            use Data::Printer;
            p $snps->{$mr->{rsid}};
            p $rsids->{$mr->{rsid}};
            print STDERR map{"$_\n"} ($a,$b,@genotypes[($i*2,$i*2+1)]);
            print STDERR "Unable to fix strandedness at $mr->{rsid} ($.)".
                "[$snps->{$mr->{rsid}}{chr}:$snps->{$mr->{rsid}}{pos}] $a/$b ".
                "vs $rsids->{$mr->{rsid}}{ref}/$rsids->{$mr->{rsid}}{alt}\n";
            if ($. == 1) { # if this is the first sample, we can flip it back
                my $cur_flip = $snps->{$mr->{rsid}}{flip} // 0;
                $snps->{$mr->{rsid}}{flip} = ! $cur_flip;
                print STDERR "Changing flip from $cur_flip to $snps->{$mr->{rsid}}{flip}\n";
                $a = flip_strand($a);
                $b = flip_strand($b);
            }
        }
        print {$ped_out}
            "\t".$a." ".$b;
    }
    print {$ped_out} "\n";
}


sub write_map_file {
    my ($map_out,$map_rsids) = @_;

    for my $mr (@{$map_rsids}) {
        next if $mr->{duplicate};
        print {$map_out} join("\t",@{$mr}{qw(chr rsid recomb pos)})."\n";
    }
}

sub read_map_file {
    my ($map) = @_;
    my @map_rsids;
    my %duplicates;
    my %duplicates_rs;
    while (<$map>) {
        chomp;
        my ($chr,$rsid,$recomb,$pos) = split /\t/;
        my $duplicate = 0;
        if (exists $duplicates{$chr}{$pos} and
            defined $duplicates{$chr}{$pos} and
            $duplicates{$chr}{$pos}
           ) {
            $duplicate = 1;
        } elsif (exists $duplicates_rs{$rsid} and
                 defined $duplicates_rs{$rsid} and
                 $duplicates_rs{$rsid}
                ) {
            $duplicate = 1;
        }
        $duplicates{$chr}{$pos} = 1;
        $duplicates_rs{$rsid} = 1;
        push @map_rsids,
           {chr => $chr,
            pos => $pos,
            rsid => $rsid,
            recomb => 0, # $recomb,
            order => $#map_rsids,
            duplicate => $duplicate,
           };
    }
    return \@map_rsids;
}

sub flip_strand {
    my ($a) = @_;
    $a = reverse split //, $a;
    $a =~ tr/ATGC/TACG/;
    return $a;
}

sub read_vcf_file {
    my ($vcf) = @_;
    my %rsids;
    while (<$vcf>) {
        if (/^#/o) {
            next;
        }
        my ($chrom,$pos,$id,$ref,$alt,@junk) = split "\t",$_;
        $rsids{$id} =
           {chr => $chrom,
            pos   => $pos,
            id    => $id,
            ref   => $ref,
            alt   => $alt,
            multiallelic => $alt =~ /,/?1:0,
           };
    }
    return \%rsids;
}

sub numeric_cmp {
    my ($a,$b) = @_;
    if (looks_like_number($a) and looks_like_number($b)) {
        return $a <=> $b;
    }
    return $a cmp $b;
}

sub read_merge_file {
    my ($merge,$map_rsids) = @_;
    my %merged_rsids;

    while (<$merge>) {
        chomp;
        my ($old,$new,undef) = split /\t/;
        next unless exists $map_rsids->{'rs'.$old};
        $merged_rsids{'rs'.$old} = 'rs'.$new;
    }
    close ($merge);
    return \%merged_rsids;
}


sub read_illumina {
    my ($fh) = @_;

    my %snps;
    my @header;
    my $in_assay;
    while (<$fh>) {
        # Skip until we read the assay section
        if (not $in_assay) {
            next unless /^\[Assay\]/;
            $in_assay=1;
            next;
        }
        last if /^\[/;
        s/\r?\n$//; # the files are usually in DOS mode
        my @r = split /,/,$_;
        # Illumina files have multiple sections; we only want the
        # Assay section.
        if (not @header) {
            @header = @r and next;
        }
        my %r;
        @r{@header} = @r;
        ## http://www.sciencedirect.com/science/article/pii/S0168952512000704

        ## I think I'm finally understanding this. Illumina is always
        ## outputting the TOP allele. Thus, we only need to strand
        ## flip if the illumina probe sequence is BOT and the refseq
        ## is + or TOP and -, respectively.

        ## http://gengen.openbioinformatics.org/en/latest/tutorial/coding/
        $r{flip} = 0;
        if ($r{IlmnStrand} eq 'TOP' and $r{RefStrand} eq '-') {
            $r{flip} = 1;
        } elsif ($r{IlmnStrand} eq 'BOT' and $r{RefStrand} eq '+') {
            $r{flip} = 1;
        }
        ## however, illumina sometimes doesn't actually correctly
        ## match the TOP to the reference, so we need to fix that
        ## later.
        $r{genom_seq} = $r{TopGenomicSeq};
        $r{genom_seq} =~ s/\[.+//;
        $snps{$r{Name}} =
           {strand => $r{IlmnStrand},
            source_strand => $r{SourceStrand},
            ref_strand => $r{RefStrand},
            flip => $r{flip},
            genom_seq => $r{genom_seq},
            chr => $r{Chr},
            pos => $r{MapInfo},
            name => $r{Name},
            info => $r{IlmnID},
            alleles => $r{SNP},
           };
    }
    return \%snps;
}

sub open_compressed_file {
    my ($file) = @_;
    my $fh;
    my $mode = '<:encoding(UTF-8)';
    my @opts;
    if ($file =~ /\.gz$/) {
        $mode = '-|:encoding(UTF-8)';
        push @opts,'gzip','-dc';
    }
    if ($file =~ /\.xz$/) {
        $mode = '-|:encoding(UTF-8)';
        push @opts,'xz','-dc';
    }
    if ($file =~ /\.bz2$/) {
        $mode = '-|:encoding(UTF-8)';
        push @opts,'bzip2','-dc';
    }
    open($fh,$mode,@opts,$file);
    return $fh;
}


__END__
