#!/usr/bin/perl
# remap_map_and_ped.pl remaps map files and ped files using a VCF
# and is released under the terms of the GNU GPL version 3, or any
# later version, at your option. See the file README and COPYING for
# more information.
# Copyright 2017 by Don Armstrong <don@donarmstrong.com>.


use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

=head1 NAME

remap_map_and_ped.pl - remaps map files and ped files using a VCF

=head1 SYNOPSIS

remap_map_and_ped.pl [options]

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

remap_map_and_ped.pl

=cut

use Scalar::Util qw(looks_like_number);

use vars qw($DEBUG);

my %options = (debug           => 0,
               help            => 0,
               man             => 0,
              );

GetOptions(\%options,
           'vcf=s',
           'ped=s',
           'map=s',
           'merge=s',
           'ped_out|ped-out=s',
           'map_out|map-out=s',
           
           'debug|d+','help|h|?','man|m');

pod2usage() if $options{help};
pod2usage({verbose=>2}) if $options{man};

$DEBUG = $options{debug};

my @USAGE_ERRORS;
for (qw(vcf ped map ped_out map_out)) {
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

my $merge = open_compressed_file($options{merge}) or
    die "Unable to open $options{merge} for reading: $!";

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

# read map file
my $map_rsids = read_map_file($map);
close($map);

# read merge file
my $merged_rsids = read_merge_file($merge,{map {($_->{rsid},1)} @{$map_rsids}});
close($merge);

# merge snps
for my $mr (@{$map_rsids}) {
    if (exists $merged_rsids->{$mr->{rsid}}) {
        $mr->{rsid} = $merged_rsids->{$mr->{rsid}};
    }
}

# fix map chr and pos in map file
for my $mr (@{$map_rsids}) {
    if (exists $rsids->{$mr->{rsid}}) {
        $mr->{chr} = $rsids->{$mr->{rsid}}{chr};
        $mr->{pos} = $rsids->{$mr->{rsid}}{pos};
    }
}

# reorder map file
@{$map_rsids} =
    sort {numeric_cmp($a->{chr},$b->{chr}) ||
              $a->{pos} <=> $b->{pos} ||
              $a->{order} <=> $b->{order}
          } @{$map_rsids};

# write map file
write_map_file($map_out,$map_rsids);
close($map_out);

# read and reorder ped file

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
    for my $mr (@{$map_rsids}) {
        print {$ped_out}
            "\t".$genotypes[$mr->{order}*2]." ".$genotypes[$mr->{order}*2+1];
    }
    print {$ped_out} "\n";
}
close($ped_out);
close($ped);

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

sub write_map_file {
    my ($map_out,$map_rsids) = @_;

    for my $mr (@{$map_rsids}) {
        print {$map_out} join("\t",@{$mr}{qw(chr rsid recomb pos)})."\n";
    }
}

sub numeric_cmp {
    my ($a,$b) = @_;
    if (looks_like_number($a) and looks_like_number($b)) {
        return $a <=> $b;
    }
    return $a cmp $b;
}

sub read_map_file {
    my ($map) = @_;
    my @map_rsids;
    while (<$map>) {
        chomp;
        my ($chr,$rsid,$recomb,$pos) = split /\t/;
        push @map_rsids,
           {chr => $chr,
            pos => $pos,
            rsid => $rsid,
            recomb => 0, # $recomb,
            order => scalar @map_rsids,
           };
    }
    return \@map_rsids;
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
           };
    }
    return \%rsids;
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
