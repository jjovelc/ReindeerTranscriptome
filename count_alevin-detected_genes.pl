#!/usr/bin/perl

=head
This script reads a count matrix and a gene names files from
a alevin directory and extract the list of genes that have 
counts greater than zero.
=cut

use strict;
use warnings;

# Ensure that the correct number of arguments is provided
if (@ARGV != 1) {
    die "Usage: $0 <dir>\n";
}

my $dir = $ARGV[0];
my $outfile = $dir . "_genes_detected.txt";
open my $out_fh, ">", $outfile;
my $compressed_file = "$dir/alevin/quants_mat.mtx.gz";

if (-e $compressed_file) {
    system("gunzip", $compressed_file) == 0
        or die "Failed to gunzip $compressed_file: $!";
}

my $count_mtx = "$dir/alevin/quants_mat.mtx";
my $gene_list = "$dir/alevin/quants_mat_cols.txt";

# Open the count matrix file
open my $counts_fh, "<", $count_mtx or die "Could not open '$count_mtx': $!";

# Open the gene list file
open my $genes_fh, "<", $gene_list or die "Could not open '$gene_list': $!";

# Skip the first two lines of the count matrix file (header)
<$counts_fh>;
<$counts_fh>;

# Initialize a hash to store the summed abundances per gene
my %counts;

while (my $line = <$counts_fh>) {
    chomp $line;
    my ($cell, $gene, $abundance) = split(/\t/, $line);

    # Sum the abundance per gene
    $counts{$gene} += $abundance;
}

close $counts_fh;

# Initialize a hash to map gene names to indices
my %genes;
my $counter = 0;

while (my $gene = <$genes_fh>) {
    chomp $gene;
    $counter++;
    $genes{$counter} = $gene;
}

close $genes_fh;

# Print the gene name and its corresponding summed abundance
foreach my $key (sort keys %counts) {

    if (exists $genes{$key}){ 
        my $gene_name = $genes{$key};
        my $gene_count = $counts{$key};
        print $out_fh "$gene_name\n";
    } else {
        print "Gene: $key not found \n";
    }  
}

close $out_fh;

exit;
