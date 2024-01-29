


#!/usr/bin/perl

use warnings;
use strict;

my $infile = $ARGV[0];

open my $in_fh, "<", $infile or die "Could not open file $infile: $!";

my $counter = 100;

while (my $line = <$in_fh>){
    chomp $line;

    # Skip the line if it doesn't start with '>'
    next unless $line =~ /^>/;

    my ($transc, @rest) = split / /, $line;
    $transc =~ s/>//;  # Remove '>' from the beginning of the ID

    my $rest_string = join ' ', @rest;

    if ($transc =~ /^ENSBTAT/){
        if ($rest_string =~ /gene:([^ ]+)/) {
            my $gene_ensbt = $1;
            if ($rest_string =~ /gene_symbol:([^ ]+)/) {
                my $gene_symbol = $1;
                print "$transc\t$gene_symbol\n";
            } else {
                print "$transc\t$gene_ensbt\n";
            }
        }
    } else {
        if ($rest_string =~ /protein:([^ ]+) gene_symbol:([^ ]+)/) {
            $counter++;
            my $protein = $1;
            $protein =~ s/_.*$//;
            my $gene_symbol = $2;
            print "$transc\t$protein\n";
        } else {
            print STDERR "Debug: Line not matching pattern: $line\n";
        }
    }
}

close $in_fh;
exit;
