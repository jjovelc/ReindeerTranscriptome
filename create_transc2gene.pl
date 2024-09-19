#!/usr/bin/perl

=head
	Newly-assembled transcripts will have IDs that do not conform
	with names of transcripts assgined by ENSEMBL. This script will
	create a pseudo-id. The user indicates the last ID of the reference
	transcriptome, in this example ENSBTAT00000087508 and pseudo-ids
	will start at the following number ENSBTAT00000087509. Although
	this is not necessary, make naming look better.

	An matching list of original and pseudo IDs will be save on file
	bt-rt_transcriptome_renaming.txt; changed name if desired/neede.

	The header of the FASTA file containing the newly-assembled reads should
	be in the following format:

	TRINITY_DN170886_c3_g1_i1_Zswim6 added protein:ZSWM6_MOUSE gene_symbol:Zswim6 description:Zinc finger SWIM domain-containing protein 6

	... from where the gene name will be parsed.

	Results will be printed to stdouput.

	Renamed transcripts will look like:

	ENSBTAT00000098082.1 added protein:ZNT5_MOUSE gene_symbol:Slc30a5 description:Proton-coupled zinc antiporter SLC30A5	
	
=cut


use warnings;
use strict;

my $infile = $ARGV[0];
my $outfile = "bt-rt_transcriptome_renaming.txt";

open my $in_fh, "<", $infile or die "Could not open file $infile: $!";
open my $out_fh, ">", $outfile or die "Could not open file $outfile: $!";

my $last_id = "ENSBTAT00000087508";

while (my $line = <$in_fh>){
    chomp $line;

    # Grab transcript ID 
    my ($transc) = $line =~ /^(\S+)/;

    # Check if transcript name does not match the desired pattern
    if ($transc !~ /^ENSBTAT/) {
        # Increment the numeric part of the last_id
        my ($prefix, $numeric_part) = $last_id =~ /^(\D+)(\d+)$/;
        $numeric_part++;

        # Generate the new transcript ID
        $last_id = sprintf("%s%011d.1", $prefix, $numeric_part);  # Add version .1

        # Replace transcript ID with the new one
        $transc = $last_id;
    }

    # Get gene symbol
    my ($gene_symbol) = $line =~ /gene_symbol:(\S+)/;
    $gene_symbol = ($line =~ /gene:(\S+)/)[0] unless $gene_symbol;

    # Print the new transcript ID along with the gene symbol
    print "$transc\t$gene_symbol\n";
}


close $in_fh;
exit;
