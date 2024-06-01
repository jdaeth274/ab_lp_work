#!/usr/bin/perl

use strict;
use Bio::SeqIO;

my $usage ="perl filter_contigs.pl <min_length> <input_file> <output_file> \n";

# Getting the input variables
my $min_len = shift or die $usage;
my $input_seq = shift or die $usage;
my $output_seq = shift or die $usage;

#print "This is the minimum length: $min_len \n";
# Reading the input fasta file
my $seqio_in = Bio::SeqIO->new(-file => $input_seq, 
                             -format => "fasta" );

# Creating the output fasta file                             
my $seqio_out = Bio::SeqIO->new(-file => ">$output_seq", 
                             -format => "fasta" );

# Saving sequences to the output if length >= min_len     
while ( my $seq = $seqio_in->next_seq ) {
    if ( $seq->length  >=  $min_len ) {
        $seqio_out->write_seq($seq);
    }
}