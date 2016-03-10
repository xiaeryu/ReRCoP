#!/usr/bin/perl

######################################################
## This script takes in the contigs in multiple fasta file format and concatenates the contigs to form one long sequence.
## Input: 1> contig fasta file 2> header line of the output fasta file 3> output fasta file
## Output file: fasta file
######################################################

if(!defined($ARGV[2])){
	die "Usage: perl $0 <contig fasta file> <header of the output fasta file> <output fasta file>\n";
}

my $header=$ARGV[1];

open(OUT,"> $ARGV[2]") or die "Cannot write to the output fasta file:$!\n";
print OUT ">$header\n";
open(INPUT,$ARGV[0]) or die "Cannot open input contig fasta file:$!\n";
while(<INPUT>){
	chomp;
	unless(/^>/){
		print OUT "$_";
	}
}
print OUT "\n";
close INPUT;
close OUT;
