use strict;
use warnings;

# Converts the annotation file into fasta and gtf files that can be added 
# to the end of an existing genome fasta/gtf. You may obtain the annotation file from:
# https://www.thermofisher.com/order/catalog/product/4456740

my @fasta_lines = ();
my @gtf_lines = ();
open (my $ifh, "annotation.txt") or die $!;
<$ifh>;  # skip the header line in tha annotation

while (<$ifh>) {
	chomp;  # do all the important stuff
	my @record = split(/\t/);
	my $sequence = $record[4];
	$sequence =~ s/\s+//g;  # get rid of any preceeding/tailing white space
	$sequence = $sequence."NNNN";  # add some buffer to the end of the sequence
	my $name = $record[0];
	my $genbank = $record[1];
	push(@fasta_lines, ">$name\n$sequence\n");
	
    # the sequences are all '+' stranded. and the gtf files are 1 indexed.
	push(@gtf_lines, "$name\tERCC\tgene\t1\t" . (length($sequence)-2) . "\t.\t+\t." .
	                 "\tgene_id \"$name-$genbank\"; transcript_id \"$name-$genbank\"; " .
	                 "exon_number \"1\"; gene_name \"ERCC $name-$genbank\"\n");
	push(@gtf_lines, "$name\tERCC\texon\t1\t" . (length($sequence)-2) . "\t.\t+\t." .
	                 "\tgene_id \"$name-$genbank\"; transcript_id \"$name-$genbank\"; " .
	                 "exon_number \"1\"; gene_name \"ERCC $name-$genbank\"\n");
} close($ifh);

# write output fasta and gtf files.
open(my $ofh, ">", "ercc.fasta") or die $!;
foreach my $line (@fasta_lines) {
	print $ofh $line;
} close ($ofh);

open($ofh, ">", "ercc.gtf") or die $!;
foreach my $line (@gtf_lines) {
	print $ofh $line;
} close ($ofh);