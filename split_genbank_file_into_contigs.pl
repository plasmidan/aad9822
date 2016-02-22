#!/bin/env perl -w 
use strict;

###############################################   split_genbank_file_into_contigs_metatranscriptome.pl   ###########################################
######	genbank file with multiple contigs into separate genbanks for each contig
######	run example: split_genbank_file_into_contigs_metatranscriptome.pl [align_file_dir] [genome_fasta] [output_dir] [is_uniq_only?]
######################################################################################################################################################

#command line inputs
my ($gbk_file,$output_dir) = @ARGV;

open(GB,$gbk_file) or die "gbk\n";
my $contig_name;
my $contig_lines; 
while(<GB>) { 
	if (/^LOCUS/) {
		my @line = split(/\s+/);
		$contig_name = $line[1]; 
	}
	if (/^\/\//) {
		open(OUT,'>',"$output_dir/$contig_name.gbk") or die "can;t write contig $contig_name\n";
		print OUT $contig_lines;
		undef $contig_lines;
		undef $contig_name;
	}
	else { 
		$contig_lines .= $_;
	}
}