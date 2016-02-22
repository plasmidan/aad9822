#!/bin/env perl -w 
use strict;

#######################################   calculate_stemloop_polyU_from_gene_table.pl  ##########################################################
######	Uses the mapped primary 3' ends table and measures the polyU length and whether there is a stem-loop fold
######	run example: perl calculate_stemloop_polyU_from_gene_table.pl termseq_genes_table
#################################################################################################################################################

#command line inputs
my ($term_table) = @ARGV;

open(TERM,$term_table) or die "missing term table $term_table\n";
chomp (my $header = <TERM>);
$header .= "\t#ofUs\t" . "stem_loop\n";
print $header;
while (<TERM>) { 
	chomp; 
	my @line = split(/\t/); 
	my $table_line = $_;
	my ($up_seq,$fold) = @line[6,8];
	# next;
	my $num_of_Us = &count_U_stretch($up_seq);
	my $stem_loop = &is_stem_loop($fold);
	$table_line .= "\t$num_of_Us\t" . "$stem_loop\n";
	print $table_line;
}




###########################################################################################
############################# Additional Functions ########################################
###########################################################################################

sub is_stem_loop { 
	my ($fold) = @_;
	$fold =~ m/(\(+\.{0,2}\(+\.{0,2}\(+)(\.+)(\)+\.{0,2}\)+\.{0,2}\)+).{0,10}$/; 
	if ($1 and $2 and $3 ) { 
		my $stem_diff = abs(length($1) - length($3));
		my $loop_stem_diff_1 = length($2) / length($1);
		my $loop_stem_diff_2 = length($2) / length($3);
		if ($stem_diff <= 4 and $loop_stem_diff_1 < 2 and $loop_stem_diff_2 < 2) { 
			my $stem_loop = $1 . $2 . $3;
			return $stem_loop;
		}
	}
	return 'no_stemloop';
}

sub count_U_stretch { 
	my ($seq) = @_;
	my $tail = substr($seq,length($seq)-8,8); 
	my @nucs = split(//,$tail); 
	my $number_of_U_residues = 0;
	foreach my $nuc (@nucs) { 
		$number_of_U_residues++ if ($nuc eq 'U');
	}
	return $number_of_U_residues;
}



