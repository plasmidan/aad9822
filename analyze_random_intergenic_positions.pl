#!/bin/env perl -w 
use strict;

####################################### random_site_terminator_analysis.pl #################################
####	generate a list of random positions in a given genome and folds the 40bp upstream
####    Measures polyU and stem-loop length and stability
####	run example: perl random_site_terminator_analysis.pl genes_file
############################################################################################################

#Modules
use genes_scripts; 

#command line inputs
my ($contig,$genes_file,$number_of_sites) = @ARGV;
$contig = "AL009126" unless ($contig); 
$number_of_sites = 50 unless($number_of_sites);

#hard coded parameters
my $upstream_span = 40;

#### Collect genome sequence and genes data
my $genome_str = get_genome_sequence($contig);
my $genome_length = length($genome_str);

#### collect the intergenic regions
my @genes_array = @{genes_scripts->get_genes($genes_file,$genome_length)};
my @intergenic_genome; 
my $total_int_lengh = 0;
foreach my $gene_ref (@genes_array) { 
	my ($gene_to,$dn_int_len) = ($gene_ref->get_to,$gene_ref->get_int_dnstream);
	for (my $i=$gene_to;$i<=$gene_to+$dn_int_len;$i++) { 
		push @intergenic_genome,$i;
		$total_int_lengh++;
	}
}

open(OUT,'>',"./random_sites_${number_of_sites}_sites") or die "can't write output file\n";
print OUT "site\tstrand\tseq\tfold\tenergy\t#ofUs\tis_stem_loop?\n";
for(my $i=0;$i<$number_of_sites;$i++) { 
	my $rand_site = @intergenic_genome[int(rand($#intergenic_genome))]; #select intergenic position randomly
	my ($site_seq,$strand);	
	if (rand(1) > 0.5) { #strand of random site
		$strand = '+';
		$rand_site+=$upstream_span if ($rand_site-$upstream_span < 0);
		$site_seq = substr($genome_str,$rand_site-$upstream_span,$upstream_span); 
	}
	else {
		$strand = '-';	
		$rand_site-=$upstream_span if ($rand_site+$upstream_span > $genome_length);
		$site_seq = substr($genome_str,$rand_site,$upstream_span); 
		$site_seq = reverse $site_seq; 
		$site_seq =~ tr/ATGC/TACG/;
	}	
	my @rna_fold_output = `echo $site_seq | RNAfold -p`; #take centroid fold
	my $rna = $rna_fold_output[0]; chomp $rna;
	$rna_fold_output[1] =~ /^(.+)\s+\((.+)\)/;
	my ($fold,$energy) = ($1,$2); 
	$site_seq =~ tr/T/U/;
	
	my $num_of_Us = &count_U_stretch($site_seq);
	my $stem_loop = &is_stem_loop($fold,$energy);
	
	print OUT "$rand_site\t$strand\t$rna\t$fold\t$energy\t$num_of_Us\t$stem_loop\n";
}
close OUT;

###########################################################################################
############################# Additional Functions ########################################
###########################################################################################

sub get_genome_sequence { 
	my ($contig) = @_; 
	open(FA,"/home/labs/sorek/repos/genomes/$contig.fasta") or die "Fasta file missing...\n";
	my $genome_str; 
	<FA>;
	while (<FA>) {
		chomp;
		$genome_str .= $_; 
	}
	return $genome_str;
}	

sub is_stem_loop { 
	my ($fold,$energy) = @_;
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


