#!bin/env perl -w 
use strict; 

#######################################   prepare_terminator_polyU_vs_energy_mat.pl  ##########################################################
######	Prepare a matrix for the joint distribution of stem-loop stability and polyU length
######	run example: perl prepare_terminator_polyU_vs_energy_mat.pl termseq_genes_table
#################################################################################################################################################

#command line inputs
my ($term_table,$is_random) = @ARGV;

my @pU_groups = (0,1,2,3,4,5,6,7,8);
my @stability_groups; 
my $stab_array_len = 0;
for( my $i = 0; $i >= -20; $i-=1) { 
	push @stability_groups,$i;
	$stab_array_len++;
}

open(TERM,$term_table) or die "can't open table $term_table\n";
<TERM>;
my $total_num_of_terms = 0;

my %output_hash; 
foreach (@pU_groups) { 	
	$output_hash{$_} = [(0) x $stab_array_len];
}
my $site;
SITES:
while (<TERM>) { 
	$total_num_of_terms++;
	my @line = split(/\t/); 
	my ($energy,$pU_len); 
	if (!$is_random) { 
		($energy,$pU_len) = @line[9,10];
		$site = $line[1];
	}
	else { 
		($energy,$pU_len) = @line[4,5];
	}
	for (my $i = 0;$i < $stab_array_len; $i++) { 	
		if ($i+1 == $stab_array_len) {
			$output_hash{$pU_len}->[$i]++;
			next SITES;
		}
		if ($energy <= $stability_groups[$i] and $energy >= $stability_groups[$i+1]) {
			$output_hash{$pU_len}->[$i]++;
			next SITES;
		}
	}
}

#report table
print "polyU_length";
foreach (@stability_groups) { my $field = "dG_${_}"; print "\t$field";}; print "\n"; 
foreach my $pU_len (@pU_groups) { 
	my @stab_counts = @{$output_hash{$pU_len}};
	print "pU_${pU_len}"; 
	foreach my $stab_grp_count (@stab_counts) { 
		my $normed = $stab_grp_count / $total_num_of_terms;
		print "\t$normed";
	}
	print "\n";
}

