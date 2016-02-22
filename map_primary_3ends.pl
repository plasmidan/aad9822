#!/bin/env perl -w 
use strict; 

##################################################     primary_3ends.pl      ####################################################################
######	Map primary 3' ends (predicted terminators) using term-seq data
######	run example: perl primary_3ends.pl AL009126.genes exp_path_file max_allowed_3utr min_coverage min_reps
#################################################################################################################################################

#Modules
use genes_scripts; 

#command line inputs
my ($genes_file,$exp_file,$max_3UTR,$min_cov,$min_reps) = @ARGV;
$max_3UTR = 150 if (!$max_3UTR); 
$min_cov = 4 if (!$min_cov); 
$min_reps = 3 if (!$min_reps);

#hard coded parameters
my $min_avg_cov = 4;
my $contig = "AL009126";
my $is_normalize = 1;

#### Collect genome sequence and genes data
$exp_file =~ m/(.+)\.exp/; my $output_file_name = $1; 
my $genome_str = &get_genome_sequence($contig);
my $genome_length = length($genome_str);
my @genes_array = @{genes_scripts->get_genes($genes_file,$genome_length)};

#experiment files
my ($file_paths_ref,$sample_names_ref,$exp_names_ref,$exp_names_and_count_href) = &parse_exp_file($exp_file,"novo");
my @file_paths = @{$file_paths_ref};

#collect term-seq data from novoalign output files
my %pillar_hash; 
my $c = 1;
foreach my $novo_file (@file_paths) {
	my $experiment = $exp_names_ref->[$c-1]; 
	open(NOVO,$novo_file) or die "can't open pillar file $novo_file\n";
	while(<NOVO>) { 
		chomp;
		my @line = split(/\t/);
		my ($read_seq,$map_type,$map_pos,$map_strand) = @line[2,4,8,9]; 
		next unless ($map_type eq 'U'); # Uniquely mapped only
		$map_strand = ($map_strand eq 'F') ? '+' : '-';
		$map_pos += length($read_seq)-1 if ($map_strand eq '-');
		if (!$pillar_hash{$map_pos}) { 
			$pillar_hash{$map_pos} = [(0) x ($#file_paths+2)]; 
			$pillar_hash{$map_pos}->[0] = $map_strand;
			$pillar_hash{$map_pos}->[$c]++;				
		}
		else { 
			$pillar_hash{$map_pos}->[$c]++;
		}
	}
	close NOVO;
	$c++;
}

open(OUT3P,'>',"./$output_file_name.3p") or die "can't write 3p file...\n";
open(TABLE,'>',"./$output_file_name.gene_table") or die "can't write gene table file...\n";
print TABLE "locus\tmax_pos\tstrand\tmax_cov\tlength_3UTR\tis_within_coding\tterm_seq\tterm_seq_dnstream\tfold\tenergy\n";

#map sites to genes and output to files
my @pillar_positions =  sort {$a<=>$b} keys %pillar_hash; 
foreach my $gene_obj (@genes_array) { 
	my ($loc,$fr,$to,$strand) = $gene_obj->get_basic_gene_data;
	my ($fr_3UTR,$to_3UTR,$span_3UTR) = (0,0,0);
	if ($strand eq '+') { 
		($fr_3UTR,$to_3UTR) = ($to,$to + $max_3UTR); 
	}
	else { 
		($fr_3UTR,$to_3UTR) = ($fr - $max_3UTR,$fr); 
	}
	my @pillars_in_bounds = grep {$_ >= $fr_3UTR and $_ <= $to_3UTR} @pillar_positions;	
	#filter sites
	@pillars_in_bounds = grep {&is_min_reps($pillar_hash{$_},$exp_names_and_count_href,$min_cov,$min_reps) and $pillar_hash{$_}->[0] ne $strand} @pillars_in_bounds;
	next unless (@pillars_in_bounds);
	### choose most covered site
	my ($max_pos,$max_cov,$length_3UTR) = &get_max_covered_pos(\@pillars_in_bounds,\%pillar_hash,$gene_obj,$exp_names_and_count_href,$min_reps); 
	next unless ($max_cov >= $min_avg_cov);
	my $is_within_coding = ($gene_obj->get_int_dnstream + 10 < $length_3UTR) ? 'within_coding' : 'intergenic'; # allow 10nt invasion
	if ($is_within_coding eq 'within_coding') {
		next unless ($max_cov >= $min_avg_cov + 10); # allow 10nt invasion into dnstream gene
	}	
	my $term_sequence = &get_seq($genome_str,$max_pos,$strand,40);
	my $term_seq_dnstream = &get_seq($genome_str,$max_pos,$strand,-20);
	chomp (my @rna_fold_output = `echo $term_sequence | RNAfold `);### fold sequence using RNAfold
	$rna_fold_output[1] =~ /^(.+)\s+\((.+)\)/;
	my ($fold,$energy) = ($1,$2);
	print TABLE "$loc\t$max_pos\t$strand\t$max_cov\t$length_3UTR\t$is_within_coding\t$term_sequence\t$term_seq_dnstream\t$fold\t$energy\n"; 
	print OUT3P "$contig\t$max_pos\t$pillar_hash{$max_pos}->[0]\t$max_cov\t-1\n";
}
close OUT3P;
close TABLE;


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

sub parse_exp_file {
	my ($exp_file_path,$type) = @_;
	$type = (!$type) ? 's2' : 'novo'; 
	my @exp_project_dirs;
	my @exp_run_dirs;
	my @sample_names;
	my @exp_names;
	my %exp_names_and_count_hash;
	my %current_exp_hash;
	my $first_exp;
	my $c=0;
	open(EXP,$exp_file_path) or die "can't open experiment file $exp_file_path\n";
	while (my $sample_line = <EXP>) { 
		chomp $sample_line; 
		my ($lib_type,$project_dir,$run_dir,$sample_name,$exp_name) = split(/\t/,$sample_line);
		next unless ($lib_type eq 'end3'); 
		die "wrong number of inputs in 3end line should be: [end3] [project_dir] [run_dir] [sample_name] [experiment_name]...\n" if (
		!$lib_type or !$project_dir or !$run_dir or !$sample_name or !$exp_name);
		$first_exp = $exp_name if (!$first_exp); 
		if ($exp_name ne $first_exp) { 
			$c++; 
			$first_exp = $exp_name;
		}
		push @exp_project_dirs,$project_dir;	
		push @exp_run_dirs,$run_dir;
		push @sample_names,$sample_name;
		push @exp_names,$exp_name;
		$exp_names_and_count_hash{$c}++;
	}

	my @file_paths;
	foreach (my $i=0; $i<=$#exp_names; $i++) { 
		my $exp_type = $exp_names[$i];
		my $file_path = "/home/labs/sorek/repos/projects/$exp_project_dirs[$i]/$exp_run_dirs[$i]/intermediate/$exp_run_dirs[$i].$type"; 
		push @file_paths,$file_path;
	}
	return (\@file_paths,\@sample_names,\@exp_names,\%exp_names_and_count_hash);
}


sub get_seq { 
	my ($genome_str,$pos,$strand,$span) = @_;
	my $seq; 
	if ($strand eq '+')	{ 
		if ($span > 0) {
			$seq = substr($genome_str,$pos-1-$span,$span);
		}
		else { 
			$span *= -1;
			$seq = substr($genome_str,$pos-1,$span);
		}
	}
	else { 
		if ($span > 0) {	
			$seq = substr($genome_str,$pos-1,$span);
		}
		else {
			$span *= -1;
			$seq = substr($genome_str,$pos-1-$span,$span);
		}
		$seq = reverse $seq; $seq =~ tr/ATGC/TACG/;
	};	
	$seq =~ tr/T/U/;
	return $seq;
}

sub average {
	my ($values_ref) = @_;
	my @values = @{$values_ref};
	my $count;
	my $total = 0;
	foreach (@values) { 
		next if ($_ eq 'NA' or $_ eq '-' or $_ eq '+');
		$count++;
		$total += $_;
	}

	return $count ? $total / $count : 0;
}

sub get_max_covered_pos {
	my ($array_ref,$hash_ref,$gene_obj,$exp_names_and_count_href,$min_cov,$min_reps) = @_;
	my @array = @{$array_ref};
	my $max_pos = $array[0];
	my $max_cov = &average($hash_ref->{$array[0]});
	foreach (@array) { 
		my $avg = &average($hash_ref->{$_});
		if ($avg >= $max_cov) { 
			$max_pos = $_;
			$max_cov = $avg;
		}
	}
	$max_cov = int $max_cov;
	my $length_3UTR = ($gene_obj->get_st eq '+') ? ($max_pos - $gene_obj->get_to + 1) : ($gene_obj->get_fr - $max_pos + 1); 
	
	my @experiment_num_array = sort {$a<=>$b} keys %{$exp_names_and_count_href}; 
	my @max_pos_array = @{$hash_ref->{$max_pos}}; shift @max_pos_array;
	my $pos_in_array = 0;
	my ($max_exp,$exp_max_cov) = (0,0);
	
	foreach my $exp_num (@experiment_num_array) { 
		my $num_of_samples = $exp_names_and_count_href->{$exp_num};
		my @array_slice = @max_pos_array[$pos_in_array..($pos_in_array+$num_of_samples-1)];
		my $exp_mean_cov = int &average(\@array_slice);
		$exp_max_cov = $exp_mean_cov if ($exp_mean_cov > $exp_max_cov); 
		$pos_in_array += $num_of_samples;
	}
	return ($max_pos,$exp_max_cov,$length_3UTR);
}

sub is_min_reps {
	my ($array_ref,$exp_names_and_count_href,$min_cov,$min_reps) = @_;
	my @array = @{$array_ref}; shift @array;
	my @experiment_num_array = sort {$a<=>$b} keys %{$exp_names_and_count_href};
	
	my $pos_in_array = 0;
	foreach my $exp_num (@experiment_num_array) {
		my $num_of_samples = $exp_names_and_count_href->{$exp_num};
		my $rep_count = 0;
		for (my $i=$pos_in_array;$i<=($pos_in_array+$num_of_samples-1);$i++) {
			$rep_count++ if ($array[$i] >= $min_cov);
		}
		if ($rep_count >= $min_reps) {
			return 1;
		}
		else { 
			return '';
		}
	}
}