package end3_scripts;

sub new_3end { 
	my ($class,$pos,$strand,$coverage,$experiment,$sample) = @_; 
	my $self = {		
		_pos => $pos,
		_st => $strand,
		_cov => $coverage, 
		_normalized_cov => 0, 
		_type => '3end',
		_exp => $experiment,
		_sample => $sample
	};
	bless ($self,$class);
	return $self;
}

sub set_exp_coverage { 
	my ($self,$exp_type,$cov) = @_;
	$self->{"_$exp_type"} = $cov;
}

sub get_type { $_[0]->{_type}};
sub get_pos { $_[0]->{_pos}};
sub get_st { $_[0]->{_st}};
sub get_cov { $_[0]->{_cov}};
sub get_exp { $_[0]->{_exp}};
sub get_sample { $_[0]->{_sample}};


sub get_basic_pill_data { 
	my ($self) =@_; 
	my $pos = $self->get_pos;
	my $strand = $self->get_st;
	my $cov = $self->get_cov;
	return ($pos,$strand,$cov);
}

sub parse_exp_file {	
	my ($self,$exp_file_path,$type) = @_; 	
	$type = (!$type) ? 's2' : 'novo'; 
	my @exp_project_dirs;
	my @exp_run_dirs;
	my @sample_names;
	my @exp_names;
	my %exp_names_and_count_hash;
	my $c=0;
	my %current_exp_hash;
	my $first_exp;
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
	my @s2_file_paths;
	foreach ($i=0; $i<=$#exp_names; $i++) { 
		my $exp_type = $exp_names[$i];
		my $s2_path = "/home/labs/sorek/repos/projects/$exp_project_dirs[$i]/$exp_run_dirs[$i]/intermediate/$exp_run_dirs[$i].$type"; 
		`sort -k3,3 -n $s2_path > $s2_path.sorted` if ($type eq 's2');
		push @s2_file_paths,$s2_path;
	}
	return (\@s2_file_paths,\@sample_names,\@exp_names,\%exp_names_and_count_hash);
}


sub map_pillars_to_genes { 
	my ($self,$end3_data_hash_ref,$genes_array_ref,$options_ref) = @_;
	my ($min_cov,$min_rep,$max_5utr_len,$max_3utr_len) = (
	$options_ref->{'3p_min_avg_cov'},$options_ref->{'3p_min_rep'},$options_ref->{max_5utr_length},$options_ref->{max_3utr_length} );
	my @genes_array = @{$genes_array_ref};
	my %end3_data_hash = %{$end3_data_hash_ref};
	my $multiple_map_hash_ref;
	my @experiments = keys %end3_data_hash;
	foreach my $exp (@experiments) { 
		my @samples_of_experiment = keys $end3_data_hash{$exp};
		foreach my $sample (@samples_of_experiment) { 
			my @pillar_positions = sort {$a <=> $b} keys $end3_data_hash{$exp}->{$sample}; 
			my $gene_idx = 0;
			PILLARS: for (my $pillar_idx=0;$pillar_idx<= $#pillar_positions;$pillar_idx++) { 
				my $pill_obj = $end3_data_hash{$exp}->{$sample}->{$pillar_positions[$pillar_idx]};
				my ($pill_pos,$pill_st,$pill_cov) = $pill_obj->get_basic_pill_data;	
				my $gene_obj = $genes_array[$gene_idx];
				my $next_gene_obj; 
				$next_gene_obj = $genes_array[$gene_idx+1];
				$gene_obj->save_experiments(\@experiments);				
				my $return_message;
				$return_message = end3_scripts->map_3end_pill_to_gene($pill_obj,$gene_obj,$exp,$sample,$options_ref);				
				if ($return_message eq 'next_gene') { 
					$gene_idx++;
					redo PILLARS;
				}
				elsif ($return_message eq 'next_pill') { 
					next PILLARS; 
				}
				elsif ($return_message eq 'compare2nextGene' and $next_gene_obj) {
					($return_message,$multiple_map_hash_ref) = end3_scripts->map_3end_pill_to_gene($pill_obj,$next_gene_obj,$exp,$sample,$options_ref,$gene_obj,$multiple_map_hash_ref);
					next PILLARS; 
				}
			}
		}
	}
	return $multiple_map_hash_ref;
}

sub map_3end_pill_to_gene { 
	my ($self,$pill_obj,$gene_obj,$exp,$sample,$options_ref,$previous_gene_obj,$multiple_map_hash_ref) = @_; 
	my ($pill_pos,$pill_strand,$pill_cov) = $pill_obj->get_basic_pill_data;
	my ($gene_loc,$gene_fr,$gene_to,$gene_strand) = $gene_obj->get_basic_gene_data;
	my $tss = $gene_obj->get_tss_obj;
	my ($min_cov,$min_rep,$max_5utr_len,$max_3utr_len) = (
	$options_ref->{'3p_min_avg_cov'},$options_ref->{'3p_min_rep'},$options_ref->{max_5utr_length},$options_ref->{max_3utr_length} );
	my $return_message;
	my $is_antisense = ($pill_strand ne $gene_strand) ? '_sense' : '_antisense';
	if ($gene_strand eq '+') { 
		if ($pill_pos < $gene_fr - $max_5utr_len) {
			$return_message = 'next_pill';
		}
		elsif ($pill_pos >= $gene_fr - $max_5utr_len and $pill_pos <= $gene_fr) { 
				$multiple_map_hash_ref->{$pill_pos} = 1 if ($previous_gene_obj and $options_ref->{multiple_mapping} == 1);
				$gene_obj->feed_gene_seq_3end_data($pill_obj,$exp,$sample,$is_antisense,'_5utr');
				$return_message = 'compare2nextGene';
		}
		elsif ( $pill_pos > $gene_fr and $pill_pos <= $gene_to ) { 
			$multiple_map_hash_ref->{$pill_pos} = 1 if ($previous_gene_obj and $options_ref->{multiple_mapping} == 1);
			$gene_obj->feed_gene_seq_3end_data($pill_obj,$exp,$sample,$is_antisense,'_coding');
			$return_message	= 'compare2nextGene';
		}
		elsif ( $pill_pos > $gene_to and $pill_pos < $gene_to + $max_3utr_len ) { 
			$multiple_map_hash_ref->{$pill_pos} = 1 if ($previous_gene_obj and $options_ref->{multiple_mapping} == 1);
			$gene_obj->feed_gene_seq_3end_data($pill_obj,$exp,$sample,$is_antisense,'_3utr');
			$return_message	= 'compare2nextGene';	
		}
		else { 
			$return_message = 'next_gene';
		}
	}
	else { 
		if ($pill_pos < $gene_fr - $max_3utr_len) { 
			$return_message = 'next_pill';
		}
		elsif ($pill_pos < $gene_fr and $pill_pos >= $gene_fr - $max_3utr_len) {
			$multiple_map_hash_ref->{$pill_pos} = 1 if ($previous_gene_obj and $options_ref->{multiple_mapping} == 1);
			$gene_obj->feed_gene_seq_3end_data($pill_obj,$exp,$sample,$is_antisense,'_3utr');
			$return_message	= 'compare2nextGene';
		}
		elsif ($pill_pos >= $gene_fr and $pill_pos < $gene_to) {
				$multiple_map_hash_ref->{$pill_pos} = 1 if ($previous_gene_obj and $options_ref->{multiple_mapping} == 1);
				$gene_obj->feed_gene_seq_3end_data($pill_obj,$exp,$sample,$is_antisense,'_coding');				
				$return_message = 'next_pill';
		}
		elsif ( $pill_pos >= $gene_to and $pill_pos <= $gene_to + $max_5utr_len ) { 
			$multiple_map_hash_ref->{$pill_pos} = 1 if ($previous_gene_obj and $options_ref->{multiple_mapping} == 1);
			$gene_obj->feed_gene_seq_3end_data($pill_obj,$exp,$sample,$is_antisense,'_5utr');
			$return_message	= 'compare2nextGene';
		}
		else  {
			$return_message	= 'next_gene';
		}
	}
	if ($multiple_map_hash_ref) { 
		return ($return_message,$multiple_map_hash_ref);
	}
	else {
		return $return_message; 
	}
}

sub collect_3end_data_into_hash {
	my ($self,$exp_file_path,$min_cov,$min_reps) = @_; 
	my ($s2_paths_ref,$sample_names_ref,$exp_names_ref) = end3_scripts->parse_exp_file($exp_file_path);
	my @s2_paths = @{$s2_paths_ref}; ;
	my @sample_names = @{$sample_names_ref};
	my @exp_names = @{$exp_names_ref};
	
	my %end3_data_hash;
	foreach ($i=0; $i<=$#exp_names; $i++) { 
		my $exp_type = 'end3';
		my ($s2,$exp,$sample) = ($s2_paths[$i],$exp_names[$i],$sample_names[$i]);;		
		open (S2,"$s2") or die "can't open S2 file in 3'end $s2\n";
		while (<S2>) { 
			my @line = split(/\t/); my ($pos,$coverage,$strand) = @line[2,4,6];			
				$end3_data_hash{$exp}->{$sample}->{$pos} = end3_scripts->new_3end($pos,$strand,$coverage,$exp,$sample); 
		}
		close S2;	
	}
	return \%end3_data_hash;
}

1