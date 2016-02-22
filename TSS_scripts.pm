package TSS_scripts;

sub new_5end {
	my ($class,$pos,$strand) = @_; 
	my $self = {		
		_pos => $pos,
		_st => $strand,
		_TAP => 0, 
		_noTAP => 0, 
		_ratio => 0,
		_type => '5end'
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
sub get_TAP { $_[0]->{_TAP}};
sub get_noTAP { $_[0]->{_noTAP}};
sub get_ratio { $_[0]->{_ratio}};

sub get_basic_pill_data { 
	my ($self) =@_; 
	my $pos = $self->get_pos;
	my $strand = $self->get_st;
	my $tap = $self->get_TAP;
	my $notap = $self->get_noTAP;
	my $ratio = $self->get_ratio;
	return ($pos,$strand,$tap,$notap,$ratio);
}

sub set_ratio { $_[0]->{_ratio}=$_[1]};
sub mark_as_within_coding { $_[0]->{_within_coding}=$_[1]};
sub is_within_coding { $_[0]->{_within_coding};};

sub make_5end_object_hash {
	my ($self,$exp_file_path) = @_; 
	my ($end5_s2_paths_ref,$end5_exp_names_ref) = common_scripts->parse_exp_file($exp_file_path,'end5');
	my @experiment_project_dirs;
	my @end5_s2_paths = @{$end5_s2_paths_ref}; ;
	my @experiment_names = @{$end5_exp_names_ref};
	die "too many 5'end experiments - only one is currently allowed\n" if ($#experiment_names > 1);
	my %end5_position_hash; 	
	foreach ($i=0; $i<=$#experiment_names; $i++) { 
		print "experiment -> $experiment_names[$i]\n";
		my $exp_type = $experiment_names[$i];
		my $s2_path = $end5_s2_paths[$i];
		my $entry = 1; $entry = 2 if ($experiment_names[$i] eq 'noTAP');
		open (S2,"$s2_path") or die "can't open S2 file in 5'end $s2_path\n";
		while (<S2>) { 
			my @line = split(/\t/); my ($pos,$coverage,$strand) = @line[2,4,6];			
			if (!$end5_position_hash{$pos}) {
				$end5_position_hash{$pos} = TSS_scripts->new_5end($pos,$strand,$exp_type,$coverage);
				$end5_position_hash{$pos}->set_exp_coverage($exp_type,$coverage);
				
			}
			else { 
				$end5_position_hash{$pos}->set_exp_coverage($exp_type,$coverage);
			}	
		}
		close S2;	
	}
	while ( my ($pos,$pill_obj) = each %end5_position_hash ) {
		my ($strand,$tap,$notap,$ratio) = ($pill_obj->get_st,$pill_obj->get_TAP,$pill_obj->get_noTAP,$pill_obj->get_ratio);
		
		if ($notap == 0) { 
			$ratio = $tap;
		}
		else { 
			$ratio = $tap / $notap;
		}	
		$pill_obj->set_ratio($ratio);
		my $a = $pill_obj->get_ratio;
	}
	return \%end5_position_hash;
}

sub map_pillars_to_genes { 
	my ($self,$end5_position_hash_ref,$genes_array_ref,$min_cov,$min_5end_ratio,$max_5utr_len) = @_;
	my @genes_array = @{$genes_array_ref};
	my %end5_position_hash = %{$end5_position_hash_ref};
	my @pillar_positions = sort {$a <=> $b} keys %end5_position_hash; 
	my $gene_idx = 0;
	PILLARS: for (my $pillar_idx=0;$pillar_idx<= $#pillar_positions;$pillar_idx++) { 
		my $pill_obj = $end5_position_hash{$pillar_positions[$pillar_idx]};
		my ($pill_pos,$pill_st,$TAP,$noTAP,$ratio) = $pill_obj->get_basic_pill_data;
		my $gene_obj = $genes_array[$gene_idx];
		my $next_gene_obj; 
		$next_gene_obj = $genes_array[$gene_idx+1];
		my $return_message = TSS_scripts->map_pill_to_gene('5end',$pill_obj,$gene_obj,$max_5utr_len,$next_gene_obj);
		if ($return_message eq 'next_gene') { 
			$gene_idx++;
			redo PILLARS;
		}
		elsif ($return_message eq 'next_pill') { 
			next PILLARS; 
		}
		elsif ($return_message eq 'compare2nextGene' and $next_gene_obj) {
			TSS_scripts->map_pill_to_gene('5end',$pill_obj,$next_gene_obj,$max_5utr_len);
			next PILLARS; 
		}
	}
}

sub map_pill_to_gene { 
	my ($self,$type,$pill_obj,$gene_obj,$max_5utr_len,$next_gene_obj) = @_; 
	my ($pill_pos,$pill_strand,$pill_type) = ($pill_obj->get_pos,$pill_obj->get_st,$pill_obj->get_type);
	my ($gene_fr,$gene_to,$gene_strand) = ($gene_obj->get_fr,$gene_obj->get_to,$gene_obj->get_st);
	my $UTR_fr = ($gene_strand eq '+') ? ($gene_fr-$max_5utr_len) : ($gene_to+$max_5utr_len);
	my ($gene_dn_int,$gene_up_int) = ($gene_obj->get_int_dnstream,$gene_obj->get_int_upstream);
	my $return_message;
	my $is_antisense = ($pill_strand eq $gene_strand) ? '_sense' : '_antisense';
	
	if ($gene_strand eq '+') { 
		if ($pill_pos < $gene_fr - $max_5utr_len) {
			$return_message = 'next_pill';
		}
		elsif ($pill_pos >= $gene_fr - $max_5utr_len and $pill_pos <= $gene_fr) {
				$gene_obj->feed_gene_seq_data($pill_obj,'_5end',$is_antisense,'_5utr') if ($is_antisense eq '_sense');
				$return_message = 'next_pill';
		}
		elsif ( $pill_pos > $gene_fr and $pill_pos <= $gene_to ) { 
			$gene_obj->feed_gene_seq_data($pill_obj,'_5end',$is_antisense,'_coding');
			$return_message	= 'compare2nextGene';			
		}
		else { 
			$return_message = 'next_gene';
		}
	}
	else {
		if ($pill_pos < $gene_fr) { 
			$return_message = 'next_pill';
		}
		elsif ($pill_pos >= $gene_fr and $pill_pos < $gene_to ) { 
			$gene_obj->feed_gene_seq_data($pill_obj,'_5end',$is_antisense,'_coding');
			$return_message = 'compare2nextGene';
		}
		elsif ($pill_pos >= $gene_to and $pill_pos <= $gene_to + $max_5utr_len) { 
			my $next_gene_to = $pill_pos+10;
			if ($next_gene_obj) { 
				$next_gene_to = $next_gene_obj->get_to; 
			}			
			if ($next_gene_to > $pill_pos) {
				$gene_obj->feed_gene_seq_data($pill_obj,'_5end',$is_antisense,'_5utr');
			}	
			$return_message = 'compare2nextGene';
		}
		else { 
			$return_message = 'next_gene';
		}
	}
	return $return_message;
}

1 