package sequencing_data_class;
use genes_scripts;

sub new_seq_data_obj {
	my ($class) = @_; 
	my $self = {
				_5end => {
						_sense => {
							_5utr => {},
							_coding => {},
							_3utr => {}
						}, 
						_antisense => {
							_5utr => {},
							_coding => {},
							_3utr => {}						
						}
						},
				_3end => {} 
				};
	bless ($self,$class);
	return $self;	
}

sub feed_new_5end_data { 
	my ($self,$pill_obj,$data_type,$is_antisense,$gene_segment) = @_;	
	my $pill_pos = $pill_obj->get_pos;
	$self->{$data_type}->{$is_antisense}->{$gene_segment}->{$pill_pos} = $pill_obj;
	return $self;
}

sub get_data { 
	my ($self,$data_type,$is_antisense,$gene_segment) = @_; 
	my $data = $self->{$data_type}->{$is_antisense}->{$gene_segment};
	return $data;	
}

sub determine_TSS_and_5utr { 
	my ($self,$gene_obj,$options_ref) = @_; 
	my ($min_cov,$min_ratio,$is_use_intergenic) = ($options_ref->{TSS_min_cov},$options_ref->{min_TAPnoTAP_ratio},$options_ref->{use_intergenic});
	my $utr_5ends = $self->get_data('_5end','_sense','_5utr');
	$gene_obj->save_5utr_len($options_ref->{max_5utr_length});
	if (keys %{$utr_5ends}) {
		my @utr_5ends_positions = sort {$a<=>$b} keys %{$utr_5ends};
		my ($tss,$tss_pos);
		my ($tss_cov,$tss_ratio)  = (0,0);
		foreach $pill_pos (@utr_5ends_positions) {
			my $pill = $utr_5ends->{$pill_pos};
			my ($pill_cov,$pill_ratio) = ($pill->get_TAP,$pill->get_ratio);
			if ($pill_ratio >= $min_ratio and $pill_cov >= $min_cov) {
				if ($pill_cov > $tss_cov) {
					($tss,$tss_pos,$tss_cov,$tss_ratio) = ($pill,$pill_pos,$pill_cov,$pill_ratio);
				}
			}			
		}
		if ($tss) {
			if ($is_use_intergenic == 1) {
				my ($gene_st,$gene_fr,$gene_to,$upstream_int) = ($gene_obj->get_st,$gene_obj->get_fr,$gene_obj->get_to,$gene_obj->get_int_upstream);
				my $_5utr_len = ($gene_st eq '+') ? ($gene_fr - $tss_pos) : ($tss_pos - $gene_to); 
				if ($_5utr_len > abs($upstream_int)) {
					$min_cov = $min_cov + 10;
					if ($tss_cov > $min_cov) {
						$tss->mark_as_within_coding('tss_within_coding');
					}
					else { 
						return;
					}
				}
				else { 
					$tss->mark_as_within_coding('tss_intergenic');
				}
			}
			$gene_obj->save_tss_obj($tss);
			$gene_obj->save_5utr_len($_5utr_len);			
		}		
	}
}

sub is_legit_withicoding_tss { 
	my ($self,$gene_obj) = @_; 
	my $tss = $gene_obj->get_tss_obj;
	return $gene_obj if (!$tss);
	return $gene_obj if ($tss->is_within_coding eq 'tss_intergenic');
	my $tss_pos = $tss->get_pos;
	my $tss_cov = $tss->get_TAP;
	my $gene_data = $gene_obj->get_sequencing_data_obj; 
	my @gene_data = $gene_obj->get_basic_gene_data;
	my $upstream_gene_obj = $gene_obj->get_upstream_gene;
	my $upstream_gene_data = $upstream_gene_obj->get_sequencing_data_obj; 
	my $upstream_strand = ($gene_obj->get_st eq $upstream_gene_obj->get_st) ? '_sense' : '_antisense';
	my %upstream_pillars = %{$upstream_gene_data->get_data('_5end',$upstream_strand,'_coding')};
	my @upstream_gene_coding_positions = keys %upstream_pillars;
	return $gene_obj unless (@upstream_gene_coding_positions); 
	my ($local_sum,$number_of_local_pills,$local_avg,$ratio) = (0,0,0,0);
	foreach my $pill_pos (@upstream_gene_coding_positions) { 
		if ($pill_pos >= $tss_pos-250 and $pill_pos <= $tss_pos+250) { 
			$local_sum += ($upstream_pillars{$pill_pos})->get_TAP if (($upstream_pillars{$pill_pos})->get_st eq $tss->get_st);
			$number_of_local_pills++; 
		}
	}
	$local_avg = $local_sum / 500;
	$ratio = $tss_cov / $local_avg if ($local_avg > 0);
	return $gene_obj unless ($tss_cov > 30 and $ratio > 100);
	$gene_obj->coding_tss_legit;
	return $gene_obj;
}

sub feed_new_3end_data { 
	my ($self,$pill_obj,$exp,$sample,$is_antisense,$gene_segment) = @_;	
	my $pill_pos = $pill_obj->get_pos;
	my $data = $self->{_3end}->{$gene_segment}->{$is_antisense}->{$pill_pos}->{$exp}->{$sample} = $pill_obj->get_cov;
	return $self;
}

sub calc_all_exp_avg_cov { 
	my ($self) = @_;
	my @segments = keys %{$self->{_3end}};
	foreach $seg (@segments) { 
		my @sense_anti = keys %{$self->{_3end}->{$seg}};
		foreach $strand (@sense_anti) { 
			my @pillar_positions = keys %{$self->{_3end}->{$seg}->{$strand}};
			foreach my $pill_pos (@pillar_positions) { 
				my @experiments = keys %{$self->{_3end}->{$seg}->{$strand}->{$pill_pos}};
				foreach my $exp (@experiments) { 
					my @samples = values %{$self->{_3end}->{$seg}->{$strand}->{$pill_pos}->{$exp}};
					my ($reps,$sum,$avg) = (0,0,0);
					foreach my $sample_cov (@samples) { 
						if ($sample_cov > 0) { 
							$reps++;
							$sum+=$sample_cov;
						}
					}
					$avg = $sum / $reps unless ($reps == 0);
					$self->{_3end}->{$seg}->{$strand}->{$pill_pos}->{$exp}->{_avg_cov} = $avg;
					$self->{_3end}->{$seg}->{$strand}->{$pill_pos}->{$exp}->{_repeats} = $reps;			
				}	
			}
		}	
	}
	return $self;
}	

sub get_exp_avg_cov_and_rep { 
	my ($self,$seg,$is_antisense,$pill_pos,$exp) = @_;
	my $avg_cov = $self->{_3end}->{$seg}->{$is_antisense}->{$pill_pos}->{$exp}->{_avg_cov}; 
	my $rep = $self->{_3end}->{$seg}->{$is_antisense}->{$pill_pos}->{$exp}->{_repeats}; 
	return ($avg_cov,$rep);
}

sub collect_r5UTR_candidates {
	my ($self,$gene_obj,$options_ref) = @_;
	my ($gene_loc,$gene_fr,$gene_to,$gene_st) = $gene_obj->get_basic_gene_data; my $gene_desc = $gene_obj->get_desc;
	return $gene_obj if (!$self->{_3end}->{_5utr}->{_sense});
	my %pillar_positions = %{$self->{_3end}->{_5utr}->{_sense}};
	my $tss = $gene_obj->get_tss_obj;
	return $gene_obj if (!$tss and $options_ref->{is_tss_only_5utrs} == 1);
	$tss_pos = (!$tss) ? 'NA' : ($tss->get_pos);
	my $is_within_coding_tss = ($tss) ? $tss->is_within_coding : 'Unknown_TSS';
	if ($is_within_coding_tss eq 'tss_within_coding') { 
		return $gene_obj if (!$gene_obj->is_coding_tss_legit);
	}
	my $min_element_len = $options_ref->{min_r5utr_len};
	return $gene_obj if (!($gene_obj->get_experiments));
	my @experiments = @{$gene_obj->get_experiments};
	my %dom_3end_per_exp;
	foreach my $exp (@experiments) { $dom_3end_per_exp{$exp}->{coverage} = 0; $dom_3end_per_exp{$exp}->{length} = 0; };
	my $is_dom = 0;
	PILLARS: foreach my $pill_pos (sort {$a<=>$b} keys %pillar_positions) { 
		my $r5utr_len = 0;
		if ($tss) { 
			$r5utr_len = ($gene_st eq '+') ? ($pill_pos - $tss_pos + 1) : ($tss_pos - $pill_pos + 1);
		}
		elsif (!$tss and $options_ref->{is_tss_only_5utrs} == 0) {
			$upstream_int_len = abs($gene_obj->get_int_upstream); 
			$r5utr_len = ($gene_st eq '+') ? ($pill_pos - $gene_fr - $upstream_int_len) : ($gene_to + $upstream_int_len - $pill_pos);
		}
		else { 
			next PILLARS;
		}
		next PILLARS unless ($r5utr_len >= $min_element_len);
		my @exps = keys %{$self->{_3end}->{_5utr}->{_sense}->{$pill_pos}};		
		foreach my $exp (@exps) { 
			if ($self->{_3end}->{_5utr}->{_sense}->{$pill_pos}->{$exp}) { 
				my ($avg,$rep) = $self->get_exp_avg_cov_and_rep('_5utr','_sense',$pill_pos,$exp);
				next PILLARS unless ($avg >= $options_ref->{'3p_min_avg_cov'} and $rep >= $options_ref->{'3p_min_rep'}); 
				if (
					($avg > $dom_3end_per_exp{$exp}->{coverage}) or 
					($avg == $dom_3end_per_exp{$exp}->{coverage} and $r5utr_len > $dom_3end_per_exp{$exp}->{length}) ) {
						$dom_3end_per_exp{$exp}->{coverage} = $avg;
						$dom_3end_per_exp{$exp}->{length} = $r5utr_len;
						$dom_3end_per_exp{$exp}->{pos} = $pill_pos;					
						$dom_3end_per_exp{$exp}->{rep} = $rep;	
						$is_dom=1;
				}
			}
		}
	}
	return $gene_obj unless ($is_dom > 0 );
	$gene_obj->save_r5utr_info(\%dom_3end_per_exp);
	return $gene_obj;
}


sub set_multiple_map { 
	my ($self,$pill_pos,$seg,$is_antisense,$previous_gene_obj) = @_; 
	$self->{_3end}->{$seg}->{$is_antisense}->{$pill_pos}->{_multiple_mapping} = $previous_gene_obj;
}

sub get_3end_experiments { 
	my ($self) = @_; 
	my @data = keys %{$self->{_3end}};
	return \@data;	
}

sub get_3end_samples { 
	my ($self,$exp) = @_; 
	my @data = keys %{$self->{_3end}->{$exp}};
	return \@data;	
}

1	