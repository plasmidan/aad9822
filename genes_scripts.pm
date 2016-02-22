package genes_scripts; 

sub new_gene { 
	my ($class,@args) = @_;
	my $self = {
		_fr => $args[0],
		_to => $args[1],
		_st => $args[2],
		_desc => $args[3],
		_locus => $args[4],
		_name => $args[5],
	};
	bless ($self,$class);
	return $self;
}

sub get_genes { 
	my ($self,$genes_file,$genome_len) = @_;
	return if (!$genes_file or !$genome_len); 
	my @genes_features;
	if ($genes_file) { 
		my @genes_lines = `cat $genes_file | grep -v riboswitch`; 
		for (my $i = 0; $i<= $#genes_lines;$i++) { 
			my $line = $genes_lines[$i];
			my @line_sp = split(/\t/,$line); 
			my ($contig,$fr,$to,$st,$desc,$locus,$gene_name) = @line_sp[0..5,8];
			my $gene_obj = genes_scripts->new_gene($fr,$to,$st,$desc,$locus,$gene_name);
			push @genes_features,$gene_obj;
		}
	}
	
	for (my $i = 0;$i <= $#genes_features; $i++) {
		my $gene = $genes_features[$i];
		my $gene_st = $genes_features[$i]->get_st(); 
		my $gene_fr = $genes_features[$i]->get_fr();
		my $gene_to = $genes_features[$i]->get_to();		
		
		my ($dnstream_gene_fr,$dnstream_gene_to); 
		my ($upstream_gene_fr,$upstream_gene_to); 
		my ($downstream_dist,$upstream_dist);
		
		if ($i eq $#genes_features) {		
			($upstream_gene_fr,$upstream_gene_to) = ( $genes_features[$i-1]->get_fr(),$genes_features[$i-1]->get_to() );
			$upstream_dist = $gene_fr - $upstream_gene_to; 
			$downstream_dist = $genome_len - $gene_to + $genes_features[0]->get_int_upstream(); 
		}
		elsif ($i == 0)  {
			($dnstream_gene_fr,$dnstream_gene_to) = ( $genes_features[$i+1]->get_fr(),$genes_features[$i+1]->get_to() );
			$upstream_dist = $gene_fr;
			$downstream_dist = $dnstream_gene_fr - $gene_to;
		}
		else { 
			($upstream_gene_fr,$upstream_gene_to) = ( $genes_features[$i-1]->get_fr(),$genes_features[$i-1]->get_to() );
			($dnstream_gene_fr,$dnstream_gene_to) = ( $genes_features[$i+1]->get_fr(),$genes_features[$i+1]->get_to() );
			$upstream_dist = $gene_fr - $upstream_gene_to; 
			$downstream_dist = $dnstream_gene_fr - $gene_to;
		}
		$downstream_dist = $downstream_dist < 0 ? ($downstream_dist-1) : ($downstream_dist+1);  
		$upstream_dist = $upstream_dist < 0 ? ($upstream_dist-1) : ($upstream_dist+1);  
		($downstream_dist,$upstream_dist) = ($upstream_dist,$downstream_dist)if ($gene_st eq '-'); 
		$gene->save_int_upstream($upstream_dist);
		$gene->save_int_dnstream($downstream_dist);
	}
	my %genomic_feat; 
	foreach (@genes_features) { 
		$genomic_feat{$_->get_fr()} = $_;
	}
	my @keys = sort {$a <=> $b} keys (%genomic_feat);
	my @gene_objects;
	foreach (@keys) { 
		push @gene_objects,$genomic_feat{$_};
	}
	return \@gene_objects;
}

sub get_fr { $_[0]->{_fr} };
sub get_to { $_[0]->{_to} };
sub get_st { $_[0]->{_st} };
sub get_desc { $_[0]->{_desc} };
sub get_locus { $_[0]->{_locus} };
sub save_int_dnstream { $_[0]->{_int_dnstream} = $_[1]};
sub get_int_dnstream { $_[0]->{_int_dnstream} };
sub save_int_upstream { $_[0]->{_int_upstream} = $_[1]};
sub get_int_upstream { $_[0]->{_int_upstream} };
sub save_upstream_gene { $_[0]->{_upstream_gene} = $_[1]};
sub get_upstream_gene { $_[0]->{_upstream_gene} };
sub save_downstream_gene { $_[0]->{_downstream_gene} = $_[1]};
sub get_downstream_gene { $_[0]->{_downstream_gene} };
sub save_sequencing_data_obj { $_[0]->{_seq_data} = $_[1]};
sub get_sequencing_data_obj {$_[0]->{_seq_data}};
sub save_tss_obj { $_[0]->{_tss} = $_[1]};
sub get_tss_obj {$_[0]->{_tss}};
sub save_5utr_len { $_[0]->{_5utr_len} = $_[1]};
sub get_5utr_len {$_[0]->{_5utr_len}};
sub save_experiments { $_[0]->{_experiments} = $_[1]};
sub get_experiments {$_[0]->{_experiments}}
sub save_r5utr_info {$_[0]->{_r5utr_info} = $_[1]};
sub get_r5utr_info {$_[0]->{_r5utr_info}};
sub coding_tss_legit {$_[0]->{_coding_tss} = 1};
sub is_coding_tss_legit {$_[0]->{_coding_tss}};
sub save_known_r5UTR {$_[0]->{_known_r5UTR} = $_[1]};
sub get_known_r5UTR {$_[0]->{_known_r5UTR}};

sub report_seq_data { 
	my ($self,$data_type,$is_antisense,$gene_segment) = @_; 
	my @gene_data = $self->get_basic_gene_data;
	my $seq_obj = $self->get_sequencing_data_obj;
	my $pill_hash_ref = $seq_obj->get_data($data_type,$is_antisense,$gene_segment);
	if (keys %{$pill_hash_ref}) { 
		print "@gene_data\n>>>>>\n";
		my @keys = sort {$a<=>$b} keys %{$pill_hash_ref}; 
		print "@keys\n" if (@keys);
	}	
}

sub get_basic_gene_data { 
	my ($self) = @_;
	my $fr = $self->get_fr;
	my $to = $self->get_to;
	my $strand = $self->get_st;
	my $loc = $self->get_locus;
	my $name = $self->{_name};
	return ($loc,$fr,$to,$strand,$name);
}

sub feed_gene_seq_data {
	my ($self,$pill_obj,$data_type,$is_antisense,$gene_segment) = @_;
	my $data_obj = $self->get_sequencing_data_obj;
	$self->save_sequencing_data_obj( $data_obj->feed_new_5end_data($pill_obj,$data_type,$is_antisense,$gene_segment) );
}

sub feed_gene_seq_3end_data { 
	my ($self,$pill_obj,$exp,$sample,$is_antisense,$gene_segment) = @_;
	my $data_obj = $self->get_sequencing_data_obj;
	$self->save_sequencing_data_obj( $data_obj->feed_new_3end_data($pill_obj,$exp,$sample,$is_antisense,$gene_segment) );
}

sub calc_3end_avg_cov { 
	my ($self) = @_;
	my $data_obj = $self->get_sequencing_data_obj;
	$self->save_sequencing_data_obj( $data_obj->calc_all_exp_avg_cov );
}

sub collect_r5UTR_genes { 
	my ($self,$options) = @_;
	my $data_obj = $self->get_sequencing_data_obj;
	$self = $data_obj->collect_r5UTR_candidates($self,$options);
	return $self;
}

sub report_r5utr_candidates { 
	my ($self,$genome_str,$options_ref) = @_;
	my ($locus,$fr,$to,$strand,$desc,$tss) = ($self->get_locus,$self->get_fr,$self->get_to,$self->get_st,$self->get_desc,$self->get_tss_obj);
	my $dom_3end_per_exp_ref = $self->get_r5utr_info;
	return if (!$dom_3end_per_exp_ref);
	my ($tss_pos,$tss_st,$tss_tap,$tss_notap,$tss_ratio); 
	if ($tss) { 
		($tss_pos,$tss_st,$tss_tap,$tss_notap,$tss_ratio,$is_tss_within_coding) = ($tss->get_basic_pill_data,$tss->is_within_coding);
	}
	else { 
		($tss_pos,$tss_st,$tss_tap,$tss_notap,$tss_ratio,$is_tss_within_coding) = ('NA','NA','NA','NA','NA','NA');
	}
	my ($exp_with_max_cov,$pos_3end,$cov_3end,$rep_3end,$element_len) = genes_scripts->get_r5utr_data($dom_3end_per_exp_ref);
	my $seq = 'NA'; 
	if ($tss) { 
		if ($strand eq '+') { 
			$seq = substr($genome_str,$tss_pos-1,abs($tss_pos-$pos_3end)+1);
		}
		else { 
			$seq = substr($genome_str,$tss_pos - abs($tss_pos-$pos_3end),abs($tss_pos-$pos_3end)+1);
			$seq = reverse $seq;
			$seq =~ tr/ATGC/TACG/;
		}
	}
	my $is_known_regulator = ($self->get_known_r5UTR) ? ($self->get_known_r5UTR) : 'NA'; 
	my $return_str = "$locus\t$fr\t$to\t$strand\t$desc\t$is_known_regulator\t";
	$return_str .= "$tss_pos\t$tss_tap\t$tss_ratio\t$is_tss_within_coding\t";
	$return_str .= "$exp_with_max_cov\t$pos_3end\t$cov_3end\t$rep_3end\t$element_len\t";
	$return_str .= "$seq\n";
	return $return_str;
}

sub get_r5utr_data { 
	my ($self,$dom_3end_per_exp_ref) = @_; 
	my @experiments = keys %$dom_3end_per_exp_ref;
	my @string_output_array;
	my ($exp_with_max_cov,$pos,$rep,$element_len);
	my $cov = 0;
	foreach my $exp (@experiments)  { 
		if ($dom_3end_per_exp_ref->{$exp}->{coverage} > $cov) { 
			$exp_with_max_cov = $exp;
			$pos = $dom_3end_per_exp_ref->{$exp}->{pos}; 
			$cov = $dom_3end_per_exp_ref->{$exp}->{coverage}; 
			$rep = $dom_3end_per_exp_ref->{$exp}->{rep};
			$element_len = $dom_3end_per_exp_ref->{$exp}->{length};
		}
	}
	return ($exp_with_max_cov,$pos,$cov,$rep,$element_len);
}

sub is_known_regulator_rfam { 
	my ($self,$gene_array_ref,$rfam_array_ref) = @_;
	my @genes_array = @$gene_array_ref; 
	my @rfam_array = @$rfam_array_ref;
	foreach my $rfam_element_array_ref (@rfam_array) { 
		my ($rfam_fr,$rfam_to,$rfam_st,$rfam_name) = @$rfam_element_array_ref;
		if ($rfam_st eq '+') { 
	GENES: for (my $i = 0; $i <= $#genes_array; $i++) { 
				my ($gene_loc,$gene_fr,$gene_to,$gene_st) = ($genes_array[$i])->get_basic_gene_data;
				my ($ups_gene_loc,$ups_gene_fr,$ups_gene_to,$ups_gene_st) = ($genes_array[$i-1])->get_basic_gene_data;
				if ($gene_fr > $rfam_fr and $ups_gene_to <= $rfam_to and $gene_st eq $rfam_st) { 
					($genes_array[$i])->save_known_r5UTR($rfam_name) unless (($genes_array[$i])->get_known_r5UTR);
					last GENES;
				}
			}
		}
		else { 
	GENES: for (my $i = 0; $i <= $#genes_array; $i++) { 
				my ($gene_loc,$gene_fr,$gene_to,$gene_st) = ($genes_array[$i])->get_basic_gene_data;
				last GENES if (!$genes_array[$i+1]);
				my ($dns_gene_loc,$dns_gene_fr,$dns_gene_to,$dns_gene_st) = ($genes_array[$i+1])->get_basic_gene_data;
				if ($gene_to < $rfam_fr and $dns_gene_fr >= $rfam_fr and $gene_st eq $rfam_st) { 
					($genes_array[$i])->save_known_r5UTR($rfam_name) unless (($genes_array[$i])->get_known_r5UTR);
					last GENES;
				}			
			}
		}			
	}
	return \@genes_array;
}

sub collect_rfam_elements_from_gff { 
	my ($self,$rfam_gff_path) = @_; 
	open(RFAM,$rfam_gff_path) or die "can't open Rfam gff file $rfam_gff_path\n";
	my @rfam_array;
	while (<RFAM>) {
		next if /#/;
		chomp; 
		my @line_sp = split(/\t/); 
		$line_sp[8] =~ m/Alias=(.+);/;
		my @element_array = @line_sp[3,4,6];
		push @element_array,$1;
		push @rfam_array,\@element_array;
	}
	close RFAM;
	return \@rfam_array;
}	

sub collect_ncRNAs_from_genesfile_BSUonly { 
	my ($self,$genes_file) = @_; 
	open(GENES,$genes_file) or die "can't open Rfam gff file $rfam_gff_path\n";
	my @ncRNAs_array;
	while (<GENES>) {
		next unless /BSU_misc_RNA/;
		chomp; 
		my ($contig,$fr,$to,$st,$desc,$ncRNA_name) = split(/\t/); 
		my @array = ($fr,$to,$st,$desc);
		push @ncRNAs_array,\@array;
	}
	close GENES;
	return \@ncRNAs_array;
}	
	
1