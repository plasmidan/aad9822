package common_scripts; 

sub make_run_directory { 
	my ($self,$run_name,$output_dir) = @_;
	chomp (my $time_stamp = `date +"%m-%d-%y-%T"`); 
	$time_stamp =~ tr/:/_/;
	$run_name .= "_$time_stamp/";
	my $run_dir = $output_dir . "$run_name";
	my $make_dir_cmd = "mkdir -p $run_dir";
	system($make_dir_cmd);
	return $run_dir;
}

sub parse_exp_file {  	
	my ($self,$exp_file_path,$exp_type) = @_; 
	my @exp_project_dirs;
	my @exp_run_dirs;
	my @sample_names;
	my @exp_names;
	
	open(EXP,$exp_file_path) or die "can't open experiment file $exp_file_path\n";
	while (my $sample_line = <EXP>) { 
		chomp $sample_line; 
		my ($lib_type,$project_dir,$run_dir,$sample_name,$exp_name) = split(/\t/,$sample_line);
		next unless ($lib_type eq $exp_type);
		push @exp_project_dirs,$project_dir;	
		push @exp_run_dirs,$run_dir;
		push @sample_names,$sample_name;
		push @exp_names,$exp_name;
	}

	my @s2_file_paths;
	foreach ($i=0; $i<=$#exp_names; $i++) { 
		my $exp_type = $exp_names[$i];
		my $s2_path = "/home/labs/sorek/repos/projects/$exp_project_dirs[$i]/$exp_run_dirs[$i]/intermediate/$exp_run_dirs[$i].s2"; 
		`sort -k3,3 -n $s2_path > $s2_path.sorted`;
		push @s2_file_paths,$s2_path;
	}
	return (\@s2_file_paths,\@sample_names,\@exp_names);
}

sub move_file_to_output_dir { 
	my ($self,$file_path,$output_dir) = @_; 
	my $cmd = "mv $file_path $output_dir";
	system($cmd);
}
sub copy_file_to_output_dir { 
	my ($self,$file_path,$output_dir) = @_; 
	my $cmd = "cp $file_path $output_dir";
	system($cmd);
}

sub get_genome_sequence { 
	my ($self,$contig) = @_; 
	open(FA,"/home/labs/sorek/repos/genomes/$contig.fasta") or die "Fasta file not provided in get_genome_seq subroutine\n";
	my $genome_str; 
	<FA>;
	while (<FA>) {
		chomp;
		$genome_str .= $_; 
	}
	return $genome_str;
}	

sub median { 
	my ($self,$array_ref) = @_; 
	my $count = scalar @$array_ref; 
	my @array = sort { $a <=> $b } @$array_ref; 
	if ($count % 2) { 
	return $array[int($count/2)]; 
	} else { 
	return ($array[$count/2] + $array[$count/2 - 1]) / 2; 
	} 
}

sub average {
	my ($self,$values_ref) = @_;
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

sub std_dev {
	my ($self,$average, $values_ref) = @_;
	my @values = @{$values_ref};
	my $count = scalar @values;
	my $std_dev_sum = 0;
	foreach (@values) { 
		next if ($_ eq 'NA');
		$std_dev_sum += ($_ - $average) ** 2; 
	}
	return $count ? sqrt($std_dev_sum / $count) : 0;
}

return 1