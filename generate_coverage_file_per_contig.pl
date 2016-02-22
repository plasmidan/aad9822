#!/bin/env perl -w 
use strict; 

###############################################   generate_coverage_file_per_contig_metatranscriptome.pl   ###########################################
######	from RNA-seq mapped reads to base resolution coverage map per contig
######	run example: perl generate_coverage_file_per_contig_metatranscriptome.pl [align_file_dir] [genome_fasta] [output_dir] [is_uniq_only?]
######################################################################################################################################################

#command line inputs
my ($align_file_dir,$genome_fasta,$output_dir,$is_unique_only) = @ARGV;
warn "unique reads only\n" if ($is_unique_only);

#measure contig lengths
open(FA,$genome_fasta) or die "can't open genome fasta $genome_fasta\n";
my %contig_lengths;
while(<FA>) { 
	/^>(\S+)/; 
	my $contig = $1; 
	my $seq = <FA>;
	my $contig_len  = length($seq);
	$contig_lengths{$contig} = $contig_len;
}
close FA;

#collect align file names from dir
my @align_files = `ls $align_file_dir`;
my $c = 0;
foreach my $align_file (@align_files) { 
	$c++; 
	warn "$c contigs done\n" if ($c % 1000 == 0);
	chomp $align_file; 
	$align_file =~ /^(\S+)\.(\S+)\.align/;
	my ($contig_name,$condition) = ($1,$2);
	next unless ($contig_lengths{$contig_name});
	my $genome_len = $contig_lengths{$contig_name};
	open(ALIGN,"$align_file_dir/$align_file") or die "can't open align $align_file_dir/$align_file\n";
	my @forward_array = @{[(0) x ($genome_len)]};
	my @reverse_array = @{[(0) x ($genome_len)]};
	while(<ALIGN>) { 
		chomp;
		my ($type,$fr,$to,$st) = split(/\t/);
		next if ($is_unique_only and $type eq "R");
		($fr,$to) = ($to,$fr) if ($fr > $to); 
		$fr = $genome_len if ($fr > $genome_len); 
		$to = $genome_len if ($to > $genome_len); 
		if ($st eq "+") { 
			for (my $i=$fr-1;$i<=$to-1;$i++) {$forward_array[$i]++;};
		}
		else { 
			for (my $i=$fr-1;$i<=$to-1;$i++) {$reverse_array[$i]++;};
		}
	}
	#output
	my $output_name = ($is_unique_only) ? "$output_dir/$contig_name.$condition.unique.coverage" : "$output_dir/$contig_name.$condition.coverage";
	open(OUT,'>',$output_name) or die "can't write out $contig_name\n";
	for (my $i = 0; $i<= $#forward_array; $i++) { 
		my ($for,$rev,$merge) = ( $forward_array[$i],$reverse_array[$i],($forward_array[$i]+$reverse_array[$i]) ); 
		print OUT "$contig_name\t$merge\t$for\t$rev\n";
	}
	close OUT;
}