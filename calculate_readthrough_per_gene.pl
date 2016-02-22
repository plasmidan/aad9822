#!/bin/env perl -w 
use strict;

###############################################   calculate_readthrough_per_gene_metatranscriptome.pl    ###########################################################
######	Use RNA-seq coverage to calcaulte the potential readthrough for each gene in the metatranscriptome experiment
######	run example: perl calculate_readthrough_per_gene_metatranscriptome.pl [merged_reads_per_gene_file] [directory_of_contig_specific_coverage_files] [is_uniq_only?]
####################################################################################################################################################################

#command line inputs
my ($reads_per_gene_file,$coverage_dir,$is_unique) = @ARGV; #genes_file for strand info 
warn "using unique reads only\n" if ($is_unique);

#hard coded parameters
my ($min_5utr_len,$max_5utr) = (150,200);
my $control_suffix = "control_ssWT";
my $antibiotic_suffix = "lincomycin_ssWT";

my (%hash,%genes_hash);
open(GENES,$reads_per_gene_file) or die "can't open reads_per_gene_file $reads_per_gene_file\n";
while(<GENES>) { 
	my @line = split(/\t/); 
	$genes_hash{$line[3]} .= $_;
}

my $c=0;
foreach my $contig (keys %genes_hash) { 
	my $cont_file = ($is_unique) ? "./$coverage_dir/$contig.$control_suffix.unique.coverage" : "./$coverage_dir/$contig.$control_suffix.coverage";
	my $ant_file = ($is_unique) ? "./$coverage_dir/$contig.$antibiotic_suffix.unique.coverage" : "./$coverage_dir/$contig.$antibiotic_suffix.coverage";

	open(COV_CONT,$cont_file) or die "$cont_file\n";
	open(COV_ANT,$ant_file) or next;
	
	my (@cov_array_control,@cov_array_antibiotic); 
	my $contig_len = 0;
	while (my $cont_line = <COV_CONT>) {
		my $ant_line = <COV_ANT>; 
		my @cont_line_sp = split(/\t/,$cont_line); 
		my @ant_line_sp = split(/\t/,$ant_line); 
		
		die "$contig $contig_len\n" if (!exists $cont_line_sp[1]);
		push @cov_array_control,$cont_line_sp[1];
		push @cov_array_antibiotic,$ant_line_sp[1];
		$contig_len++;
	}

	my @contig_genes = @{ &process_genes_line($genes_hash{$contig},$contig_len,$min_5utr_len,$max_5utr)};
	
	GENES:
	foreach my $gene (@contig_genes) { 
		if ($gene->{int_size} < $min_5utr_len ) {  #look only at genes that have an intergenic region > set minimum value 
			print "$gene->{control_reads}\t$gene->{antibiotics_reads}\t0\t0\t$contig\t";
			print "$gene->{id}\t$gene->{fr}\t$gene->{to}\t$gene->{st}\t$gene->{desc}\t";
			print "0\t0\t0\t0\n"; 		
			next GENES;
		}
		my ($gene_cov_control,$gene_cov_antibiotic,$reg_cov_control,$reg_cov_antibiotic) = (0,0,0,0);
		my $utr_fr = ($gene->{st} eq "+") ? $gene->{utr_fr} : $gene->{to};
		my $utr_to = ($gene->{st} eq "+") ? $gene->{fr} : $gene->{utr_fr};		
		for (my $i=$gene->{fr}-1;$i<=$gene->{to}-1;$i++) { #coverage over gene
			$gene_cov_control+= $cov_array_control[$i];
			$gene_cov_antibiotic+= $cov_array_antibiotic[$i];
		}	
		for (my $i=$utr_fr-1;$i<=$utr_to-2;$i++) { #coverage over regulatory region
			$reg_cov_control+= $cov_array_control[$i];
			$reg_cov_antibiotic+= $cov_array_antibiotic[$i];
		}			
		#calculate the average coverage and readthrough
		my ($gene_control_avg,$gene_antibiotic_avg,$reg_control_avg,$reg_antibiotic_avg) = (0,0,0,0);
		($gene_control_avg,$gene_antibiotic_avg) = ( $gene_cov_control / abs($gene->{fr}-$gene->{to}+1),$gene_cov_antibiotic / abs($gene->{fr}-$gene->{to}+1) );
		($reg_control_avg,$reg_antibiotic_avg) = ( $reg_cov_control / ($utr_to-$utr_fr+1),$reg_cov_antibiotic / ($utr_to-$utr_fr+1) );
		my ($control_readthrough,$antibiotic_readthrough) = (0,0);
		$control_readthrough = int 100* $gene_control_avg / ($reg_control_avg+1);
		$antibiotic_readthrough = int 100* $gene_antibiotic_avg / ($reg_antibiotic_avg+1);		
		$control_readthrough = 100 if ($control_readthrough > 100);
		$antibiotic_readthrough = 100 if ($antibiotic_readthrough > 100);	
		$gene_control_avg =~ /(\d+\.*\d{0,1})/; $gene_control_avg = $1;
		$gene_antibiotic_avg =~ /(\d+\.*\d{0,1})/; $gene_antibiotic_avg = $1;
		$reg_control_avg =~ /(\d+\.*\d{0,1})/; $reg_control_avg = $1;
		$reg_antibiotic_avg =~ /(\d+\.*\d{0,1})/; $reg_antibiotic_avg = $1;
		#output
		print "$gene->{control_reads}\t$gene->{antibiotics_reads}\t$control_readthrough\t$antibiotic_readthrough\t$contig\t";
		print "$gene->{id}\t$gene->{fr}\t$gene->{to}\t$gene->{st}\t$gene->{desc}\t";
		print "$gene_control_avg\t$reg_control_avg\t$gene_antibiotic_avg\t$reg_antibiotic_avg\n"; 	
	}
	close COV_CONT;
	close COV_ANT;
}


###########################################################################################
############################# Additional Functions ########################################
###########################################################################################

sub process_genes_line 
{ 
	my ($gene_lines,$contig_len,$min_5utr_len,$max_5utr) = @_;
	my @lines = split(/\n/,$gene_lines);
	my @genes;
	foreach my $line (@lines) {
		my ($control_reads,$antibiotics_reads,$gene_id,$contig,$fr,$to,$st,$desc) = split(/\t/,$line);
		my $gene->{fr} = $fr;
		$gene->{to} = $to;
		$gene->{st} = $st;
		$gene->{id} = $gene_id;
		$gene->{desc} = $desc;
		$gene->{control_reads} = $control_reads;
		$gene->{antibiotics_reads} = $antibiotics_reads;
		push @genes,$gene;
	}
	my @out_genes;
	for (my $i=0;$i<=$#genes;$i++) {
		next if ( $genes[$i]->{fr} - $min_5utr_len < 0 and $genes[$i]->{st} eq "+" );
		next if ( $genes[$i]->{to} + $min_5utr_len > $contig_len and $genes[$i]->{st} eq '-' );
		if ($genes[$i]->{st} eq "+") { 
			$genes[$i]->{int_size} = ($i == 0) ? $genes[$i]->{fr} : ($genes[$i]->{fr} - $genes[$i-1]->{to} + 1);
			$genes[$i]->{utr_fr} = ( $genes[$i]->{int_size} > $max_5utr ) ? ($genes[$i]->{fr} - $max_5utr + 1) : ($genes[$i]->{fr} - $genes[$i]->{int_size} + 1);
		}
		else { 
			$genes[$i]->{int_size} = ($i == $#genes) ? ($contig_len - $genes[$i]->{to} + 1 ) : ($genes[$i+1]->{fr} - $genes[$i]->{to});
			$genes[$i]->{utr_fr} = ( $genes[$i]->{int_size} > $max_5utr ) ? ($genes[$i]->{to} + $max_5utr + 1) : ($genes[$i]->{to} + $genes[$i]->{int_size});
		}
		push @out_genes,$genes[$i];
	}
	return \@out_genes;
}