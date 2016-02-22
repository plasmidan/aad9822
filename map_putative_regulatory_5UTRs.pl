#!/bin/env perl -w 

###########################################       map_putative_regulatory_5UTRs.pl    ###########################################################
######	combine term-seq and TSS data to map putative regulatory 5'UTRs
######	run example: script running using run_map_putative_r5UTR.pl and an configuration file
#################################################################################################################################################

use strict; 

### Modules 
use common_scripts;
use genes_scripts; 
use TSS_scripts;
use end3_scripts;
use sequencing_data_class;
use Getopt::Long;

#hard coded parameters
my %options =  (output_dir => './',run_name => 'run',make_3p_table_output => 0,'3p_min_rep' => 1,'3p_min_avg_cov' => 2,
				TSS_min_cov => 2, '5UTR_max_len' => 300, min_TAPnoTAP_ratio => 0,use_intergenic => 1,max_3utr_length=>150,
				is_tss_only_5utrs => 0, min_r5utr_len => 70, multiple_mapping => 1, is_dump_genes=>0,report_WT_readthrough => 0,
				report_r5utr_candidates => 0
		);
#values set from configuration		
GetOptions(\%options,"run_name=s","output_dir|o=s","exp_file_path|exp=s","3p_min_rep:i","3p_min_avg_cov:i",
			"make_3p_table_output|3p:i","contig=s","conf_file=s","conf_section=s","make_comparative_3p_table:i",
			"normalize_comparative_table=s","genes_file_path=s","use_5end_data:i","TSS_min_cov:i","max_5utr_length:i",
			"min_TAPnoTAP_ratio:i","use_intergenic:i","max_3utr_length:i","is_tss_only_5utrs:i","min_r5utr_len:i",
			"multiple_mapping:i","is_dump_genes:i","max_5utr_within_coding:i","within_coding_min_TAP:i","report_r5utr_candidates:1",
			"rfam_gff_file_path=s","collect_BSU_specific_ncRNA:i","report_WT_readthrough:i"
			);
		
#make run directory
$options{output_dir} = common_scripts->make_run_directory($options{run_name},$options{output_dir});
print "\ncopying configuration, experiment and Log files to new directory at $options{output_dir}\n\n";
common_scripts->copy_file_to_output_dir( $options{conf_file},$options{output_dir} ); 
common_scripts->copy_file_to_output_dir( $options{'exp_file_path'},$options{output_dir} ); 
common_scripts->move_file_to_output_dir( 'run_log.txt',$options{output_dir} ); 
common_scripts->copy_file_to_output_dir( $options{genes_file_path},$options{output_dir} ); 
open(LOG,'>>',"$options{output_dir}run_log.txt") or die "can't write to log file\n"; 
`echo $options{conf_section} > $options{output_dir}conf_section.txt`; # save which section was used

# Collect genome sequence and genes data
my $genome_str = common_scripts->get_genome_sequence($options{contig});
my $genome_length = length($genome_str);
my $genes_array_ref = genes_scripts->get_genes($options{genes_file_path},$genome_length);
foreach my $gene_obj (@{$genes_array_ref}) { 
	$gene_obj->save_sequencing_data_obj( sequencing_data_class->new_seq_data_obj() );
}

# Annotate genes that are adjacent to an Rfam regulatory element
if ($options{rfam_gff_file_path}) {
	my $rfam_array_ref = genes_scripts->collect_rfam_elements_from_gff($options{rfam_gff_file_path});
	$genes_array_ref = genes_scripts->is_known_regulator_rfam($genes_array_ref,$rfam_array_ref);
}

# Annotate genes that are adjacent to ncRNA (Bacillus specific as the ncRNA are usually not annoatated in gbk in other genomes)
if ($options{collect_BSU_specific_ncRNA}) {
	my $ncRNA_array_ref = genes_scripts->collect_ncRNAs_from_genesfile_BSUonly($options{genes_file_path});
	$genes_array_ref = genes_scripts->is_known_regulator_rfam($genes_array_ref,$ncRNA_array_ref);
}

#collect 5' end sequencing data
my $end5_object_hash_ref = TSS_scripts->make_5end_object_hash( $options{exp_file_path} );
TSS_scripts->map_pillars_to_genes( $end5_object_hash_ref,$genes_array_ref,$options{TSS_min_cov},$options{min_TAPnoTAP_ratio},$options{max_5utr_length} ); 
foreach my $gene_obj (@{$genes_array_ref}) { 
	($gene_obj->get_sequencing_data_obj)->determine_TSS_and_5utr($gene_obj,\%options);
}
for (my $i = 0;$i < @$genes_array_ref; $i++) { 
	if ($i > 0) { 
		$genes_array_ref->[$i]->save_upstream_gene($genes_array_ref->[$i-1]);
	}
	if ($i <@$genes_array_ref) { 
		$genes_array_ref->[$i]->save_downstream_gene($genes_array_ref->[$i+1]);
	}
	$genes_array_ref->[$i] = sequencing_data_class->is_legit_withicoding_tss($genes_array_ref->[$i]);
}

#collect term-seq sequencing data
my $end3_data_hash_ref = end3_scripts->collect_3end_data_into_hash( $options{exp_file_path},$options{'3p_min_avg_cov'},$options{'3p_min_avg_rep'} );
my $multiple_map_hash_ref = end3_scripts->map_pillars_to_genes($end3_data_hash_ref,$genes_array_ref,\%options);

for (my $i=0; $i< @$genes_array_ref;$i++) { 
	$genes_array_ref->[$i]->calc_3end_avg_cov;
	$genes_array_ref->[$i] = $genes_array_ref->[$i]->collect_r5UTR_genes(\%options);
}

# output list of putative regulatory 5'UTRs
if ($options{report_r5utr_candidates} == 1) { 
	open(CAND,'>',"$options{output_dir}r5UTR_candidate_table.txt");
	print CAND "Locus_tag\tfrom\tto\tst\tdesc\tis_known_ncRNA?\tTSS_pos\tTAP\tRatio\tis_tss_within_coding\tDiscovered_in\tr5utr_3end_pos\t3'_cov\trepeated_in\tr5UTR_len\tSeq\n";
	for (my $i=0; $i< @$genes_array_ref;$i++) { 
		my $output_str = ($genes_array_ref->[$i])->report_r5utr_candidates($genome_str,\%options);
		print CAND "$output_str" if ($output_str);
	}
}
