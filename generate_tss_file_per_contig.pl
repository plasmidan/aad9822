#!/bin/env perl -w 
use strict; 

###############################################   generate_tss_file_per_contig_metatranscriptome.pl   ###########################################
######	from 5' end mapped reads to position specific 5' sites per contig
######	run example: perl generate_tss_file_per_contig_metatranscriptome.pl [contig_list] [align_file_dir] [output_dir] [is_uniq_only?]
######################################################################################################################################################

#command line inputs
my ($contig_list,$align_file_dir,$output_dir,$is_unique) = @ARGV;  
warn "using unique only reads\n" if ($is_unique);

my @align_files = `ls $align_file_dir`; 
@align_files = grep {$_ =~ /(tap|notap)/} @align_files;
open(LIST,$contig_list) or die "contig list\n";
while(my $contig = <LIST>) { 
	warn "$contig";
	chomp $contig;
	my @contig_files = grep {/$contig.+\.align/} @align_files; 
	foreach my $align_file (@contig_files) { 
		chomp $align_file;
		$align_file =~ /^(\S+).align/;
		my $file_out_name = ($is_unique) ? "$1.unique.s2" : "$1.s2"; 
		open(ALN,"$align_file_dir/$align_file") or die "can't open align $align_file for $contig\n";
		my %pills;
		while(<ALN>) {
			chomp;
			my ($type,$pos,$to,$st) = split(/\t/);
			next if ($type eq "R" and $is_unique);
			$pills{$pos}->{$st}++;
		}
		close ALN; 
		open(OUT,'>',"$output_dir/$file_out_name") or die "can't write to tss s2 $file_out_name for $contig\n";
		next unless (%pills);
		foreach my $site (sort {$a<=>$b} keys %pills) { 
			foreach my $st (keys %{$pills{$site}}) { 
				my $cov = $pills{$site}->{$st}; 
				print OUT "$contig\t$site\t$st\t$cov\t-1\n";
			}
		}
		close OUT;
	}
}