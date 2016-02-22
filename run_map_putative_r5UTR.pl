#!/bin/env perl -w 
use strict;

############################################    run script for map_putative_regulatory_5UTRs.pl    ##############################################
##### uses a configuration file to initatie a run of map_putative_regulatory_5UTRs.pl (see example .ini file)
##### run example: perl run_map_putative_r5UTR.pl ex_conf.ini bsubtilis_LB
#################################################################################################################################################

use Config::IniFiles;

my ($path_to_conf,$section_in_conf) = @ARGV;
my $config_obj = new Config::IniFiles(-file => $path_to_conf);
my @section_parameters = $config_obj->Parameters($section_in_conf);

my $commands = "perl map_putative_regulatory_5UTRs.pl --conf_file $path_to_conf --conf_section $section_in_conf";

foreach my $parameter (@section_parameters) { 
	my $parameter_value = $config_obj->val($section_in_conf,$parameter);
	$commands .= " --$parameter  $parameter_value ";
}

warn "Running script as:\n\n$commands\n";
my $save_command_run = "echo $commands > ./run_log.txt";
system($save_command_run); 
system($commands);












# my $conf_obj = new config_manager_dd($cfg_path);
