#!/usr/bin/env perl

# Calls detect_potential_cross_contamination.pl or other helper scripts depending on
# user input.

use strict;
use warnings;

# option entered by user
my $DETECT_POTENTIAL_CROSS_CONTAMIANTION_OPTION = "detect_cross_contam";

# script run in response to option entered by user
my $DETECT_POTENTIAL_CROSS_CONTAMINATION_SCRIPT_FILE_PATH = "detect_potential_cross_contamination.pl";


# if no command line arguments supplied, prints options
if(!scalar @ARGV) # no command line arguments supplied
{
	print STDOUT "\n";
	print STDOUT "Usage: polyphonia.pl <command> [options]\n\n";
	print STDOUT "\tdetect_cross_contam\tDetect potential cross-contamination.\n";
	print STDOUT "\n\n";
	exit;
}

# collects all command line arguments after the first one
my @arguments_to_pass_along = @ARGV[1..$#ARGV];
my $arguments_to_pass_along_string = join(" ", @arguments_to_pass_along);

# passes command line arguments to appropriate helper script
my $command = $ARGV[0];
if($command eq $DETECT_POTENTIAL_CROSS_CONTAMIANTION_OPTION)
{
	# passes arguments to detect potential cross-contamination helper script
	exec("$DETECT_POTENTIAL_CROSS_CONTAMINATION_SCRIPT_FILE_PATH $arguments_to_pass_along_string");
}
else
{
	# command not recognized
	print STDERR "Error: command ".$command." not recognized. Exiting.\n";
	die;
}

# June 14, 2021
