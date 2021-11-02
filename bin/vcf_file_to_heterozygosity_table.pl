#!/usr/bin/env perl

# Reads in vcf file output produced by LoFreq call and prints human-readable heterozygosity
# table. Optionally filters output by read depth and minor allele readcount and frequency.

# See https://csb5.github.io/lofreq/commands/#call for information on LoFreq call.

# Output columns, tab-separated:
# - Reference
# - Position (relative to reference, 1-indexed)
# - Major allele
# - Major allele readcount
# - Major allele frequency
# - Minor allele
# - Minor allele readcount
# - Minor allele frequency

# Usage:
# perl vcf_file_to_heterozygosity_table.pl [vcf file output by LoFreq]
# [1 to filter output, 0 if not]

# Prints to console. To print to file, use
# perl vcf_file_to_heterozygosity_table.pl [vcf file output by LoFreq]
# [1 to filter output, 0 if not] > [output file path]


use strict;
use warnings;
use Math::Round;

my $vcf_file = $ARGV[0]; # output of LoFreq, uncompressed or .gz
my $filter_output = $ARGV[1]; # if 1, filters output using thresholds for including a heterozygous position; if 0, includes all positions with heterozygosity in output


# thresholds for including a heterozygous position
my $MINIMUM_MINOR_ALLELE_READCOUNT = 10;
my $MINIMUM_MINOR_ALLELE_FREQUENCY = 0.03; # 3%
my $MINIMUM_READ_DEPTH = 100;

# columns in input vcf file
my $REFERENCE_COLUMN = 0; # 0-indexed
my $POSITION_COLUMN = 1;
my $REFERENCE_ALLELE_COLUMN = 3;
my $ALTERNATIVE_ALLELES_COLUMN = 4;
my $INFO_COLUMN = 7;

# in info column of input vcf file
my $INFO_DELIMITER = ";";

# in both input and output
my $DELIMITER = "\t";
my $NEWLINE = "\n";
my $NO_DATA = "NA";
my $NO_DATA_INPUT_VALUE = ".";

# for printing allele frequencies
my $DECIMALS_TO_ROUND_TO = 6; # same as in vcf files
my $NUMBER_ALLELES_TO_PRINT_PER_POSITION = 2;


# verifies that input file exists
if(!$vcf_file or !-e $vcf_file)
{
	print STDERR "Error: vcf file $vcf_file does not exist. Exiting.\n";
	die;
}

# unzips input file if it is compressed
my $created_unzipped_vcf_file = 0;
if($vcf_file =~ /^(.*)[.]gz$/) # assumes file is zipped if ends in .gz, not otherwise
{
	$vcf_file = $1;
	
	# does not unzip if there is already an unzipped file
	if(-e $vcf_file)
	{
		print STDERR "Warning: using existing unzipped file $vcf_file\n";
	}
	else
	{
		`gzip -dk $vcf_file`; # unzips input file (does not delete zipped file)
		$created_unzipped_vcf_file = 1; # marks unzipped file for deletion once we're done
	}
}

# reads in vcf file
my %chr_to_position_to_allele_to_readcount = (); # key: assembly reference -> key: position -> key: allele -> value: readcount
my %chr_to_position_to_read_depth = (); # key: assembly reference -> key: position -> value: read depth
open VCF_FILE, "<$vcf_file" || die "Could not open $vcf_file to read; terminating =(\n";
while(<VCF_FILE>) # for each line in the file
{
	chomp;
	my $line = $_;
	if($line =~ /\S/ and $line !~ /^#/) # not empty line, not header line
	{
		my @items_in_row = split($DELIMITER, $line);
		my $info = $items_in_row[$INFO_COLUMN];
		
		# retrieves position information
		my $assembly_reference = $items_in_row[$REFERENCE_COLUMN];
		my $position = $items_in_row[$POSITION_COLUMN];
		
		# retrieves allele identities
		my $allele_ref = $items_in_row[$REFERENCE_ALLELE_COLUMN]; # higher-frequency allele
		my $allele_alt = $items_in_row[$ALTERNATIVE_ALLELES_COLUMN]; # lower-frequency allele
		
		# retrieves read depth
		my $read_depth = 0;
		if($info =~ /DP=(\d+)$INFO_DELIMITER/)
		{
			$read_depth = $1;
		}
		
		# retrieves allele readcounts
		my $allele_ref_readcount = 0;
		my $allele_alt_readcount = 0;
		if($info =~ /DP4=(\d+),(\d+),(\d+),(\d+)$/)
		{
			$allele_ref_readcount = $1 + $2;
			$allele_alt_readcount = $3 + $4;
		}
		
		# prints warning if allele is not one character
		if($allele_ref !~ /^\w$/ or $allele_alt !~ /^\w$/)
		{
			print STDERR "Warning: not a base substitution (alleles ".$allele_ref
				." and ".$allele_alt.") at position ".$position." in vcf file:\n\t".$vcf_file."\n";
		}
		else # only saves base substitutions for printing
		{
			# saves allele readcounts if they are non-zero
			if($allele_ref_readcount > 0)
			{
				$chr_to_position_to_allele_to_readcount{$assembly_reference}{$position}{$allele_ref} = $allele_ref_readcount;
			}
			if($allele_alt_readcount > 0)
			{
				$chr_to_position_to_allele_to_readcount{$assembly_reference}{$position}{$allele_alt} = $allele_alt_readcount;
			}
		}
		
		# saves read depth
		$chr_to_position_to_read_depth{$assembly_reference}{$position} = $read_depth;
	}
}
close VCF_FILE;


# prints human-readable row describing each position
for my $assembly_reference(sort keys %chr_to_position_to_allele_to_readcount)
{
	for my $position(sort {$a <=> $b} keys %{$chr_to_position_to_allele_to_readcount{$assembly_reference}})
	{
		# checks that this position has at least 2 alleles with non-zero readcount
		if(scalar keys %{$chr_to_position_to_allele_to_readcount{$assembly_reference}{$position}} >= 2)
		{
			# saves start of line
			my $output_line = $assembly_reference.$DELIMITER.$position;
			
			# retrieves read depth at this position
			my $read_depth = $chr_to_position_to_read_depth{$assembly_reference}{$position};
		
			# saves info on the two alleles with non-zero readcount, from greatest readcount to lowest
			my $number_alleles_printed = 0;
			my $minor_allele_readcount = -1;
			my $minor_allele_frequency = -1;
			for my $allele(sort {$chr_to_position_to_allele_to_readcount{$assembly_reference}{$position}{$b}
					<=> $chr_to_position_to_allele_to_readcount{$assembly_reference}{$position}{$a}}
				keys %{$chr_to_position_to_allele_to_readcount{$assembly_reference}{$position}})
			{
				if($number_alleles_printed < $NUMBER_ALLELES_TO_PRINT_PER_POSITION)
				{
					# retrieves allele readcount and frequency
					my $allele_readcount = $chr_to_position_to_allele_to_readcount{$assembly_reference}{$position}{$allele};
					my $allele_frequency = round_value($allele_readcount / $read_depth, $DECIMALS_TO_ROUND_TO);
			
					# adds allele info to output line
					$output_line .= $DELIMITER.$allele;
					$output_line .= $DELIMITER.$allele_readcount;
					$output_line .= $DELIMITER.$allele_frequency;
					
					$number_alleles_printed++;
					
					# saves info on minor allele for filtering
					if($number_alleles_printed == 2)
					{
						$minor_allele_readcount = $allele_readcount;
						$minor_allele_frequency = $allele_frequency;
					}
				}
			}
			
			# prints output line for this position if we aren't filtering or if position passes thresholds
			if(!$filter_output or ($read_depth >= $MINIMUM_READ_DEPTH
				and $minor_allele_readcount >= $MINIMUM_MINOR_ALLELE_READCOUNT
				and $minor_allele_frequency >= $MINIMUM_MINOR_ALLELE_FREQUENCY))
			{
				print $output_line.$NEWLINE;
			}
		}
	}
}


# deletes unzipped vcf file if we created it
if($created_unzipped_vcf_file)
{
	`rm $vcf_file`;
}


# input: floating point value (example: 3.14159)
# output: value rounded to nth decimal point (example: 3.1)
sub round_value
{
	my $value = $_[0];
	my $number_decimals = $_[1];
	return sprintf("%.".$number_decimals."f", $value);
}


# May 11, 2020
# June 1, 2021
# October 29, 2021
# November 2, 2021
