#!/usr/bin/env perl

# reads in vcf file and prints human-readable heterozygosity table

# includes positions with with all of:
#   base substitution (indels not included)
#   minor allele readcount >= 10 reads
#   minor allele frequency >= 3%

# output columns, tab-separated:
#   reference
#   position (relative to reference, 1-indexed)
#   major allele
#   major allele readcount
#   major allele frequency
#   minor allele
#   minor allele readcount
#   minor allele frequency

use strict;
use warnings;

my $vcf_file = $ARGV[0]; # output of GATK or LoFreq, uncompressed or .gz


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

# thresholds for including a heterozygous position
my $MINIMUM_MINOR_ALLELE_READCOUNT = 10;
my $MINIMUM_MINOR_ALLELE_FREQUENCY = 0.03; # 3%


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
open VCF_FILE, "<$vcf_file" || die "Could not open $vcf_file to read; terminating =(\n";
while(<VCF_FILE>) # for each line in the file
{
	chomp;
	my $line = $_;
	if($line =~ /\S/ and $line !~ /^#/) # not empty line, not header line
	{
		my @items_in_row = split($DELIMITER, $line);
		my $info = $items_in_row[$INFO_COLUMN];
		
		# checks if line refers to position of heterozygosity and retrieves allele frequency
		my $line_is_position_of_heterozygosity = 0;
		my $allele_1_frequency;
		if($info =~ /ABHet=([\d.]+)$INFO_DELIMITER/) # heterozygous call in GATK file
		{
			$allele_1_frequency = $1;
			$line_is_position_of_heterozygosity = 1;
		}
		elsif($info =~ /ABHom=([\d.]+)$INFO_DELIMITER/) # homozygous call in GATK file
		{
			# ignore
		}
		elsif($info =~ /AF=([\d.]+)$INFO_DELIMITER/) # allele frequency in LoFreq file
		{
			$allele_1_frequency = 1 - $1;
			$line_is_position_of_heterozygosity = 1;
		}
		
		# retrieves all other valuable information from line
		if($line_is_position_of_heterozygosity)
		{
			# retrieves allele identities
			my $allele_1 = $items_in_row[$REFERENCE_ALLELE_COLUMN]; # higher-frequency allele
			my $allele_2 = $items_in_row[$ALTERNATIVE_ALLELES_COLUMN]; # lower-frequency allele
			
			# swaps higher- and lower-frequency alleles if necessary
			if($allele_1_frequency < 0.5)
			{
				$allele_1_frequency = 1 - $allele_1_frequency;
				
				# swaps alleles
				$allele_1 = $items_in_row[$ALTERNATIVE_ALLELES_COLUMN];
				$allele_2 = $items_in_row[$REFERENCE_ALLELE_COLUMN];
			}
			my $allele_2_frequency = 1 - $allele_1_frequency;
			
			if($items_in_row[$ALTERNATIVE_ALLELES_COLUMN] !~ /^\w$/
				and $items_in_row[$ALTERNATIVE_ALLELES_COLUMN] ne $NO_DATA_INPUT_VALUE)
			{
				print STDERR "Error: more than two alleles. This script isn't written to handle that. Exiting.\n";
				die;
			}
			
			# retrieves read depth
			my $read_depth = 0;
			if($info =~ /DP=(\d+)$INFO_DELIMITER/)
			{
				$read_depth = $1;
			}
			
			# calculated read count of each allele
			# can also get from AD values
			my $allele_1_readcount = round($allele_1_frequency * $read_depth);
			my $allele_2_readcount = round($allele_2_frequency * $read_depth);
			
			# retrieves other information
			my $assembly_reference = $items_in_row[$REFERENCE_COLUMN];
			my $position = $items_in_row[$POSITION_COLUMN];
	
			if( # major and minor allele are each one base
				$allele_1 =~ /^\w$/ and $allele_2 =~ /^\w$/)
			
				# this position's minor allele frequency and readcount pass our thresholds
# 				and $allele_2_readcount >= $MINIMUM_MINOR_ALLELE_READCOUNT
# 				and $allele_2_frequency >= $MINIMUM_MINOR_ALLELE_FREQUENCY)
			{
				# prints human-readable row describing this position
				print $assembly_reference.$DELIMITER;
				print $position.$DELIMITER;
			
				print $allele_1.$DELIMITER;
				print $allele_1_readcount.$DELIMITER;
				print $allele_1_frequency.$DELIMITER;
			
				print $allele_2.$DELIMITER;
				print $allele_2_readcount.$DELIMITER;
				print $allele_2_frequency.$NEWLINE;
			}
		}
	}
}
close VCF_FILE;

# deletes unzipped vcf file if we created it
if($created_unzipped_vcf_file)
{
	`rm $vcf_file`;
}

# input: floating point value (example: 3.14159)
# output: value rounded to first decimal point (example: 3.1)
use Math::Round;
sub round_value
{
	my $value = $_[0];
	return sprintf("%.1f", $value);
}

# May 11, 2020
# June 1, 2021
