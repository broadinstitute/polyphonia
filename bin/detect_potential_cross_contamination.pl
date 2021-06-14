#!/usr/bin/env perl

# Detects potential cross-contamination by comparing consensus-level and minor alleles.

use strict;
use warnings;

# tools used:
use Math::Round;
use Parallel::ForkManager; # download here: https://metacpan.org/pod/Parallel::ForkManager
my $LOFREQ_EXECUTABLE_FILE_PATH = "lofreq";
my $MAFFT_EXECUTABLE_FILE_PATH = "mafft";
my $VCF_TO_HETEROZYGOSITY_TABLE_SCRIPT_FILE_PATH = "vcf_file_to_heterozygosity_table.pl";
my $PLATE_VISUALIZATION_FILE_PATH = "/opt/polyphonia/visualize_potential_cross_contamination.R";

# plate map input file:
my $PLATE_MAP_SAMPLE_COLUMN = 0;
my $PLATE_MAP_POSITION_COLUMN = 1;

# for reading plate map:
my %LETTER_TO_NEXT_LETTER = (
	"A" => "B", "B" => "C", "C" => "D", "D" => "E", "E" => "F", "F" => "G", "G" => "H",
	"H" => "I", "I" => "J", "J" => "K", "K" => "L", "L" => "M", "M" => "N", "N" => "O",
	"O" => "P", "P" => "Q", "Q" => "R", "R" => "S", "S" => "T", "T" => "U", "U" => "V",
	"V" => "W", "W" => "X", "X" => "Y", "Y" => "Z", "Z" => "");
my %LETTER_TO_PREVIOUS_LETTER = (
	"A" => "", "B" => "A", "C" => "B", "D" => "C", "E" => "D", "F" => "E", "G" => "F",
	"H" => "G", "I" => "H", "J" => "I", "K" => "J", "L" => "K", "M" => "L", "N" => "M",
	"O" => "N", "P" => "O", "Q" => "P", "R" => "W", "S" => "R", "T" => "S", "U" => "T",
	"V" => "U", "W" => "V", "X" => "W", "Y" => "X", "Z" => "Y");

# intermediate file heterozygosity table columns:
my $HETEROZYGOSITY_TABLE_POSITION_COLUMN = 1; # (0-indexed)
my $HETEROZYGOSITY_TABLE_MAJOR_ALLELE_COLUMN = 2;
my $HETEROZYGOSITY_TABLE_MAJOR_ALLELE_READCOUNT_COLUMN = 3;
my $HETEROZYGOSITY_TABLE_MAJOR_ALLELE_FREQUENCY_COLUMN = 4;
my $HETEROZYGOSITY_TABLE_MINOR_ALLELE_COLUMN = 5;
my $HETEROZYGOSITY_TABLE_MINOR_ALLELE_READCOUNT_COLUMN = 6;
my $HETEROZYGOSITY_TABLE_MINOR_ALLELE_FREQUENCY_COLUMN = 7;

# options for output
my $ROUND_PERCENTAGES = 1;
my $ALLELE_LIST_SEPARATOR = "; ";
my $PERCENTAGES_LIST_SEPARATOR = ", ";
my $PLATE_MAP_POSITION_LIST_SEPARATOR = ", ";
my $FILE_LIST_SEPARATOR = ", ";

my $DELIMITER = "\t";
my $NEWLINE = "\n";
my $NO_DATA = "NA";


# default values for command line arguments
my $DEFAULT_MAXIMUM_ALLOWED_MISMATCHES = 1;
my $DEFAULT_MINIMUM_GENOME_COVERAGE = 0.98;

my $DEFAULT_MINIMUM_MINOR_ALLELE_READCOUNT = 10; # ignore minor alleles with readcount <10 (does not consider the position to have heterozygosity)
my $DEFAULT_MINIMUM_MINOR_ALLELE_FREQUENCY = 0.03; # 3%

my $DEFAULT_OUTPUT_FILE_NAME = "potential_cross-contamination.txt";
my $DEFAULT_OVERWRITE = 0; # false
my $DEFAULT_CORES_TO_USE = 1;
my $DEFAULT_VERBOSE = 1; # true

my $DEFAULT_COMPARE_DIRECT_NEIGHBORS = 1; # true
my $DEFAULT_COMPARE_DIAGONAL_NEIGHBORS = 0; # false
my $DEFAULT_COMPARE_ROW = 0; # false
my $DEFAULT_COMPARE_COLUMN = 0; # false
my $DEFAULT_COMPARE_WHOLE_PLATE_MAP = 0; # false


# generates default directory and output file from current working directory
my $default_temp_intermediate_files_directory = retrieve_current_working_directory(); # "CURRENT WORKING DIRECTORY"; # current working directory
my $default_output_file = $default_temp_intermediate_files_directory.$DEFAULT_OUTPUT_FILE_NAME;


# if no command line arguments supplied, prints options
if(!scalar @ARGV) # no command line arguments supplied
{
	print STDOUT "\nDetects potential cross-contamination.\n";
	print STDOUT "Usage: polyphonia detect_cross_contam [options]\n";
	print STDOUT "\n";
	
	print STDOUT "OPTIONS:\n";
	print STDOUT "- Reference (required):\n";
	print STDOUT "\t-f | --ref FILE\t\t\tReference fasta file [null]\n";
	print STDOUT "\n";
	
	print STDOUT "- Consensus genomes (aligned or not aligned, not both; at least one file required):\n";
	print STDOUT "\t-c | --consensus FILE(S)\tUnaligned consensus genome or genomes [null]\n";
	print STDOUT "\t-a | --consensus-aligned FILE\tConsensus genomes pre-aligned to reference as fasta alignment; reference provided by --ref must be first [null]\n";
	print STDOUT "\t-g | --min-covered FLOAT\tMinimum proportion genome covered for a sample to be included [".$DEFAULT_MINIMUM_GENOME_COVERAGE."]\n";
	print STDOUT "\n";
	
	print STDOUT "- Within-sample diversity (any combination; at least one file required):\n";
	print STDOUT "\t-b | --bam FILE(S)\t\tAligned and trimmed reads as bam file(s); must use reference provided by --ref [null]\n";
	print STDOUT "\t-v | --vcf FILE(S)\t\tVCF file(s) output by LoFreq or GATK; must use reference provided by --ref [null]\n";
	print STDOUT "\t-h | --het FILE(S)\t\tTab-separated heterozygosity summary tables; see documentation for format [null]\n";
	
	print STDOUT "\t-e | --min-readcount INT\tMinimum minor allele readcount for position to be considered heterozygous [".$DEFAULT_MINIMUM_MINOR_ALLELE_READCOUNT."]\n";
	print STDOUT "\t-i | --min-maf FLOAT\t\tMinimum minor allele frequency for position to be considered heterozygous [".$DEFAULT_MINIMUM_MINOR_ALLELE_FREQUENCY."]\n";
	print STDOUT "\n";
	
	print STDOUT "- Plate map and neighbors (any combination, all optional):\n";
	print STDOUT "\t-m | --map FILE(S)\t\tOptional plate map (tab-separated, no header: sample name, plate position (e.g., A8)); provides substantial speed-up [null]\n";
	print STDOUT "\t-n | --direct BOOL\t\tCompare direct plate neighbors (left, right, top, bottom) [".int_to_bool_string($DEFAULT_COMPARE_DIRECT_NEIGHBORS)."]\n";
	print STDOUT "\t-d | --diagonal BOOL\t\tCompare diagonal plate neighbors (top-right, bottom-right, top-left, bottom-left) [".int_to_bool_string($DEFAULT_COMPARE_DIAGONAL_NEIGHBORS)."]\n";
	print STDOUT "\t-w | --row BOOL\t\t\tCompare samples in the same row (e.g., row A) [".int_to_bool_string($DEFAULT_COMPARE_ROW)."]\n";
	print STDOUT "\t-l | --column BOOL\t\tCompare samples in the same column (e.g., column 8) [".int_to_bool_string($DEFAULT_COMPARE_COLUMN)."]\n";
	print STDOUT "\t-t | --plate BOOL\t\tCompare all samples in the same plate map [".int_to_bool_string($DEFAULT_COMPARE_WHOLE_PLATE_MAP)."]\n";
	print STDOUT "\n";
	
	print STDOUT "- Misc:\n";
	print STDOUT "\t-a | --max-mismatches INT\tMaximum allowed bases in contaminating sample consensus not matching contaminated sample alleles [".$DEFAULT_MAXIMUM_ALLOWED_MISMATCHES."]\n";
	print STDOUT "\t-p | --cores INT\t\tOptional number of cores to use for preprocessing in parallel [".$DEFAULT_CORES_TO_USE."]\n";
	print STDOUT "\t-u | --verbose BOOL\t\tPrint progress to STDOUT [".int_to_bool_string($DEFAULT_VERBOSE)."]\n";
	print STDOUT "\t-i | --directory DIRECTORY\tPath of directory to store intermediate and temp files [".$default_temp_intermediate_files_directory."]\n";
	print STDOUT "\t-o | --output FILE\t\tOutput file path [".$default_output_file."]\n";
	print STDOUT "\t-j | --overwrite FILE\t\tOverwrite output, intermediate, and temp files at input paths [".int_to_bool_string($DEFAULT_OVERWRITE)."]\n";
	print STDOUT "\n\n";
	exit;
}


# command line arguments set to default values
my $reference_genome_file = "";
my @consensus_genome_files = ();
my $consensus_genomes_aligned_file = "";
my $minimum_genome_coverage = $DEFAULT_MINIMUM_GENOME_COVERAGE;
my $maximum_allowed_mismatches = $DEFAULT_MAXIMUM_ALLOWED_MISMATCHES;

my @aligned_and_trimmed_bam_files = ();
my @vcf_files = ();
my @heterozygosity_tables = ();

my $temp_intermediate_directory = $default_temp_intermediate_files_directory;
my $output_file_path = $default_output_file;

my $cores_to_use = $DEFAULT_CORES_TO_USE;
my $overwrite = $DEFAULT_OVERWRITE;
my $verbose = $DEFAULT_VERBOSE;

my $minimum_minor_allele_readcount = $DEFAULT_MINIMUM_MINOR_ALLELE_READCOUNT;
my $minimum_minor_allele_frequency = $DEFAULT_MINIMUM_MINOR_ALLELE_FREQUENCY;

my @plate_map_files = ();
my $compare_direct_neighbors = $DEFAULT_COMPARE_DIRECT_NEIGHBORS;
my $compare_diagonal_neighbors = $DEFAULT_COMPARE_DIAGONAL_NEIGHBORS;
my $compare_row = $DEFAULT_COMPARE_ROW;
my $compare_column = $DEFAULT_COMPARE_COLUMN;
my $compare_whole_plate_map = $DEFAULT_COMPARE_WHOLE_PLATE_MAP;

# reads in command line arguments
my $argument;
my $argument_index;
for($argument_index = 0; $argument_index <= $#ARGV; $argument_index++)
{
	$argument = $ARGV[$argument_index];
	
	my $input;
	if(($input = read_in_input_file_argument("-f", "--ref")) ne "-1")
	{
		$reference_genome_file = $input;
	}
	elsif(($input = read_in_input_file_argument("-o", "--output")) ne "-1")
	{
		$output_file_path = $input;
	}
	elsif(($input = read_in_input_file_argument("-i", "--directory")) ne "-1")
	{
		$temp_intermediate_directory = $input;
	}
	elsif(($input = read_in_boolean_argument("-j", "--overwrite")) != -1)
	{
		$overwrite = $input;
	}
	elsif(($input = read_in_positive_integer_argument("-a", "--max-mismatches")) != -1)
	{
		$maximum_allowed_mismatches = $input;
	}
	elsif(($input = read_in_positive_integer_argument("-p", "--cores")) != -1)
	{
		$cores_to_use = $input;
	}
	elsif(($input = read_in_input_file_argument("-a", "--consensus-aligned")) ne "-1")
	{
		$consensus_genomes_aligned_file = $input;
	}
	elsif(($input = read_in_input_files_argument("-c", "--consensus")) ne "-1")
	{
		push(@consensus_genome_files, @$input);
	}
	elsif(($input = read_in_positive_float_argument("-g", "--min-covered")) != -1)
	{
		$minimum_genome_coverage = $input;
	}
	elsif(($input = read_in_input_files_argument("-b", "--bam")) ne "-1")
	{
		push(@aligned_and_trimmed_bam_files, @$input);
	}
	elsif(($input = read_in_input_files_argument("-v", "--vcf")) ne "-1")
	{
		push(@vcf_files, @$input);
	}
	elsif(($input = read_in_input_files_argument("-h", "--het")) ne "-1")
	{
		push(@heterozygosity_tables, @$input);
	}
	elsif(($input = read_in_positive_integer_argument("-e", "--min-readcount")) != -1)
	{
		$minimum_minor_allele_readcount = $input;
	}
	elsif(($input = read_in_positive_float_argument("-i", "--min-maf")) != -1)
	{
		$minimum_minor_allele_frequency = $input;
	}
	elsif(($input = read_in_input_files_argument("-m", "--map")) ne "-1")
	{
		push(@plate_map_files, @$input);
	}
	elsif(($input = read_in_boolean_argument("-n", "--direct")) != -1)
	{
		$compare_direct_neighbors = $input;
	}
	elsif(($input = read_in_boolean_argument("-d", "--diagonal")) != -1)
	{
		$compare_diagonal_neighbors = $input;
	}
	elsif(($input = read_in_boolean_argument("-w", "--row")) != -1)
	{
		$compare_row = $input;
	}
	elsif(($input = read_in_boolean_argument("-l", "--column")) != -1)
	{
		$compare_column = $input;
	}
	elsif(($input = read_in_boolean_argument("-t", "--plate")) != -1)
	{
		$compare_whole_plate_map = $input;
	}
	elsif(($input = read_in_boolean_argument("-u", "--verbose")) != -1)
	{
		$verbose = $input;
	}
	else
	{
		print STDERR "Error: argument ".$argument." not recognized. Exiting.\n";
		die;
	}
}


# verifies that all necessary command line arguments have been included
if(!$reference_genome_file)
{
	print STDERR "Error: reference genome not provided. Exiting.\n";
	die;
}
if(!scalar @consensus_genome_files and $consensus_genomes_aligned_file)
{
	print STDERR "Error: no consensus genome files provided. Exiting.\n";
	die;
}
if(scalar @consensus_genome_files and $consensus_genomes_aligned_file)
{
	print STDERR "Error: consensus genomes must be either aligned or not aligned, not both. Exiting.\n";
	die;
}
if(!scalar @vcf_files and !scalar @aligned_and_trimmed_bam_files and !scalar @heterozygosity_tables)
{
	print STDERR "Error: no minor allele files provided. Exiting.\n";
	die;
}
if(!scalar @plate_map_files and !$compare_direct_neighbors and !$compare_diagonal_neighbors
	and !$compare_row and !$compare_column and !$compare_whole_plate_map)
{
	print STDERR "Warning: plate map(s) supplied but plate map options all set to False. "
		."Plate map(s) will not be used.\n";
	@plate_map_files = ();
}
if($minimum_genome_coverage < 0 or $minimum_genome_coverage > 1)
{
	print STDERR "Error: minimum genome coverage is not between 0 and 1. Exiting.\n";
	die;
}


# verifies that all input files exist
if($reference_genome_file and !-e $reference_genome_file)
{
	print STDERR "Error: reference genome file does not exist:\n\t"
		.$reference_genome_file."\nExiting.\n";
	die;
}
foreach my $consensus_genome_file(@consensus_genome_files)
{
	if($consensus_genome_file and !-e $consensus_genome_file)
	{
		print STDERR "Error: consensus genome file does not exist:\n\t"
			.$consensus_genome_file."\nExiting.\n";
		die;
	}
}
if($consensus_genomes_aligned_file and !-e $reference_genome_file)
{
	print STDERR "Error: aligned consensus genomes file does not exist:\n\t"
		.$consensus_genomes_aligned_file."\nExiting.\n";
	die;
}
foreach my $aligned_and_trimmed_bam_file(@aligned_and_trimmed_bam_files)
{
	if($aligned_and_trimmed_bam_file and !-e $aligned_and_trimmed_bam_file)
	{
		print STDERR "Error: aligned and trimmed bam file does not exist:\n\t"
			.$aligned_and_trimmed_bam_file."\nExiting.\n";
		die;
	}
}
foreach my $vcf_file(@vcf_files)
{
	if($vcf_file and !-e $vcf_file)
	{
		print STDERR "Error: vcf file does not exist:\n\t".$vcf_file."\nExiting.\n";
		die;
	}
}
foreach my $heterozygosity_table(@heterozygosity_tables)
{
	if($heterozygosity_table and !-e $heterozygosity_table)
	{
		print STDERR "Error: heterozygosity table does not exist:\n\t"
			.$heterozygosity_table."\nExiting.\n";
		die;
	}
}
foreach my $plate_map_file(@plate_map_files)
{
	if($plate_map_file and !-e $plate_map_file)
	{
		print STDERR "Error: plate map file does not exist:\n\t"
			.$plate_map_file."\nExiting.\n";
		die;
	}
}


# prepares directory for temporary and intermediate files
# adds / to end of directory path if it isn't already there
if($temp_intermediate_directory !~ /\/$/) # if doesn't end in /
{
	$temp_intermediate_directory .= "/"; # adds / to the end
}

# creates directory for temporary and intermediate files if it doesn't already exist
if(-e $temp_intermediate_directory and -d $temp_intermediate_directory)
{
	# directory already exists
	print STDERR "Warning: adding temporary and intermediate files to already existing directory:\n\t"
		.$temp_intermediate_directory."\n";
}
elsif(-e $temp_intermediate_directory)
{
	# directory already exists and is a file
	print STDERR "Error: temporary intermediate directory already exists and is a file:\n\t"
		.$temp_intermediate_directory."\nExiting.\n";
	die;
}
else
{
	# directory doesn't exist
	# create directory and all necessary parent directories
	`mkdir -p $temp_intermediate_directory`;
}


# prints input files and options entered
# reference
if($verbose)
{
	print STDOUT "\n";
	print STDOUT "REFERENCE:\n\t".$reference_genome_file."\n";
}
# print STDOUT "\n" if $verbose;
# print STDOUT "REFERENCE:\n\t".$reference_genome_file."\n" if $verbose;

# consensus genome files
print STDOUT "CONSENSUS GENOMES:\n" if $verbose;
foreach my $consensus_genome_file(@consensus_genome_files)
{
	print STDOUT "\t".$consensus_genome_file."\n" if $verbose;
}
if($consensus_genomes_aligned_file)
{
	print STDOUT "\tpre-aligned: ".$consensus_genomes_aligned_file."\n" if $verbose;
}

# within-sample diversity files
print STDOUT "WITHIN-SAMPLE DIVERSITY:\n" if $verbose;
foreach my $aligned_and_trimmed_bam_file(@aligned_and_trimmed_bam_files)
{
	print STDOUT "\t".$aligned_and_trimmed_bam_file."\n" if $verbose;
}
foreach my $vcf_file(@vcf_files)
{
	print STDOUT "\tpre-processed vcf files: ".$vcf_file."\n" if $verbose;
}
foreach my $heterozygosity_table(@heterozygosity_tables)
{
	print STDOUT "\tfully pre-processed heterozygosity tables: ".$heterozygosity_table."\n" if $verbose;
}

# optional plate map file(s) and related options
if(scalar @plate_map_files)
{
	print STDOUT "PLATE MAP(S):\n" if $verbose;
	foreach my $plate_map_file(@plate_map_files)
	{
		print STDOUT "\t".$plate_map_file."\n" if $verbose;
	}
	
	print STDOUT "PLATE MAP USE:\n" if $verbose;
	if($compare_whole_plate_map)
	{
		print STDOUT "Comparing all samples in the same plate map.\n" if $verbose;
	}
	else
	{
		if($compare_direct_neighbors)
		{
			print STDOUT "\tComparing direct plate neighbors (left, right, top, bottom).\n" if $verbose;
		}
		if($compare_diagonal_neighbors)
		{
			print STDOUT "\tComparing diagonal plate neighbors (top-right, bottom-right, top-left, bottom-left).\n" if $verbose;
		}
		if($compare_row)
		{
			print STDOUT "\tComparing samples in the same row (e.g., row A).\n" if $verbose;
		}
		if($compare_column)
		{
			print STDOUT "\tComparing samples in the same column (e.g., column 8).\n" if $verbose;
		}
	}
}
print STDOUT "\n" if $verbose;


# retrieves sample names from plate maps if possible
my %sample_names = (); # key: sample name -> value: 1
if(scalar @plate_map_files)
{
	print STDOUT "retrieving sample names from plate map(s)...\n" if $verbose;
	foreach my $plate_map_file(@plate_map_files)
	{
		open PLATE_MAP, "<$plate_map_file" || die "Could not open $plate_map_file to read; terminating =(\n";
		while(<PLATE_MAP>) # for each line in the file
		{
			chomp;
			my $line = $_;
			if($line =~ /\S/) # non-empty line
			{
				my @items = split($DELIMITER, $line);
				my $sample_name = $items[$PLATE_MAP_SAMPLE_COLUMN];
				my $position = uc $items[$PLATE_MAP_POSITION_COLUMN];
				
				if(is_valid_plate_position($position))
				{
					$sample_names{$sample_name} = 1;
				}
			}
		}
		close PLATE_MAP;
	}
	
	# prints number of samples remaining
	print STDOUT (keys %sample_names)." samples...\n" if $verbose;
	
	# catalogues consensus genome sample names
	print STDOUT "removing sample names without associated consensus genome...\n" if $verbose;
	my %sample_name_has_consensus_genome = (); # key: sample name -> value: 1 if sample has associated consensus genome
	foreach my $consensus_genome_fasta_file(@consensus_genome_files, $consensus_genomes_aligned_file)
	{
		open FASTA_FILE, "<$consensus_genome_fasta_file" || die "Could not open $consensus_genome_fasta_file to read; terminating =(\n";
		while(<FASTA_FILE>) # for each line in the file
		{
			chomp;
			if($_ =~ /^>(.*)/) # header line
			{
				my $sample_name = $1;
				$sample_name_has_consensus_genome{$sample_name} = 1;
			}
		}
		close FASTA_FILE;
	}
	
	# prints number of samples remaining
	print STDOUT (keys %sample_names)." samples remain...\n" if $verbose;
	
	# removes any sample names that don't have a consensus genome
	foreach my $sample_name(keys %sample_names)
	{
		if(!$sample_name_has_consensus_genome{$sample_name})
		{
			delete $sample_names{$sample_name};
		}
	}
}

# if no plate map, retrieves sample names from consensus genome fasta files
else
{
	print STDOUT "retrieving sample names from consensus genome fasta file(s)...\n" if $verbose;
	foreach my $consensus_genome_fasta_file(@consensus_genome_files, $consensus_genomes_aligned_file)
	{
		open FASTA_FILE, "<$consensus_genome_fasta_file" || die "Could not open $consensus_genome_fasta_file to read; terminating =(\n";
		while(<FASTA_FILE>) # for each line in the file
		{
			chomp;
			if($_ =~ /^>(.*)/) # header line
			{
				my $sample_name = $1;
				$sample_names{$sample_name} = 1;
			}
		}
		close FASTA_FILE;
	}
	
	# prints number of samples remaining
	print STDOUT (keys %sample_names)." samples...\n" if $verbose;
}

# records file stage for each within-sample diversity file
print STDOUT "recording stage of each within-sample diversity file...\n" if $verbose;
my %within_sample_diversity_file_to_stage = (); # key: file path -> value: "bam" or "vcf" or "het" (processed heterozygosity table)
foreach my $file_path(@heterozygosity_tables)
{
	$within_sample_diversity_file_to_stage{$file_path} = "het";
}
foreach my $file_path(@vcf_files)
{
	$within_sample_diversity_file_to_stage{$file_path} = "vcf";
}
foreach my $file_path(@aligned_and_trimmed_bam_files)
{
	$within_sample_diversity_file_to_stage{$file_path} = "bam";
}

# retrieves within-sample diversity file for each sample
print STDOUT "retrieving within-sample diversity file for each sample...\n" if $verbose;
# removes any sample names that don't have an associated within-sample diversity file
my %sample_name_to_within_sample_diversity_file = (); # key: sample name -> value: file path of within-sample diversity file (bam or vcf or heterozygosity table)
foreach my $file_path(@aligned_and_trimmed_bam_files, @vcf_files, @heterozygosity_tables) # most processed is looked at last (so that most processed replaces least processed)
{
	# trims file path to file name
	my $potential_sample_name = retrieve_file_name($file_path);

	# retrieves largest possible sample name that collides with a sample name
	# (file name sample1.ext1.ext2 has possible sample names sample1.ext1.ext2, sample1.ext1, sample1)
	my $sample_name_found = 0;
	while($potential_sample_name and !$sample_name_found)
	{
		if($sample_names{$potential_sample_name})
		{
			# potential sample name collides with a sample name from consensus genome or plate map files
			# this is our sample name
			$sample_name_to_within_sample_diversity_file{$potential_sample_name} = $file_path;
			$sample_name_found = 1;
		}
		else
		{
			$potential_sample_name = trim_off_file_extension($potential_sample_name);
		}
	}
}

# removes any sample names that don't have a within-sample diversity file
print STDOUT "removing sample names without within-sample diversity file...\n" if $verbose;
foreach my $sample_name(keys %sample_names)
{
	if(!$sample_name_to_within_sample_diversity_file{$sample_name})
	{
		delete $sample_names{$sample_name};
	}
}

# prints number of samples remaining
print STDOUT (keys %sample_names)." samples remain...\n" if $verbose;


# if a plate map is provided, removes any samples that do not have neighbors
# reads in plate map positions of all samples
my %sample_name_to_all_plate_positions = (); # key: sample name -> value: string including all plate positions the sample appears in
my %sample_name_to_all_plates = (); # key: sample name -> value: string including all plates the sample appears in

my %plate_position_to_sample_name = (); # key: plate position -> value: sample name
my %sample_name_to_plate_position = (); # key: sample name -> value: plate position

if(scalar @plate_map_files)
{
	print STDOUT "reading in plate map positions and removing sample names without plate neighbors...\n" if $verbose;
	my %sample_has_plate_neighbors = (); # key: sample name -> value: 1 if sample has plate neighbors
	foreach my $plate_map_file(@plate_map_files)
	{
		# reads in plate map
		open PLATE_MAP, "<$plate_map_file" || die "Could not open $plate_map_file to read; terminating =(\n";
		while(<PLATE_MAP>) # for each line in the file
		{
			chomp;
			my $line = $_;
			if($line =~ /\S/) # non-empty line
			{
				my @items = split($DELIMITER, $line);
				my $sample_name = $items[$PLATE_MAP_SAMPLE_COLUMN];
				my $position = uc $items[$PLATE_MAP_POSITION_COLUMN];
				
				if($position and $sample_name and $sample_names{$sample_name})
				{
					$plate_position_to_sample_name{$position} = $sample_name;
					$sample_name_to_plate_position{$sample_name} = $position;
					
					# adds plate map position to string
					if($sample_name_to_all_plate_positions{$sample_name})
					{
						$sample_name_to_all_plate_positions{$sample_name} .= $PLATE_MAP_POSITION_LIST_SEPARATOR;
					}
					$sample_name_to_all_plate_positions{$sample_name} .= $position;
					
					# adds corresponding plate map file
					if(scalar @plate_map_files > 1)
					{
						if($sample_name_to_all_plates{$sample_name})
						{
							$sample_name_to_all_plates{$sample_name} .= $FILE_LIST_SEPARATOR;
						}
						$sample_name_to_all_plates{$sample_name} .= $plate_map_file;
					}
				}
			}
		}
		close PLATE_MAP;
		
		# checks if each sample has at least one neighbor
		foreach my $plate_position(keys %plate_position_to_sample_name)
		{
			my $sample_name = $plate_position_to_sample_name{$plate_position};
			my @neighboring_samples = retrieve_samples_neighboring_plate_position($plate_position, $sample_name);
			
			if(scalar @neighboring_samples)
			{
				$sample_has_plate_neighbors{$sample_name} = 1;
			}
		}
	}
	
	# removes any sample names that don't have at least one plate neighbor
	foreach my $sample_name(keys %sample_names)
	{
		if(!$sample_has_plate_neighbors{$sample_name})
		{
			delete $sample_names{$sample_name};
		}
	}
}

# prints number of samples remaining
print STDOUT (keys %sample_names)." samples remain...\n" if $verbose;

# verifies that we still have samples to compare
print STDOUT "verifying that there are samples to compare...\n" if $verbose;
if(scalar keys %sample_names < 2)
{
	print STDERR "Error: no pairs of samples to compare. Exiting.\n";
	die;
}


# aligns consensus genomes if they aren't already aligned
if(!$consensus_genomes_aligned_file)
{
	# makes temp file with just a newline
	my $newline_temp_file = $temp_intermediate_directory."newline_temp.fasta";
	check_if_file_exists_before_writing($newline_temp_file);
	`echo "" > $newline_temp_file`;

	# puts together command to concatenate reference genome and all consensus genome fasta files
	my $all_consensus_genomes_file = $temp_intermediate_directory."all_consensus_genomes_concat.fasta";
	check_if_file_exists_before_writing($all_consensus_genomes_file);
	my $cat_command = "cat ".$reference_genome_file." ";
	foreach my $consensus_genome_file(@consensus_genome_files)
	{
		$cat_command .= $consensus_genome_file." $newline_temp_file ";
	}
	$cat_command .= "> ".$all_consensus_genomes_file;
	# cat k.txt <(echo) h.txt > new.txt
	
	# concatenates reference genome and all consensus genome fasta files
	`$cat_command`;
	
	# aligns all consensus genomes
	$consensus_genomes_aligned_file = $temp_intermediate_directory."all_consensus_genomes_MAFFT_aligned.fasta";
	check_if_file_exists_before_writing($consensus_genomes_aligned_file);
	`$MAFFT_EXECUTABLE_FILE_PATH $all_consensus_genomes_file > $consensus_genomes_aligned_file`;
	
	# removes temp files
	`rm $newline_temp_file`;
	`rm $all_consensus_genomes_file`;
}


# pre-processes within-sample diversity files
# parallelization based on https://perlmaven.com/speed-up-calculation-by-running-in-parallel
print STDOUT "pre-processing within-sample diversity files, using ".$cores_to_use." cores in parallel...\n" if $verbose;

my %updated_within_sample_diversity_files = ();
my $pm = Parallel::ForkManager -> new($cores_to_use);
$pm -> run_on_finish(
sub
{
	my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
	my $q = $data_structure_reference -> {input};
	$updated_within_sample_diversity_files{$q} = $data_structure_reference -> {updated_file};
});

foreach my $sample(keys %sample_names)
{
	my $pid = $pm -> start and next;
	my $res = process_within_sample_diversity_file_for_sample($sample);
	$pm -> finish(0, {updated_file => $res, input => $sample});
}
$pm -> wait_all_children;

# saves updated within-sample diversity files
foreach my $sample(keys %updated_within_sample_diversity_files)
{
	my $within_sample_diversity_file = $updated_within_sample_diversity_files{$sample};
	
	$sample_name_to_within_sample_diversity_file{$sample} = $within_sample_diversity_file;
	$within_sample_diversity_file_to_stage{$within_sample_diversity_file} = "het";
}


# reads in aligned consensus genomes
print STDERR "reading in aligned consensus genomes...\n";
open ALIGNED_CONSENSUS_GENOMES, "<$consensus_genomes_aligned_file" || die "Could not open $consensus_genomes_aligned_file to read; terminating =(\n";
my %sequence_name_to_consensus = (); # key: sequence name -> value: consensus sequence, including gaps froms alignment
my $reference_sequence = ""; # first sequence in alignment

my $sequence = "";
my $sample_name = "";
while(<ALIGNED_CONSENSUS_GENOMES>) # for each line in the file
{
	chomp;
	if($_ =~ /^>(.*)/) # header line
	{
		# process previous sequence
		$sequence = uc($sequence);
		if($sequence and $sample_name and $sample_names{$sample_name})
		{
			$sequence_name_to_consensus{$sample_name} = $sequence;
		}
		if(!$reference_sequence) # reference sequence is first sequence in alignment
		{
			$reference_sequence = $sequence;
		}
	
		# prepare for next sequence
		$sequence = "";
		$sample_name = $1;
	}
	else
	{
		$sequence .= $_;
	}
}
# process final sequence
if($sequence and $sample_name and $sample_names{$sample_name})
{
	$sequence_name_to_consensus{$sample_name} = uc($sequence);
}
close ALIGNED_CONSENSUS_GENOMES;


# counts unambiguous bases in reference
my $reference_sequence_length = count_unambiguous_bases_in_sequence(split(//, $reference_sequence));
if(!$reference_sequence_length)
{
	print STDERR "Error: reference sequence in consensus genome alignment contains no "
		."unambiguous (A, T, C, G) bases:\n\t".$consensus_genomes_aligned_file
		."\nExiting.\n";
	die;
}


# header line for plate-specific output file
my $plate_header_line = "";
$plate_header_line .= "well	sample".$DELIMITER;
$plate_header_line .= "contamination_source_well".$DELIMITER;
$plate_header_line .= "contamination_source_sample".$DELIMITER;
$plate_header_line .= "estimated_contamination_volume";


# prepares to process sample pairs in parallel
# parallelization based on https://perlmaven.com/speed-up-calculation-by-running-in-parallel
my %results = ();
my %plate_results = ();
$pm = Parallel::ForkManager -> new($cores_to_use);
$pm -> run_on_finish(
sub
{
	my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
	my $q = $data_structure_reference -> {input};
	$results{$q} = $data_structure_reference -> {result};
	$plate_results{$q} = $data_structure_reference -> {plate_result};
});

# if plate map(s) provided, compares all neighboring samples
if(scalar @plate_map_files)
{
	print STDOUT "comparing all neighboring samples...\n" if $verbose;
	my %sample_pair_compared = (); # keys: sample name, sample name -> value: 1 if sample pair compared
	foreach my $plate_map_file(@plate_map_files)
	{
		# reads in plate map
		%plate_position_to_sample_name = (); # key: position -> value: sample name
		open PLATE_MAP, "<$plate_map_file" || die "Could not open $plate_map_file to read; terminating =(\n";
		while(<PLATE_MAP>) # for each line in the file
		{
			chomp;
			my $line = $_;
			if($line =~ /\S/) # non-empty line
			{
				my @items = split($DELIMITER, $line);
				my $sample_name = $items[$PLATE_MAP_SAMPLE_COLUMN];
				my $plate_position = uc $items[$PLATE_MAP_POSITION_COLUMN];
				
				if($plate_position and $sample_name and $sample_names{$sample_name})
				{
					$plate_position_to_sample_name{$plate_position} = $sample_name;
				}
			}
		}
		close PLATE_MAP;
		
		# verifies that all plate positions are consistent with this plate being a 96-well plate
		my $plate_is_96_well = 1;
		foreach my $plate_position(keys %plate_position_to_sample_name)
		{
			if($plate_position =~ /^([A-Z])\s*(\d+)$/)
			{
				my $letter = $1;
				my $number = $2;
				
				if($number < 1 or $number > 12)
				{
					$plate_is_96_well = 0;
				}
				if(!($letter eq "A" or $letter eq "B" or $letter eq "C" or $letter eq "D"
					or $letter eq "E" or $letter eq "F" or $letter eq "G" or $letter eq "H"))
				{
					$plate_is_96_well = 0;
				}
			}
		}
		
		# compares all pairs of samples that are neighbors
		foreach my $plate_position(keys %plate_position_to_sample_name)
		{
			my $sample_1 = $plate_position_to_sample_name{$plate_position};
			my @neighboring_samples = retrieve_samples_neighboring_plate_position($plate_position, $sample_1);
			
			foreach my $sample_2(@neighboring_samples)
			{
				my $neighboring_plate_position = $sample_name_to_plate_position{$sample_2};
				if(!$sample_pair_compared{$sample_1}{$sample_2}) # this pair of samples not already compared
				{
					# records that we are comparing this pair (to avoid comparing it again)
					$sample_pair_compared{$sample_1}{$sample_2} = 1;
					$sample_pair_compared{$sample_2}{$sample_1} = 1;
				
					# tests sample 1 contaminating sample 2 and vice versa
					my $pid = $pm -> start and next;
					my ($res, $plate_res) = detect_potential_contamination_in_sample_pair_both_directions($sample_1, $sample_2);
					$pm -> finish(0, {result => $res, plate_result => $plate_res, input => $sample_1."__".$sample_2});
				}
			}
		}
		$pm -> wait_all_children;
		
		# prints plate map for this particular plate map
		if($plate_is_96_well and scalar keys %plate_results)
		{
			my $plate_output_file = $temp_intermediate_directory.retrieve_file_name($plate_map_file)."_potential_cross_contamination.txt";
			check_if_file_exists_before_writing($plate_output_file);
			
			open PLATE_OUT_FILE, ">$plate_output_file" || die "Could not open $plate_output_file to write; terminating =(\n";
			print PLATE_OUT_FILE $plate_header_line.$NEWLINE;

			# prints output line for each pair of samples with potential contamination
			foreach my $output_line(values %plate_results)
			{
				if($output_line)
				{
					print PLATE_OUT_FILE $output_line.$NEWLINE;
				}
			}

			# closes output file
			close PLATE_OUT_FILE;
			
			# generates visualization
			my $plate_visualization_output_file = trim_off_file_extension($plate_output_file)."_visualization";
			check_if_file_exists_before_writing($plate_visualization_output_file.".jpg");
			check_if_file_exists_before_writing($plate_visualization_output_file.".pdf");
			
			exec("Rscript $PLATE_VISUALIZATION_FILE_PATH $plate_output_file $plate_visualization_output_file");
		}
		
		# clears plate output for next plate map
		my %plate_results = ();
	}
}

# if no plate map provided, compares all pairs of samples
else
{
	print STDOUT "comparing all pairs of samples...\n" if $verbose;
	my @sample_names_array = keys %sample_names;
	for(my $index_1 = 0; $index_1 <= $#sample_names_array; $index_1++)
	{
		my $sample_1 = $sample_names_array[$index_1];
		for(my $index_2 = $index_1+1; $index_2 <= $#sample_names_array; $index_2++)
		{
			my $sample_2 = $sample_names_array[$index_2];
			
			# tests sample 1 contaminating sample 2 and vice versa
			my $pid = $pm -> start and next;
			my ($res, $plate_res) = detect_potential_contamination_in_sample_pair_both_directions($sample_1, $sample_2);
			$pm -> finish(0, {result => $res, input => $sample_1."__".$sample_2});
		}
	}
	$pm -> wait_all_children;
}


# header line for output
my $header_line = "";
$header_line .= "potential_contaminated_sample".$DELIMITER;
$header_line .= "potential_contaminated_sample_unambiguous_bases".$DELIMITER;
$header_line .= "potential_contaminated_sample_genome_covered".$DELIMITER;
$header_line .= "num_positions_with_heterozygosity".$DELIMITER;
$header_line .= "alleles_at_positions_with_heterozygosity".$DELIMITER;

$header_line .= "potential_contaminating_sample".$DELIMITER;
$header_line .= "potential_contaminating_sample_unambiguous_bases".$DELIMITER;
$header_line .= "potential_contaminating_sample_genome_covered".$DELIMITER;
$header_line .= "minor_alleles_matched".$DELIMITER;
$header_line .= "major_alleles_matched".$DELIMITER;
$header_line .= "heterozygous_positions_matched".$DELIMITER;
$header_line .= "alleles_matched".$DELIMITER;
$header_line .= "num_mismatches".$DELIMITER;
$header_line .= "mismatches".$DELIMITER;

$header_line .= "appearance_of_potential_contamination".$DELIMITER;
$header_line .= "estimated_contamination_volume".$DELIMITER;
$header_line .= "contaminating_allele_frequency_range".$DELIMITER;
$header_line .= "contaminating_allele_frequencies";

if(scalar @plate_map_files)
{
	$header_line .= $DELIMITER;
	$header_line .= "potential_contaminated_sample_plate_position".$DELIMITER;
	$header_line .= "potential_contaminating_sample_plate_position";
	if(scalar @plate_map_files > 1)
	{
		$header_line .= $DELIMITER;
		$header_line .= "potential_contaminated_sample_plate";
		$header_line .= $DELIMITER;
		$header_line .= "potential_contaminating_sample_plate";
	}
}

# creates directory for output file if it does not already exist
my $output_file_directory = retrieve_file_directory($output_file_path);
if(!-d $output_file_directory and -e $output_file_directory)
{
	# directory already exists and is a file
	print STDERR "Error: output file directory is a file:\n\t"
		.$output_file_directory."\nExiting.\n";
	die;
}
elsif(!-d $output_file_directory)
{
	# directory doesn't exist
	# create directory and all necessary parent directories
	`mkdir -p $output_file_directory`;
}

# opens output file and prints header line to output file or to console
check_if_file_exists_before_writing($output_file_path);
open OUT_FILE, ">$output_file_path" || die "Could not open $output_file_path to write; terminating =(\n";
print OUT_FILE $header_line.$NEWLINE;

# prints output line for each pair of samples with potential contamination
foreach my $output_line(values %results)
{
	if($output_line)
	{
		print OUT_FILE $output_line.$NEWLINE;
	}
}

# closes output file
close OUT_FILE;


# HELPER FUNCTIONS FOR DETECTING POTENTIAL CONTAMINATION

# compares these two samples to check for contamination in both directions
sub detect_potential_contamination_in_sample_pair_both_directions
{
	my $sample_1 = $_[0];
	my $sample_2 = $_[1];
	
	my ($result_1, $result_1_plate) = detect_potential_contamination_in_sample_pair($sample_1, $sample_2);
	my ($result_2, $result_2_plate) = detect_potential_contamination_in_sample_pair($sample_2, $sample_1);
	
	# returns output lines
	if($result_1 and $result_2)
	{
		return ($result_1.$NEWLINE.$result_2, $result_1_plate.$NEWLINE.$result_2_plate);
	}
	elsif($result_1)
	{
		return ($result_1, $result_1_plate);
	}
	elsif($result_2)
	{
		return ($result_2, $result_2_plate);
	}
}

# prepares and retrieves files for the two input samples; checks if the second sample
# potentially could have contaminated the first
sub detect_potential_contamination_in_sample_pair
{
	my $potential_contaminated_sample = $_[0];
	my $potential_contaminating_sample = $_[1];
	
	# if needed, processes within-sample diversity file for potential contaminated sample
	my $within_sample_diversity_file = process_within_sample_diversity_file_for_sample($potential_contaminated_sample);
	$sample_name_to_within_sample_diversity_file{$potential_contaminated_sample} = $within_sample_diversity_file;
	$within_sample_diversity_file_to_stage{$within_sample_diversity_file} = "het";
	
	# retrieves within-sample diversity file for potential contaminated sample
	my $potential_contaminated_within_sample_diversity_file = $sample_name_to_within_sample_diversity_file{$potential_contaminated_sample};
	
	# reads in heterozygosity table generated for within-sample diversity of potential contaminated sample
	my %minor_alleles = (); # key: position, allele -> value: 1 if a sample has >=10 reads and >=0.3 MAF at this minor allele at this position
	my %major_alleles = (); # key: position, allele -> value: 1
	my %minor_allele_frequencies = (); # key: position, allele -> value: frequency of minor allele
	my %major_allele_frequencies = (); # key: position, allele -> value: frequency of major allele
	my %positions_with_heterozygosity = (); # key: position -> value: 1 if position has minor allele
	my $number_positions_with_heterozygosity = 0;
	my $list_of_alleles_to_print = "";
	open HETEROZYGOSITY_TABLE, "<$potential_contaminated_within_sample_diversity_file" || die "Could not open $potential_contaminated_within_sample_diversity_file to read; terminating =(\n";
	while(<HETEROZYGOSITY_TABLE>) # for each line in the file
	{
		chomp;
		my $line = $_;
		
		# parses this line
		my @items = split($DELIMITER, $line);
		my $position = $items[$HETEROZYGOSITY_TABLE_POSITION_COLUMN];
		my $minor_allele = $items[$HETEROZYGOSITY_TABLE_MINOR_ALLELE_COLUMN];
		my $major_allele = $items[$HETEROZYGOSITY_TABLE_MAJOR_ALLELE_COLUMN];
		my $minor_allele_readcount = $items[$HETEROZYGOSITY_TABLE_MINOR_ALLELE_READCOUNT_COLUMN];
		my $major_allele_readcount = $items[$HETEROZYGOSITY_TABLE_MAJOR_ALLELE_READCOUNT_COLUMN];
		my $minor_allele_frequency = $items[$HETEROZYGOSITY_TABLE_MINOR_ALLELE_FREQUENCY_COLUMN];
		my $major_allele_frequency = $items[$HETEROZYGOSITY_TABLE_MAJOR_ALLELE_FREQUENCY_COLUMN];
		
		# verifies that input values all make sense
		if($position !~ /^\d+$/)
		{
			print STDERR "Warning: position is not non-zero integer"
				.$position." in heterozygosity table:\n\t"
				.$potential_contaminated_within_sample_diversity_file."\n";
		}
		if($minor_allele_readcount !~ /^\d+$/)
		{
			print STDERR "Warning: minor allele readcount is not non-zero integer"
				.$minor_allele_readcount." in heterozygosity table:\n\t"
				.$potential_contaminated_within_sample_diversity_file."\n";
		}
		if($major_allele_readcount !~ /^\d+$/)
		{
			print STDERR "Warning: major allele readcount is not non-zero integer"
				.$major_allele_readcount." in heterozygosity table:\n\t"
				.$potential_contaminated_within_sample_diversity_file."\n";
		}
		if($minor_allele_frequency !~ /^[\de.-]+$/)
		{
			print STDERR "Warning: non-numerical minor allele frequency "
				.$minor_allele_frequency." in heterozygosity table:\n\t"
				.$potential_contaminated_within_sample_diversity_file."\n";
		}
		if($major_allele_frequency !~ /^[\de.-]+$/)
		{
			print STDERR "Warning: non-numerical major allele frequency "
				.$major_allele_frequency." in heterozygosity table:\n\t"
				.$potential_contaminated_within_sample_diversity_file."\n";
		}
		if($minor_allele !~ /^[ATCG]$/)
		{
			print STDERR "Warning: minor allele is not A, T, C, or G in heterozygosity table:\n\t"
				.$potential_contaminated_within_sample_diversity_file."\n";
		}
		if($major_allele !~ /^[ATCG]$/)
		{
			print STDERR "Warning: major allele is not A, T, C, or G in heterozygosity table:\n\t"
				.$potential_contaminated_within_sample_diversity_file."\n";
		}
	
		# only includes positions with minor allele readcount >= 10, minor allele frequency >= 3%
		# assumes that major allele frequency = 100% - minor allele frequency
		if($minor_allele_readcount >= $minimum_minor_allele_readcount
			and $minor_allele_frequency >= $minimum_minor_allele_frequency)
		{
			if($positions_with_heterozygosity{$position})
			{
				print STDERR "Warning: position appears in more than one line in "
					."heterozygosity table:\n\t".$potential_contaminated_within_sample_diversity_file."\n";
			}
			else
			{
				$positions_with_heterozygosity{$position} = 1;
				$number_positions_with_heterozygosity++;
			}
			
			$minor_alleles{$position}{$minor_allele} = 1;
			$major_alleles{$position}{$major_allele} = 1;
			
			$minor_allele_frequencies{$position}{$minor_allele} = $minor_allele_frequency;
			$major_allele_frequencies{$position}{$major_allele} = $major_allele_frequency;
		
			if($list_of_alleles_to_print)
			{
				$list_of_alleles_to_print .= $ALLELE_LIST_SEPARATOR;
			}
			$list_of_alleles_to_print .= prepare_allele_to_print($position, $major_allele.$minor_allele);
		}
	}
	close HETEROZYGOSITY_TABLE;
	
	# verifies that this sample is worth looking at
# 	if($number_positions_with_heterozygosity < $MINIMUM_MINOR_ALLELES_FOUND_IN_CONSENSUS)
# 	{
# 		return;
# 	}

	# retrieves consensus-level alleles for potential contaminated sample and contaminating sample
	my $potential_contaminated_consensus = $sequence_name_to_consensus{$potential_contaminated_sample};
	my $potential_contaminating_consensus = $sequence_name_to_consensus{$potential_contaminating_sample};
	
	# breaks consensus genomes into alleles
	my @potential_contaminated_consensus_values = split(//, $potential_contaminated_consensus);
	my @potential_contaminating_consensus_values = split(//, $potential_contaminating_consensus);
	
	# counts unambiguous bases in contaminated and contaminating consensus genomes
	my $potential_contaminated_consensus_unambig_bases = count_unambiguous_bases_in_sequence(@potential_contaminated_consensus_values);
	my $potential_contaminating_consensus_unambig_bases = count_unambiguous_bases_in_sequence(@potential_contaminating_consensus_values);
	my $potential_contaminated_consensus_percent_covered = $potential_contaminated_consensus_unambig_bases / $reference_sequence_length;
	my $potential_contaminating_consensus_percent_covered = $potential_contaminating_consensus_unambig_bases / $reference_sequence_length;
	
	# exits if either sample does not have minimum genome coverage
	if($potential_contaminated_consensus_percent_covered < $minimum_genome_coverage
		or $potential_contaminating_consensus_percent_covered < $minimum_genome_coverage)
	{
		return;
	}
	
	# counts matches and mismatches at positions with heterozygosity
	my $number_mismatches = 0; # positions where contaminating sample allele does not match minor, major, or consensus sequence allele in contaminated sequence
	my $mismatches_string = "";
	
	my $minor_alleles_matched = 0; # positions where contaminating sample allele matches contaminated sample minor allele
	my $major_alleles_matched = 0; # positions where contaminating sample allele matches contaminated sample major allele
	my $list_of_alleles_matched = ""; # alleles matched, at positions with heterozygosity only, for printing
	my @matched_allele_frequencies = ();
	
	foreach my $position(sort {$a <=> $b} keys %positions_with_heterozygosity)
	{
		my $nucleotide_at_position = $potential_contaminating_consensus_values[$position - 1];
		if(is_unambiguous_base($nucleotide_at_position))
		{
			if($minor_alleles{$position}{$nucleotide_at_position}
				or $major_alleles{$position}{$nucleotide_at_position})
			{
				# saves matched allele to string to print
				if($list_of_alleles_matched)
				{
					$list_of_alleles_matched .= $ALLELE_LIST_SEPARATOR;
				}
				$list_of_alleles_matched .= prepare_allele_to_print($position, $nucleotide_at_position);
				
				# increments count of minor or major alleles matched
				if($minor_alleles{$position}{$nucleotide_at_position})
				{
					$minor_alleles_matched++;
					push(@matched_allele_frequencies, $minor_allele_frequencies{$position}{$nucleotide_at_position});
				}
				elsif($major_alleles{$position}{$nucleotide_at_position})
				{
					$major_alleles_matched++;
					push(@matched_allele_frequencies, $major_allele_frequencies{$position}{$nucleotide_at_position});
				}
			}
			else
			{
				# this base does not match minor or major allele in contaminated sequence
				$number_mismatches++;
				
				# verifies that this potential contamination scenario is still worth looking at
				if($number_mismatches > $maximum_allowed_mismatches)
				{
					return;
				}
				
				# saves mismatch to string to print
				if($mismatches_string)
				{
					$mismatches_string .= $ALLELE_LIST_SEPARATOR;
				}
				$mismatches_string .= prepare_allele_to_print($position, $nucleotide_at_position);
			}
		}
	}

	# verifies that this potential contamination scenario is still worth looking at
# 	if($minor_alleles_matched < $MINIMUM_MINOR_ALLELES_FOUND_IN_CONSENSUS)
# 	{
# 		return;
# 	}
	
	# counts mismatches in rest of genome
	foreach my $position(1..$#potential_contaminating_consensus_values + 1)
	{
		my $potential_contaminated_consensus_value = $potential_contaminated_consensus_values[$position - 1];
		my $potential_contaminating_consensus_value = $potential_contaminating_consensus_values[$position - 1];
		
		if( # not a location with heterozygosity
			!$positions_with_heterozygosity{$position}
		
			 # mismatch between contaminated and contaminating consensus sequences
			and $potential_contaminated_consensus_value ne $potential_contaminating_consensus_value
			
			# bases in both contaminating and contaminated sequence are A, T, C, or G
			and is_unambiguous_base($potential_contaminating_consensus_value)
			and is_unambiguous_base($potential_contaminated_consensus_value))
		{
			# this base does not match consensus genome
			$number_mismatches++;
			
			# verifies that this potential contamination scenario is still worth looking at
			if($number_mismatches > $maximum_allowed_mismatches)
			{
				return;
			}
			
			# saves mismatch to string to print
			if($mismatches_string)
			{
				$mismatches_string .= $ALLELE_LIST_SEPARATOR;
			}
			$mismatches_string .= prepare_allele_to_print($position, $potential_contaminating_consensus_value);
		}
	}
	
	# summarizes matched allele frequencies
	my $median_frequency = 1; # default 100% if only consensus sequences match (no heterozygous positions)
	my $frequency_range_min = $NO_DATA;
	my $frequency_range_max = $NO_DATA;
	my $matched_allele_frequency_string = "";
	
	if(scalar @matched_allele_frequencies)
	{
		# calculates median
		@matched_allele_frequencies = sort @matched_allele_frequencies;
		if(scalar @matched_allele_frequencies == 1)
		{
			$median_frequency = $matched_allele_frequencies[0];
		}
		elsif(scalar @matched_allele_frequencies % 2 == 0) # even
		{
			# average of the two center values
			$median_frequency = ($matched_allele_frequencies[(scalar @matched_allele_frequencies)/2 + 1]
				+ $matched_allele_frequencies[(scalar @matched_allele_frequencies)/2]) / 2;
		}
		else # odd
		{
			$median_frequency = $matched_allele_frequencies[(scalar @matched_allele_frequencies)/2 + 1]
		}
		
		# generates range
		$frequency_range_min = $matched_allele_frequencies[0];
		$frequency_range_max = $matched_allele_frequencies[(scalar @matched_allele_frequencies) - 1];
	
		# generates list
		foreach my $allele_frequency(@matched_allele_frequencies)
		{
			if($matched_allele_frequency_string)
			{
				$matched_allele_frequency_string .= $PERCENTAGES_LIST_SEPARATOR;
			}
			$matched_allele_frequency_string .= prepare_percentage_to_print($allele_frequency);
		}
	}
	
	# summarizes contamination type (minor or consensus-level)
	my $contamination_type = $NO_DATA;
	if($minor_alleles_matched and !$major_alleles_matched)
	{
		$contamination_type = "minor alleles";
	}
	elsif($minor_alleles_matched and $major_alleles_matched)
	{
		$contamination_type = "minor and consensus-level";
	}
	elsif(!$minor_alleles_matched)
	{
		$contamination_type = "consensus-level";
	}
	
	# if we've gotten this far, this sequence could be a contaminating sequence
	# prepare and return output to print
	my $output_line = "";
	my $output_line_plate = "";
	
	# adds columns about contaminated sample
	$output_line .= $potential_contaminated_sample.$DELIMITER;
	$output_line .= add_comma_separators($potential_contaminated_consensus_unambig_bases).$DELIMITER;
	$output_line .= prepare_percentage_to_print($potential_contaminated_consensus_percent_covered).$DELIMITER;
	$output_line .= $number_positions_with_heterozygosity.$DELIMITER;
	$output_line .= $list_of_alleles_to_print.$DELIMITER;
	
	# adds columns about contaminating sample
	$output_line .= $potential_contaminating_sample.$DELIMITER;
	$output_line .= add_comma_separators($potential_contaminating_consensus_unambig_bases).$DELIMITER;
	$output_line .= prepare_percentage_to_print($potential_contaminating_consensus_percent_covered).$DELIMITER;
	
	# adds columns about positions with and without matched alleles
	$output_line .= $minor_alleles_matched.$DELIMITER;
	$output_line .= $major_alleles_matched.$DELIMITER;
	
	if($number_positions_with_heterozygosity)
	{
		my $percent_heterozygous_positions_matched = (($minor_alleles_matched+$major_alleles_matched)/$number_positions_with_heterozygosity);
		$output_line .= prepare_percentage_to_print($percent_heterozygous_positions_matched).$DELIMITER;
	}
	else
	{
		$output_line .= $NO_DATA.$DELIMITER;
	}
	
	$output_line .= $list_of_alleles_matched.$DELIMITER;
	$output_line .= $number_mismatches.$DELIMITER;
	$output_line .= $mismatches_string.$DELIMITER;
	
	# adds columns for contaminating allele frequencies
	$output_line .= $contamination_type.$DELIMITER;
	$output_line .= prepare_percentage_to_print($median_frequency).$DELIMITER; # approximate contamination volume
	if($frequency_range_min eq $NO_DATA or $frequency_range_min == $frequency_range_max)
	{
		$output_line .= prepare_percentage_to_print($frequency_range_min).$DELIMITER;
	}
	else
	{
		$output_line .= prepare_percentage_to_print($frequency_range_min)
			." - ".prepare_percentage_to_print($frequency_range_max).$DELIMITER;
	}
	
	$output_line .= $matched_allele_frequency_string;
	
	# if plate map(s) provided, adds columns for plate map positions
	if(scalar @plate_map_files)
	{
		# for main output file
		$output_line .= $DELIMITER;
		$output_line .= $sample_name_to_all_plate_positions{$potential_contaminated_sample}.$DELIMITER;
		$output_line .= $sample_name_to_all_plate_positions{$potential_contaminating_sample};
		if(scalar @plate_map_files > 1)
		{
			$output_line .= $DELIMITER;
			$output_line .= $sample_name_to_all_plates{$potential_contaminated_sample};
			$output_line .= $DELIMITER;
			$output_line .= $sample_name_to_all_plates{$potential_contaminating_sample};
		}
		
		# for plate map-specific output file
		$output_line_plate .= $sample_name_to_all_plate_positions{$potential_contaminated_sample}.$DELIMITER;
		$output_line_plate .= $potential_contaminated_sample.$DELIMITER;
		$output_line_plate .= $sample_name_to_all_plate_positions{$potential_contaminating_sample}.$DELIMITER;
		$output_line_plate .= $potential_contaminating_sample.$DELIMITER;
		$output_line_plate .= $median_frequency;
	}

	return ($output_line, $output_line_plate);
}

# finishes processing within-sample diversity file for sample, starting at whatever stage
# we have input file for
# if we have bam file, runs LoFreq to get vcf file
# if we have vcf file, runs vcf_file_to_heterozygosity_table.pl to get heterozygosity table
# returns final within-sample diversity file (heterozygosity table)
sub process_within_sample_diversity_file_for_sample
{
	my $sample = $_[0];
	my $within_sample_diversity_file = $sample_name_to_within_sample_diversity_file{$sample};
	my $stage = $within_sample_diversity_file_to_stage{$within_sample_diversity_file};
	
	# if we have trimmed aligned bam file, generates vcf file
	if($stage eq "bam")
	{
		# runs LoFreq for bam -> vcf
		my $output_vcf_file = $temp_intermediate_directory.retrieve_file_name($within_sample_diversity_file."_LoFreq.vcf");
		check_if_file_exists_before_writing($output_vcf_file);
		print STDOUT "$LOFREQ_EXECUTABLE_FILE_PATH call -f $reference_genome_file -o $output_vcf_file $within_sample_diversity_file\n" if $verbose;
		`$LOFREQ_EXECUTABLE_FILE_PATH call -f $reference_genome_file -o $output_vcf_file $within_sample_diversity_file`;
		print STDOUT "\n" if $verbose;
		
		# updates within-sample diversity file saved for this sample
		$within_sample_diversity_file = $output_vcf_file;
		$stage = "vcf";
	}
	
	# if we have vcf file, generates heterozygosity table
	if($stage eq "vcf")
	{
		# runs vcf_file_to_heterozygosity_table.pl for vcf -> heterozygosity table
		my $output_heterozygosity_table = $temp_intermediate_directory.retrieve_file_name($within_sample_diversity_file)."_heterozygosity.txt";
		check_if_file_exists_before_writing($output_heterozygosity_table);
		print STDOUT "$VCF_TO_HETEROZYGOSITY_TABLE_SCRIPT_FILE_PATH $within_sample_diversity_file > $output_heterozygosity_table\n" if $verbose;
		`$VCF_TO_HETEROZYGOSITY_TABLE_SCRIPT_FILE_PATH $within_sample_diversity_file > $output_heterozygosity_table`;
		print STDOUT "\n" if $verbose;
		
		# updates within-sample diversity file saved for this sample
		$within_sample_diversity_file = $output_heterozygosity_table;
		$stage = "het";
	}
	
	# we should now have heterozygosity table
	if($stage ne "het")
	{
		print STDERR "Error: stage of within-sample diversity file not recognized: "
			.$within_sample_diversity_file_to_stage{$within_sample_diversity_file}
			.". Exiting.\n";
		die;
	}
	
	return $within_sample_diversity_file;
}

# counts number of bases in sequence that are A, T, C, or G
sub count_unambiguous_bases_in_sequence
{
	my @bases = @_; # must be all caps
	
	my $unambiguous_bases = 0;
	foreach my $base(@bases)
	{
		if(is_unambiguous_base($base))
		{
			$unambiguous_bases++;
		}
	}
	return $unambiguous_bases;
}

# returns 1 if base is A, T, C, G; returns 0 if not
# input base must be capitalized
sub is_unambiguous_base
{
	my $base = $_[0]; # must be capitalized
	if($base eq "A" or $base eq "T" or $base eq "C" or $base eq "G")
	{
		return 1;
	}
	return 0;
}


# HELPER FUNCTIONS FOR HANDLING PLATE NEIGHBORS

# retrieves samples neighboring given plate position
# assumes that %plate_position_to_sample_name and %sample_name_to_plate_position have been read in
# input: plate position (example input: H9)
# output: list of neighboring plate positions
# (example output: sample at H8, sample at H10, sample at G9, sample at I9)
sub retrieve_samples_neighboring_plate_position
{
	my $plate_position = $_[0];
	my $plate_position_sample_name = $_[1];
	
	my @neighbors = ();
	if($plate_position =~ /^([A-Z])(\s*)(\d+)$/)
	{
		my $letter = $1;
		my $number = $3;
		my $whitespace = $2;
		
		# retrieve neighboring rows and columns
		my $previous_letter = get_previous_letter($letter);
		my $next_letter = get_next_letter($letter);
		
		my $previous_number = $number - 1;
		my $next_number = $number + 1;
		
		# retrieve neighbors
		if($compare_whole_plate_map)
		{
			# retrieves all samples on plate
			foreach my $sample_name(keys %sample_name_to_plate_position)
			{
				push(@neighbors, $sample_name);
			}
		}
		else
		{
			my @neighbor_positions = ();
			if($compare_direct_neighbors)
			{
				# retrieves neighbor positions above, below, to the left, to the right
				push(@neighbor_positions, $previous_letter.$whitespace.$number); # top
				push(@neighbor_positions, $next_letter.$whitespace.$number); # bottom
				push(@neighbor_positions, $letter.$whitespace.$previous_number); # left
				push(@neighbor_positions, $letter.$whitespace.$next_number); # right
			}
			if($compare_diagonal_neighbors)
			{
				# retrieves neighbor positions at top right, top left, bottom right, bottom left
				push(@neighbor_positions, $previous_letter.$whitespace.$previous_number); # top left
				push(@neighbor_positions, $previous_letter.$whitespace.$next_number); # top right
				push(@neighbor_positions, $next_letter.$whitespace.$previous_number); # bottom left
				push(@neighbor_positions, $next_letter.$whitespace.$next_number); # bottom right
			}
			if($compare_direct_neighbors or $compare_diagonal_neighbors)
			{
				# retrieves neighbors belonging to neighbor positions
				foreach my $neighbor_position(@neighbor_positions)
				{
					if($plate_position_to_sample_name{$neighbor_position})
					{
						push(@neighbors, $plate_position_to_sample_name{$neighbor_position});
					}
				}
			}
			
			if($compare_row)
			{
				# retrieves all samples on same row (e.g., row A) on plate
				# compares this plate position's letter to letters of all plate positions
				foreach my $sample_name(keys %sample_name_to_plate_position)
				{
					my $plate_position_option = $sample_name_to_plate_position{$sample_name};
					if($plate_position_option =~ /^([A-Z])\s*\d+$/)
					{
						my $plate_position_option_letter = $1;
						if($letter eq $plate_position_option_letter
							and $plate_position_option ne $plate_position)
						{
							push(@neighbors, $sample_name);
						}
					}
				}
			}
			if($compare_column)
			{
				# retrieves all samples on same column (e.g., column 8) on plate
				# compares this plate position's number to number of all plate positions
				foreach my $sample_name(keys %sample_name_to_plate_position)
				{
					my $plate_position_option = $sample_name_to_plate_position{$sample_name};
					if($plate_position_option =~ /^([A-Z])\s*\d+$/)
					{
						my $plate_position_option_number = $1;
						if($number eq $plate_position_option_number
							and $plate_position_option ne $plate_position)
						{
							push(@neighbors, $sample_name);
						}
					}
				}
			}
		}
	}
	
	# remove duplicate neighbors, empty neighbors,
	# and neighbors identical to sample we're getting neighbors of
	my %neighbors_hash = ();
	foreach my $neighbor(@neighbors)
	{
		if($neighbor and $neighbor ne $plate_position_sample_name)
		{
			$neighbors_hash{$neighbor} = 1;
		}
	}
	@neighbors = keys %neighbors_hash;
	
	return @neighbors;
}

# returns 1 if input is valid plate position, 0 if not
# plate position is character followed by number (example A8)
sub is_valid_plate_position
{
	my $plate_position = $_[0];
	
	if($plate_position =~ /^[A-Z]\s*\d+$/)
	{
		return 1;
	}
	return 0;
}

# returns next letter
# returns empty string if Z
sub get_next_letter
{
	my $letter = $_[0];
	return $LETTER_TO_NEXT_LETTER{$letter};
}

# returns previous letter
# returns empty string if A
sub get_previous_letter
{
	my $letter = $_[0];
	return $LETTER_TO_PREVIOUS_LETTER{$letter};
}


# HELPER FUNCTIONS FOR HANDLING FILE PATHS

# retrieves file name from file path
# input: file path, ex. /Users/lakras/filepath.txt
# output: file name, ex. filepath.txt
sub retrieve_file_name
{
	my $file_path = $_[0];
	if($file_path =~ /.*\/([^\/]+)$/)
	{
		return $1;
	}
	return $file_path;
}

# retrieves directory of file file path
# input: file path, ex. /Users/lakras/filepath.txt
# output: file directory, ex. /Users/lakras/
sub retrieve_file_directory
{
	my $file_path = $_[0];
	if($file_path =~ /^(.*\/)[^\/]+$/)
	{
		return $1;
	}
	return "/";
}

# trims input file path from the right through first .
# returns empty string if no .
# example input: /Users/lakras/sample1.ext1.ext2.ext3
# example output: /Users/lakras/sample1.ext1.ext2
sub trim_off_file_extension
{
	my $file_name = $_[0];
	if($file_name =~ /^(.*)[.][^.]+$/)
	{
		return $1;
	}
	return "";
}

# checks if file exists, to be run before writing to a file
# if file already exists at provided path and $overwrite is set to true (1), prints a warning
# if file already exists at provided path and $overwrite is set to true (1), prints an error and exits
sub check_if_file_exists_before_writing
{
	my $file_path = $_[0];

	if(-e $file_path)
	{
		if($overwrite)
		{
			print STDERR "Warning";
		}
		else
		{
			print STDERR "Error";
		}
		print STDERR ": file already exists at file path to write to:\n\t".$file_path."\n";
		if(!$overwrite)
		{
			print STDERR "Exiting.\n";
			die;
		}
	}
}

# retrieves current working directory, ending in /
sub retrieve_current_working_directory
{
	my $current_working_directory = `pwd`;
	chomp $current_working_directory;
	
	if($current_working_directory !~ /\/$/) # if doesn't end in /
	{
		$current_working_directory .= "/"; # adds / to the end
	}
	return $current_working_directory;
}


# HELPER FUNCTIONS FOR PRINTING

# input: floating point value (example: 3.14159)
# output: value rounded to first decimal point (example: 3.1)
sub round_value
{
	my $value = $_[0];
	return sprintf("%.1f", $value);
}

sub prepare_percentage_to_print
{
	my $percent = $_[0];
	if($percent eq $NO_DATA)
	{
		return $NO_DATA;
	}
	if($ROUND_PERCENTAGES)
	{
		$percent = round_value(100*$percent)."%";
	}
	return $percent;
}

sub prepare_allele_to_print
{
	my $position = $_[0];
	my $bases = $_[1];
	
	return add_comma_separators($position)." ".$bases;
}

# adds comma thousands separator(s)
# from https://stackoverflow.com/questions/33442240/perl-printf-to-use-commas-as-thousands-separator
sub add_comma_separators
{
	my $value = $_[0];
	
	my $text = reverse $value;
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}




# HELPER FUNCTIONS FOR HANDING INPUT ARGUMENTS

# checks if parameter argument was supplied and, if so, reads in input file path that follows
# returns input file supplied after argument
# returns -1 if this argument was not entered
sub read_in_input_file_argument
{
	my $argument_option_1 = $_[0]; # for example: -f
	my $argument_option_2 = $_[1]; # for example: --ref
	
	if($argument eq $argument_option_1 or $argument eq $argument_option_2) # next item is input file
	{
		my $next_item = $ARGV[$argument_index+1];
		if($argument_index + 1 > $#ARGV or $next_item =~ /^-\w$/ or $next_item =~ /^--[\w-]+$/)
		{
			print STDERR "Warning: ignoring ".$argument_option_1." | ".$argument_option_2." argument with no input file.\n";
		}
		else
		{
			$argument_index++;
			return $next_item; # return input file
		}
	}
	return -1; # this argument was not entered
}

# checks if parameter argument was supplied and, if so, reads in positive integer that follows
# returns positive integer supplied after argument
# returns -1 if this argument was not entered
sub read_in_positive_integer_argument
{
	my $argument_option_1 = $_[0]; # for example: -p
	my $argument_option_2 = $_[1]; # for example: --cores
	
	if($argument eq $argument_option_1 or $argument eq $argument_option_2)
	{
		my $next_item = $ARGV[$argument_index+1];
		if($argument_index + 1 > $#ARGV or $next_item !~ /^\d+$/)
		{
			print STDERR "Warning: ignoring ".$argument_option_1." | ".$argument_option_2." argument with no int.\n";
		}
		else
		{
			$argument_index++;
			return $next_item;
		}
	}
	return -1; # this argument was not entered
}

# checks if parameter argument was supplied and, if so, reads in positive float that follows
# returns positive float supplied after argument
# returns -1 if this argument was not entered
sub read_in_positive_float_argument
{
	my $argument_option_1 = $_[0]; # for example: -p
	my $argument_option_2 = $_[1]; # for example: --cores
	
	if($argument eq $argument_option_1 or $argument eq $argument_option_2)
	{
		my $next_item = $ARGV[$argument_index+1];
		if($argument_index + 1 > $#ARGV or $next_item !~ /^[\d.]+$/)
		{
			print STDERR "Warning: ignoring ".$argument_option_1." | ".$argument_option_2." argument with no float.\n";
		}
		else
		{
			$argument_index++;
			return $next_item;
		}
	}
	return -1; # this argument was not entered
}

# checks if parameter argument was supplied and, if so, reads in boolean that follows
# accepts TRUE, FALSE, True, False, true, false, t, f, T, F, 1, 0
# returns boolean supplied after argument (as 0 or 1)
# returns -1 if this argument was not entered
sub read_in_boolean_argument
{
	my $argument_option_1 = $_[0]; # for example: -d
	my $argument_option_2 = $_[1]; # for example: --diagonal
	
	if($argument eq $argument_option_1 or $argument eq $argument_option_2)
	{
		my $next_item = $ARGV[$argument_index+1];
		if($argument_index + 1 > $#ARGV or $next_item =~ /^-\w$/ or $next_item =~ /^--[\w-]+$/)
		{
			print STDERR "Warning: ignoring ".$argument_option_1." | ".$argument_option_2." argument with no boolean.\n";
		}
		else
		{
			$argument_index++;
			$next_item = uc($next_item);
			if($next_item eq "TRUE" or $next_item eq "T" or $next_item eq "1")
			{
				return 1;
			}
			if($next_item eq "FALSE" or $next_item eq "F" or $next_item eq "0")
			{
				return 0;
			}
		}
	}
	return -1; # this argument was not entered
}

# checks if parameter argument was supplied and, if so, reads in input file paths that follow
# returns reference to list of input files supplied after argument
# returns -1 if this argument was not entered
sub read_in_input_files_argument
{
	my $argument_option_1 = $_[0]; # for example: -c
	my $argument_option_2 = $_[1]; # for example: --consensus
	
	if($argument eq $argument_option_1 or $argument eq $argument_option_2)
	{
		my $next_item = $ARGV[$argument_index+1];
		my $items_added = 0;
		my @input_files = ();
		while($argument_index + 1 <= $#ARGV and $next_item !~ /^-\w$/ and $next_item !~ /^--[\w-]+$/)
		{
			push(@input_files, $next_item);
			$items_added++;
		
			$argument_index++;
			$next_item = $ARGV[$argument_index+1];
		}
		if(!$items_added)
		{
			print STDERR "Warning: ignoring ".$argument_option_1." | ".$argument_option_2." argument with no input file.\n";
		}
		return \@input_files;
	}
	return -1; # this argument was not entered
}

# returns "TRUE" if 1 or >1, "FALSE" if 0 or less
sub int_to_bool_string
{
	my $int_representation_of_bool = $_[0];
	if($int_representation_of_bool > 0)
	{
		return "TRUE";
	}
	return "FALSE";
}

# June 1, 2021
