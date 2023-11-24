#!/usr/bin/env perl

# Detects potential cross-contamination by comparing consensus-level and minor alleles.

use strict;
use warnings;

# tools used:
use Math::Round;
use Parallel::ForkManager; # download here: https://metacpan.org/pod/Parallel::ForkManager
my $LOFREQ_EXECUTABLE_FILE_PATH = "lofreq";
my $MAFFT_EXECUTABLE_FILE_PATH = "mafft";
my $SAMTOOLS_EXECUTABLE_FILE_PATH = "samtools";
my $VCF_TO_HETEROZYGOSITY_TABLE_SCRIPT_FILE_PATH = "vcf_file_to_heterozygosity_table.pl";
my $PLATE_VISUALIZATION_FILE_PATH = "visualize_plate_map.R";

# plate map input file:
my $PLATE_MAP_SAMPLE_COLUMN = 0;
my $PLATE_MAP_POSITION_COLUMN = 1;

# standard plate layouts:
#   6-well plate    3 columns (A, B, C) x 2 rows (1, 2)
#   12-well plate   4 columns x 3 rows
#   24-well plate   6 x 4
#   48-well plate   8 x 6
#   96-well plate   12 x 8
#   384-well plate  24 x 16
#   1536-well plate 48 x 32
#   3456-well plate 72 x 48
my %STANDARD_PLATE_SIZE_TO_NUMBER_COLUMNS = (
	6 => 3, 12 => 4, 24 => 6, 48 => 8, 96 => 12, 384 => 24, 1536 => 48, 3456 => 72);

# for reading plate map:
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

# columns in read-depth tables produced by samtools:
my $READ_DEPTH_POSITION_COLUMN = 1;
my $READ_DEPTH_COLUMN = 2;

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
my $DEFAULT_MINIMUM_GENOME_COVERAGE = 0.95;
my $DEFAULT_MINIMUM_READ_DEPTH = 100;
my $DEFAULT_MINIMUM_MINOR_ALLELE_READCOUNT = 10; # ignore minor alleles with readcount <10 (does not consider the position to have heterozygosity)
my $DEFAULT_MINIMUM_MINOR_ALLELE_FREQUENCY = 0.03; # 3%

my $DEFAULT_OUTPUT_FILE_NAME = "potential_cross-contamination.txt";
my $DEFAULT_OVERWRITE = 0; # false
my $DEFAULT_CORES_TO_USE = 1;
my $DEFAULT_VERBOSE = 1; # true
my $DEFAULT_PRINT_ALL_ISNVS = 0; # false
my $DEFAULT_PRINT_ALL = 0; # false

my $DEFAULT_PLATE_SIZE = 96; # 12 columns, 8 rows
my $DEFAULT_PLATE_NUMBER_ROWS = 8; # A through H
my $DEFAULT_PLATE_NUMBER_COLUMNS = 12; # 1 through 12

my $DEFAULT_COMPARE_DIRECT_NEIGHBORS = 1; # true
my $DEFAULT_COMPARE_DIAGONAL_NEIGHBORS = 0; # false
my $DEFAULT_COMPARE_ROW = 0; # false
my $DEFAULT_COMPARE_COLUMN = 0; # false
my $DEFAULT_COMPARE_WHOLE_PLATE_MAP = 0; # false


# generates default directory and output file from current working directory
my $default_temp_intermediate_files_directory = retrieve_current_working_directory(); # current working directory
my $default_visualizations_directory = retrieve_current_working_directory(); # current working directory
my $default_output_file = $default_temp_intermediate_files_directory.$DEFAULT_OUTPUT_FILE_NAME;


# if no command line arguments supplied, prints options
if(!scalar @ARGV) # no command line arguments supplied
{
	print STDERR "\nDetects potential cross-contamination.\n";
	print STDERR "Usage: polyphonia detect_cross_contam [options]\n";
	print STDERR "\n";
	
	print STDERR "OPTIONS:\n";
	print STDERR "- Reference (required):\n";
	print STDERR "\t-f | --ref FILE\t\t\tReference fasta file [null]\n";
	print STDERR "\n";
	
	print STDERR "- Consensus genomes (aligned or not aligned, not both; at least one file required):\n";
	print STDERR "\t-c | --consensus FILE(S)\tUnaligned consensus genome or genomes [null]\n";
	print STDERR "\t-a | --consensus-aligned FILE\tFasta alignment with consensus genomes pre-aligned to reference; reference provided by --ref must appear first [null]\n";
	print STDERR "\n";
	
	print STDERR "- Within-sample diversity (any combination; at least one file required):\n";
	print STDERR "\t-b | --bam FILE(S)\t\tAligned and trimmed reads as bam file(s); must use reference provided by --ref [null]\n";
	print STDERR "\t-v | --vcf FILE(S)\t\tVCF file(s) output by LoFreq; must use reference provided by --ref [null]\n";
	print STDERR "\t-h | --het FILE(S)\t\tTab-separated heterozygosity summary tables; see documentation for format [null]\n";
	print STDERR "\n";
	
	print STDERR "- Filtering options:\n";
	print STDERR "\t-i | --min-maf FLOAT\t\tMinimum minor allele frequency for a position to be considered heterozygous [".$DEFAULT_MINIMUM_MINOR_ALLELE_FREQUENCY."]\n";
	print STDERR "\t-e | --min-readcount INT\tMinimum minor allele readcount for a position to be considered heterozygous [".$DEFAULT_MINIMUM_MINOR_ALLELE_READCOUNT."]\n";
	print STDERR "\t-r | --min-depth INT\t\tMinimum read depth for a position to be used for comparison [".$DEFAULT_MINIMUM_READ_DEPTH."]\n";
	print STDERR "\t-1 | --read-depths FILE(S)\tRead depth tables; provide alongside vcf files or heterozygosity tables if min-depth>0; see documentation for format [null]\n";
	print STDERR "\t-g | --min-covered FLOAT\tMinimum proportion genome that must be covered at minimum read depth for a sample to be included [".$DEFAULT_MINIMUM_GENOME_COVERAGE."]\n";
	print STDERR "\t-y | --max-mismatches INT\tIn flagged potential cross-contamination, maximum allowed unambiguous bases in contaminating sample consensus not matching contaminated sample alleles [".$DEFAULT_MAXIMUM_ALLOWED_MISMATCHES."]\n";
	print STDERR "\t-3 | --masked-positions STRING\t1-indexed positions to mask (e.g., 1-10,50,55-70) [null]\n";
	print STDERR "\t-4 | --masked-positions-file FILE\t1-indexed positions to mask, one per line [null]\n";
	print STDERR "\n";
	
	print STDERR "- Plate map and neighbors (any combination, all optional):\n";
	print STDERR "\t-m | --plate-map FILE(S)\tOptional plate map(s) (tab-separated, no header: sample name, plate position (e.g., A8)); provides substantial speed-up [null]\n";
	print STDERR "\t-z | --plate-size INT\t\tStandard plate size (6-well, 12-well, 24, 48, 96, 384, 1536, or 3456) [".$DEFAULT_PLATE_SIZE."]\n";
	print STDERR "\t-q | --plate-columns INT\tNumber columns in plate (e.g., 1, 2, 3, 4) [".$DEFAULT_PLATE_NUMBER_COLUMNS."]\n";
	print STDERR "\t-k | --plate-rows INT\t\tNumber rows in plate (e.g., A, B, C, D) [".$DEFAULT_PLATE_NUMBER_ROWS."]\n";
	print STDERR "\t-n | --compare-direct BOOL\tCompare direct plate neighbors (left, right, top, bottom) [".int_to_bool_string($DEFAULT_COMPARE_DIRECT_NEIGHBORS)."]\n";
	print STDERR "\t-d | --compare-diagonal BOOL\tCompare diagonal plate neighbors (top-right, bottom-right, top-left, bottom-left) [".int_to_bool_string($DEFAULT_COMPARE_DIAGONAL_NEIGHBORS)."]\n";
	print STDERR "\t-w | --compare-row BOOL\t\tCompare samples in the same row (e.g., row A) [".int_to_bool_string($DEFAULT_COMPARE_ROW)."]\n";
	print STDERR "\t-l | --compare-column BOOL\tCompare samples in the same column (e.g., column 8) [".int_to_bool_string($DEFAULT_COMPARE_COLUMN)."]\n";
	print STDERR "\t-t | --compare-plate BOOL\tCompare all samples in each plate map [".int_to_bool_string($DEFAULT_COMPARE_WHOLE_PLATE_MAP)."]\n";
	print STDERR "\n";
	
	print STDERR "- Output:\n";
	print STDERR "\t-o | --output FILE\t\tOutput file path [".$default_output_file."]\n";
	print STDERR "\t-s | --out-figures DIRECTORY\tPath of directory to store plate visualization files [".$default_visualizations_directory."]\n";
	print STDERR "\t-x | --out-temp DIRECTORY\tPath of directory to store intermediate and temporary files [".$default_temp_intermediate_files_directory."]\n";
	print STDERR "\n";
	
	print STDERR "- Misc:\n";
	print STDERR "\t-p | --cores INT\t\tOptional number of cores to use for preprocessing in parallel [".$DEFAULT_CORES_TO_USE."]\n";
	print STDERR "\t-u | --verbose BOOL\t\tPrint progress updates to STDOUT [".int_to_bool_string($DEFAULT_VERBOSE)."]\n";
	print STDERR "\t-j | --overwrite BOOL\t\tOverwrite files that already exist at output, intermediate, and temp file paths [".int_to_bool_string($DEFAULT_OVERWRITE)."]\n";
	print STDERR "\t-2 | --print-all-iSNVs BOOL\tInclude all threshold-passing samples in iSNVs visualizations, including samples without plate neighbors [".int_to_bool_string($DEFAULT_PRINT_ALL_ISNVS)."]\n";
	print STDERR "\t-0 | --print-all BOOL\t\tOutput outcomes of all comparisons (all comparisons are marked as potential cross-contamination) [".int_to_bool_string($DEFAULT_PRINT_ALL)."]\n";
	print STDERR "\n\n";
	exit;
}


# command line arguments set to default values
my $reference_genome_file = "";
my @consensus_genome_files = ();
my $consensus_genomes_aligned_file = "";
my $minimum_genome_coverage = $DEFAULT_MINIMUM_GENOME_COVERAGE;
my $maximum_allowed_mismatches = $DEFAULT_MAXIMUM_ALLOWED_MISMATCHES;
my $minimum_read_depth = $DEFAULT_MINIMUM_READ_DEPTH;
my $masked_positions_string = "";
my $masked_positions_file = "";

my @aligned_and_trimmed_bam_files = ();
my @vcf_files = ();
my @heterozygosity_tables = ();
my @read_depth_tables = ();

my $output_file_path = $default_output_file;
my $visualizations_directory = $default_visualizations_directory;
my $temp_intermediate_directory = $default_temp_intermediate_files_directory;

my $cores_to_use = $DEFAULT_CORES_TO_USE;
my $overwrite = $DEFAULT_OVERWRITE;
my $verbose = $DEFAULT_VERBOSE;
my $print_all_isnvs = $DEFAULT_PRINT_ALL_ISNVS;
my $print_all = $DEFAULT_PRINT_ALL;

my $minimum_minor_allele_readcount = $DEFAULT_MINIMUM_MINOR_ALLELE_READCOUNT;
my $minimum_minor_allele_frequency = $DEFAULT_MINIMUM_MINOR_ALLELE_FREQUENCY;

my @plate_map_files = ();
my $plate_size = $DEFAULT_PLATE_SIZE;
my $plate_size_entered = 0;
my $plate_number_rows = $DEFAULT_PLATE_NUMBER_ROWS;
my $plate_number_columns = $DEFAULT_PLATE_NUMBER_COLUMNS;
my $plate_row_count_entered = 0;
my $plate_column_count_entered = 0;
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
	elsif(($input = read_in_input_file_argument("-s", "--out-figures")) ne "-1")
	{
		$visualizations_directory = $input;
	}
	elsif(($input = read_in_input_file_argument("-x", "--out-temp")) ne "-1")
	{
		$temp_intermediate_directory = $input;
	}
	elsif(($input = read_in_boolean_argument("-j", "--overwrite")) != -1)
	{
		$overwrite = $input;
	}
	elsif(($input = read_in_positive_integer_argument("-y", "--max-mismatches")) != -1)
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
	elsif(($input = read_in_input_files_argument("-1", "--read-depths")) ne "-1")
	{
		push(@read_depth_tables, @$input);
	}
	elsif(($input = read_in_positive_integer_argument("-r", "--min-depth")) != -1)
	{
		$minimum_read_depth = $input;
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
	elsif(($input = read_in_input_files_argument("-m", "--plate-map")) ne "-1")
	{
		push(@plate_map_files, @$input);
	}
	elsif(($input = read_in_string_argument("-3", "--masked-positions")) ne "-1")
	{
		$masked_positions_string = $input;
	}
	elsif(($input = read_in_input_file_argument("-4", "--masked-positions-file")) ne "-1")
	{
		$masked_positions_file = $input;
	}
	elsif(($input = read_in_positive_integer_argument("-z", "--plate-size")) != -1)
	{
		$plate_size = $input;
		$plate_size_entered = 1;
	}
	elsif(($input = read_in_positive_integer_argument("-k", "--plate-rows")) != -1)
	{
		$plate_number_rows = $input;
		$plate_row_count_entered = 1;
	}
	elsif(($input = read_in_positive_integer_argument("-q", "--plate-columns")) != -1)
	{
		$plate_number_columns = $input;
		$plate_column_count_entered = 1;
	}
	elsif(($input = read_in_boolean_argument("-n", "--compare-direct")) != -1)
	{
		$compare_direct_neighbors = $input;
	}
	elsif(($input = read_in_boolean_argument("-d", "--compare-diagonal")) != -1)
	{
		$compare_diagonal_neighbors = $input;
	}
	elsif(($input = read_in_boolean_argument("-w", "--compare-row")) != -1)
	{
		$compare_row = $input;
	}
	elsif(($input = read_in_boolean_argument("-l", "--compare-column")) != -1)
	{
		$compare_column = $input;
	}
	elsif(($input = read_in_boolean_argument("-t", "--compare-plate")) != -1)
	{
		$compare_whole_plate_map = $input;
	}
	elsif(($input = read_in_boolean_argument("-u", "--verbose")) != -1)
	{
		$verbose = $input;
	}
	elsif(($input = read_in_boolean_argument("-2", "--print-all-isnvs")) != -1)
	{
		$print_all_isnvs = $input;
	}
	elsif(($input = read_in_boolean_argument("-0", "--print-all")) != -1)
	{
		$print_all = $input;
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
if(!scalar @consensus_genome_files and !$consensus_genomes_aligned_file)
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
	print STDERR "Error: no within-sample diversity files provided. Exiting.\n";
	die;
}
if(!scalar @plate_map_files and !$compare_direct_neighbors and !$compare_diagonal_neighbors
	and !$compare_row and !$compare_column and !$compare_whole_plate_map)
{
	print STDERR "Warning: plate map(s) supplied but plate map options all set to False. "
		."Plate map(s) will not be used.\n";
	@plate_map_files = ();
}
if($plate_number_rows < 1)
{
	print STDERR "Error: plate map rows ".$plate_number_rows." < 1.\n";
}
if($plate_number_columns < 1)
{
	print STDERR "Error: plate map columns ".$plate_number_columns." < 1.\n";
}
if($plate_size_entered and !@plate_map_files)
{
	print STDERR "Warning: plate size entered but no plate map. Ignoring plate size.\n";
	$plate_size_entered = 0;
}
if(!@plate_map_files and ($plate_column_count_entered or $plate_row_count_entered))
{
	print STDERR "Warning: plate map column or row count entered but no plate map. "
		."Ignoring plate map column and row count.\n";
	$plate_column_count_entered = 0;
	$plate_row_count_entered = 0;
}
if($plate_size_entered and $plate_column_count_entered)
{
	print STDERR "Warning: plate map column count entered as well as plate size "
		."entered. Substituting provided plate column count ".$plate_number_columns
		." in entered standard plate size ".$plate_size.".\n";
}
if($plate_size_entered and $plate_row_count_entered)
{
	print STDERR "Warning: plate map row count entered as well as plate size "
		."entered. Substituting provided plate row count ".$plate_number_rows
		." in entered standard plate size ".$plate_size.".\n";
}
if($minimum_genome_coverage < 0 or $minimum_genome_coverage > 1)
{
	print STDERR "Error: minimum genome coverage is not between 0 and 1. Exiting.\n";
	die;
}
if($minimum_read_depth < 0)
{
	print STDERR "Warning: minimum read depth is less than 0; not applying read depth "
		."filter.\n";
}
if($minimum_read_depth and (scalar @heterozygosity_tables or scalar @vcf_files)
	and !scalar @read_depth_tables)
{
	# does not apply read depth filter
	$minimum_read_depth = 0;
	print STDERR "Warning: not all input files are vcf files but no read depth files "
		."provided; not applying read depth filter.\n";
}


# verifies that all input files exist and are not empty
verify_input_file_exists_and_is_nonempty($reference_genome_file, "reference genome file", 1, 1);
foreach my $consensus_genome_file(@consensus_genome_files)
{
	verify_input_file_exists_and_is_nonempty($consensus_genome_file, "consensus genome file", 1, 0);
}
if($consensus_genomes_aligned_file)
{
	verify_input_file_exists_and_is_nonempty($consensus_genomes_aligned_file, "aligned consensus genomes file", 1, 1);
}
foreach my $aligned_and_trimmed_bam_file(@aligned_and_trimmed_bam_files)
{
	verify_input_file_exists_and_is_nonempty($aligned_and_trimmed_bam_file, "aligned and trimmed bam file", 1, 0);
}
foreach my $vcf_file(@vcf_files)
{
	verify_input_file_exists_and_is_nonempty($vcf_file, "vcf file", 1, 0);
}
foreach my $heterozygosity_table(@heterozygosity_tables)
{
	verify_input_file_exists_and_is_nonempty($heterozygosity_table, "heterozygosity table", 1, 0);
}
foreach my $plate_map_file(@plate_map_files)
{
	verify_input_file_exists_and_is_nonempty($plate_map_file, "plate map file", 1, 0);
}
foreach my $read_depth_table(@read_depth_tables)
{
	verify_input_file_exists_and_is_nonempty($read_depth_table, "read depth table", 1, 0);
}
if($masked_positions_file)
{
	verify_input_file_exists_and_is_nonempty($masked_positions_file, "masked positions file", 1, 0);
}

# option to print all iSNV counts set to FALSE if no plate maps provided
if(!scalar @plate_map_files)
{
	$print_all_isnvs = 0;
}

# retrieves dimensions of standard plate map size entered
# substitutes in entered dimensions if row or column specifically entered
if($plate_size_entered)
{
	if($STANDARD_PLATE_SIZE_TO_NUMBER_COLUMNS{$plate_size})
	{
		my $standard_plate_number_columns = $STANDARD_PLATE_SIZE_TO_NUMBER_COLUMNS{$plate_size};
		my $standard_plate_number_rows = $plate_size/$standard_plate_number_columns;
		
		if(!$plate_column_count_entered)
		{
			$plate_number_columns = $standard_plate_number_columns;
		}
		if(!$plate_row_count_entered)
		{
			$plate_number_rows = $standard_plate_number_rows;
		}
	}
	else
	{
		print STDERR "Error: plate size ".$plate_size." not recognized as standard plate "
			."size. Exiting.\n";
		die;
	}
}


# verifies that all positions in plate map fit in entered plate map
if(scalar @plate_map_files)
{
	print STDERR "verifying entered plate map positions are valid...\n" if $verbose;
	
	# retrieves all valid letters for this plate map
	my %valid_row_letters = (); # key: letter(s) designating valid row on plate -> 1
	my $first_valid_row_letter = "A";
	my $last_valid_row_letter = $first_valid_row_letter;
	$valid_row_letters{$last_valid_row_letter} = 1;
	my $letter_count = 1;
	while($letter_count < $plate_number_rows)
	{
		$last_valid_row_letter = get_next_letter($last_valid_row_letter);
		$valid_row_letters{$last_valid_row_letter} = 1;
		$letter_count++;
	}
	
	# verifies that all plate map positions are within range
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
				my $plate_position = uc $items[$PLATE_MAP_POSITION_COLUMN];
				
				if($plate_position =~ /^([A-Z]+)\s*(\d+)$/)
				{
					my $letter = $1;
					my $number = $2;
					
					# verifies that number is within range
					if($number < 1 or $number > $plate_number_columns)
					{
						print STDERR "Error: plate position ".$plate_position
							." of sample ".$sample_name." in plate map ".$plate_map_file
							." contains column number ".$number." outside allowed range of "
							.$plate_number_columns." columns (allowed column numbers 1-"
							.$plate_number_columns."). Use --plate-size or --plate-columns "
							."to set number columns in plate. Exiting.\n";
						die;
					}
					
					# verifies that letter is within range
					if(!$valid_row_letters{$letter})
					{
						print STDERR "Error: plate position ".$plate_position
							." of sample ".$sample_name." in plate map ".$plate_map_file
							." contains row letter ".$letter." outside allowed range of "
							.$plate_number_rows." rows (allowed row letters "
							.$first_valid_row_letter."-".$last_valid_row_letter
							."). Use --plate-size or --plate-rows "
							."to set number rows in plate. Exiting.\n";
						die;
					}
				}
				else
				{
					print STDERR "Warning: plate map position ".$plate_position
						." of sample ".$sample_name." on plate map ".$plate_map_file
						." could not be parsed and will be ignored. Plate map positions "
						."must be letter(s) (row number) followed by digit(s) (column "
						."number), e.g., H8.\n";
				}
			}
		}
		close PLATE_MAP;
	}
}


# prepares directories to write to
$temp_intermediate_directory = prepare_directory($temp_intermediate_directory);
if(scalar @plate_map_files) # visualizations only generated if plate maps are provided
{
	$visualizations_directory = prepare_directory($visualizations_directory);
}


# prints input files and options entered
# reference
print STDERR "\n" if $verbose;
print STDERR "REFERENCE:\n\t".$reference_genome_file."\n" if $verbose;

# consensus genome files
print STDERR "CONSENSUS GENOMES:\n" if $verbose;
foreach my $consensus_genome_file(@consensus_genome_files)
{
	print STDERR "\t".$consensus_genome_file."\n" if $verbose;
}
if($consensus_genomes_aligned_file)
{
	print STDERR "\tpre-aligned: ".$consensus_genomes_aligned_file."\n" if $verbose;
}

# within-sample diversity files
print STDERR "WITHIN-SAMPLE DIVERSITY:\n" if $verbose;
foreach my $aligned_and_trimmed_bam_file(@aligned_and_trimmed_bam_files)
{
	print STDERR "\t".$aligned_and_trimmed_bam_file."\n" if $verbose;
}
foreach my $vcf_file(@vcf_files)
{
	print STDERR "\tpre-processed vcf file: ".$vcf_file."\n" if $verbose;
}
foreach my $heterozygosity_table(@heterozygosity_tables)
{
	print STDERR "\tfully pre-processed heterozygosity table: ".$heterozygosity_table."\n" if $verbose;
}

# read depth tables
if(scalar @read_depth_tables and $minimum_read_depth > 0)
{
	print STDERR "READ DEPTH TABLES:\n" if $verbose;
	foreach my $read_depth_table(@read_depth_tables)
	{
		print STDERR "\t".$read_depth_table."\n" if $verbose;
	}
}

# optional plate map file(s) and related options
if(scalar @plate_map_files)
{
	print STDERR "PLATE MAP(S):\n" if $verbose;
	foreach my $plate_map_file(@plate_map_files)
	{
		print STDERR "\t".$plate_map_file."\n" if $verbose;
	}
	print STDERR "\t".$plate_number_rows." rows in plate\n" if $verbose;
	print STDERR "\t".$plate_number_columns." columns in plate\n" if $verbose;
	
	print STDERR "PLATE MAP USE:\n" if $verbose;
	if($compare_whole_plate_map)
	{
		print STDERR "Comparing all samples in the same plate map.\n" if $verbose;
	}
	else
	{
		if($compare_direct_neighbors)
		{
			print STDERR "\tComparing direct plate neighbors (left, right, top, bottom).\n" if $verbose;
		}
		if($compare_diagonal_neighbors)
		{
			print STDERR "\tComparing diagonal plate neighbors (top-right, bottom-right, top-left, bottom-left).\n" if $verbose;
		}
		if($compare_row)
		{
			print STDERR "\tComparing samples in the same row (e.g., row A).\n" if $verbose;
		}
		if($compare_column)
		{
			print STDERR "\tComparing samples in the same column (e.g., column 8).\n" if $verbose;
		}
	}
}

# options
print STDERR "OPTIONS:\n" if $verbose;
print STDERR "\tminimum read depth: ".$minimum_read_depth."\n" if $verbose;
print STDERR "\tminimum genome coverage: ".($minimum_genome_coverage*100)."%\n" if $verbose;
print STDERR "\tminimum minor allele readcount: ".$minimum_minor_allele_readcount."\n" if $verbose;
print STDERR "\tminimum minor allele frequency: ".($minimum_minor_allele_frequency*100)."%\n" if $verbose;
print STDERR "\tmaximum allowed mismatches: ".$maximum_allowed_mismatches."\n" if $verbose;
if($masked_positions_string)
{
	print STDERR "\tmasked positions: ".$masked_positions_string."\n" if $verbose;
}
if($masked_positions_file)
{
	print STDERR "\tmasked positions described in file: ".$masked_positions_file."\n" if $verbose;
}

# output files
print STDERR "OUTPUT:\n" if $verbose;
print STDERR "\toutput file: ".$output_file_path."\n" if $verbose;
print STDERR "\tplate visualization files: ".$visualizations_directory."\n" if $verbose;
print STDERR "\tintermediate and temporary files: ".$temp_intermediate_directory."\n" if $verbose;

# output options
if($print_all)
{
	print STDERR "\tPrinting all comparisons.\n" if $verbose;
}
if($print_all_isnvs)
{
	print STDERR "\tPrinting iSNVs for all threshold-passing samples in plate visualization file, including samples without plate neighbors.\n" if $verbose;
}
else
{
	print STDERR "\tPrinting iSNVs only for threshold-passing samples with plate neighbors.\n" if $verbose;
}
print STDERR "\n" if $verbose;


# retrieves masked positions from masked positions string argument
my %position_is_masked = ();
if($masked_positions_string)
{
	# split string on ,
	my @position_items = split(",", $masked_positions_string);
	foreach my $position_item(@position_items)
	{
		# fills in ranges
		if($position_item =~ /(\d+)\-(\d+)/)
		{
			my $first_position = $1;
			my $last_position = $2;
			
			for(my $position = $first_position; $position <= $last_position; $position++)
			{
				$position_is_masked{$position} = 1;
			}
		}
		
		# adds single position
		else
		{
			$position_is_masked{$position_item} = 1;
		}
	}
}
if($masked_positions_file)
{
	open MASKED_POSITIONS, "<$masked_positions_file" || die "Could not open $masked_positions_file to read; terminating =(\n";
	while(<MASKED_POSITIONS>) # for each line in the file
	{
		chomp;
		my $position = $_;
		if($position)
		{
			$position_is_masked{$position} = 1;
		}
	}
	close MASKED_POSITIONS;
}


# retrieves sample names from plate maps if possible
my %sample_names = (); # key: sample name to include in comparisons -> value: 1
if(scalar @plate_map_files)
{
	print STDERR "retrieving sample names from plate map(s)...\n" if $verbose;
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
	print_number_samples_remaining_and_exit_if_none();
	
	# catalogues consensus genome sample names
	print STDERR "removing samples without associated consensus genome...\n" if $verbose;
	my %sample_has_consensus_genome = (); # key: sample name -> value: 1 if sample has associated consensus genome
	foreach my $consensus_genome_fasta_file(@consensus_genome_files, $consensus_genomes_aligned_file)
	{
		if($consensus_genome_fasta_file)
		{
			open FASTA_FILE, "<$consensus_genome_fasta_file" || die "Could not open $consensus_genome_fasta_file to read; terminating =(\n";
			while(<FASTA_FILE>) # for each line in the file
			{
				chomp;
				if($_ =~ /^>(.*)/) # header line
				{
					my $sample_name = $1;
					$sample_has_consensus_genome{$sample_name} = 1;
				}
			}
			close FASTA_FILE;
		}
	}
	
	# prints number of samples remaining
	print_number_samples_remaining_and_exit_if_none();
	
	# removes samples that don't have a consensus genome
	foreach my $sample_name(keys %sample_names)
	{
		if(!$sample_has_consensus_genome{$sample_name})
		{
			delete $sample_names{$sample_name};
		}
	}
}

# if no plate map, retrieves sample names from consensus genome fasta files
else
{
	print STDERR "retrieving sample names from consensus genome fasta file(s)...\n" if $verbose;
	foreach my $consensus_genome_fasta_file(@consensus_genome_files, $consensus_genomes_aligned_file)
	{
		if($consensus_genome_fasta_file)
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
	}
	
	# prints number of samples remaining
	print_number_samples_remaining_and_exit_if_none();
}

# records file stage for each within-sample diversity file
print STDERR "recording stage of each within-sample diversity file...\n" if $verbose;
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
print STDERR "retrieving within-sample diversity file for each sample...\n" if $verbose;
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

# removes samples that don't have a within-sample diversity file
print STDERR "removing samples without within-sample diversity file...\n" if $verbose;
foreach my $sample_name(keys %sample_names)
{
	if(!$sample_name_to_within_sample_diversity_file{$sample_name})
	{
		delete $sample_names{$sample_name};
	}
}

# prints number of samples remaining
print_number_samples_remaining_and_exit_if_none();

print STDERR "removing samples with non-existent within-sample diversity file...\n" if $verbose;
foreach my $sample_name(keys %sample_names)
{
	if(!-e $sample_name_to_within_sample_diversity_file{$sample_name})
	{
		delete $sample_names{$sample_name};
	}
}

# prints number of samples remaining
print_number_samples_remaining_and_exit_if_none();


# retrieves read depth file for each sample
my %sample_name_to_read_depth_file = (); # key: sample name -> value: file path of read depth file
if($minimum_read_depth and scalar @read_depth_tables)
{
	print STDERR "retrieving read depth table for each sample...\n" if $verbose;
	foreach my $file_path(@read_depth_tables)
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
				$sample_name_to_read_depth_file{$potential_sample_name} = $file_path;
				$sample_name_found = 1;
			}
			else
			{
				$potential_sample_name = trim_off_file_extension($potential_sample_name);
			}
		}
	}
}

# removes samples that don't have a read depth file or bam file
if($minimum_read_depth)
{
	print STDERR "removing samples without read depth table or bam file to generate it from...\n" if $verbose;
	foreach my $sample_name(keys %sample_names)
	{
		if(!$sample_name_to_read_depth_file{$sample_name}
			and $within_sample_diversity_file_to_stage{$sample_name_to_within_sample_diversity_file{$sample_name}} ne "bam")
		{
			delete $sample_names{$sample_name};
		}
	}
	
	# prints number of samples remaining
	print_number_samples_remaining_and_exit_if_none();
}


# reads in reference sequence
open REFERENCE_FASTA, "<$reference_genome_file" || die "Could not open $reference_genome_file to read; terminating =(\n";
my $reference_sequence = "";
my $reference_sequence_name = "";
while(<REFERENCE_FASTA>) # for each line in the file
{
	chomp;
	if($_ =~ /^>(.*)/) # header line
	{
		if($reference_sequence_name)
		{
			print STDERR "Error: more than one sequence in reference sequence. Exiting.\n";
			die;
		}
		$reference_sequence_name = $1;
	}
	else # sequence
	{
		$reference_sequence .= $_;
	}
}
close REFERENCE_FASTA;

# counts unambiguous bases in reference
my $reference_sequence_length = count_unambiguous_bases_in_sequence(split(//, $reference_sequence));
if(!$reference_sequence_length)
{
	print STDERR "Error: reference sequence contains no unambiguous (A, T, C, G) bases:\n\t"
		.$reference_genome_file."\nExiting.\n";
	die;
}


# reads in consensus genomes; removes samples that do not have sufficiently complete genomes
my %sequence_name_to_consensus = (); # key: sequence name -> value: consensus sequence, including gaps froms alignment
if($minimum_genome_coverage)
{
	print STDERR "reading in consensus genomes...\n";
	foreach my $consensus_genome_fasta_file(@consensus_genome_files, $consensus_genomes_aligned_file)
	{
		if($consensus_genome_fasta_file)
		{
			open FASTA_FILE, "<$consensus_genome_fasta_file" || die "Could not open $consensus_genome_fasta_file to read; terminating =(\n";
			my $sequence = "";
			my $sample_name = "";
			while(<FASTA_FILE>) # for each line in the file
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
			close FASTA_FILE;
		}
	}
	
	print STDERR "removing samples without at least ".($minimum_genome_coverage*100)
		."% coverage ("."of ".$reference_sequence_length." total bases)...\n" if $verbose;
	remove_samples_without_minimum_genome_coverage();
	
	# prints number of samples remaining
	print_number_samples_remaining_and_exit_if_none();
	
	# clears consensus genomes
	%sequence_name_to_consensus = ();
}


# if a plate map is provided, reads in plate map positions of all samples
# removes samples that do not have neighbors
my %sample_name_to_all_plate_positions = (); # key: sample name -> value: string including all plate positions the sample appears in
my %sample_name_to_all_plates = (); # key: sample name -> value: string including all plates the sample appears in

my %plate_position_to_sample_name = (); # key: plate map file -> plate position -> value: sample name
my %sample_name_to_plate_position = (); # key: plate map file -> sample name -> value: plate position

if(scalar @plate_map_files)
{
	print STDERR "reading in plate map positions...\n" if $verbose;
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
					$plate_position_to_sample_name{$plate_map_file}{$position} = $sample_name;
					$sample_name_to_plate_position{$plate_map_file}{$sample_name} = $position;
					
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
	}
	
	# removes samples that don't have at least one plate neighbor
	if(!$print_all_isnvs)
	{
		print STDERR "removing samples without plate neighbors...\n" if $verbose;
		remove_samples_without_plate_neighbors();
	
		# prints number of samples remaining
		print_number_samples_remaining_and_exit_if_none();
	}
}


# generates and reads in read depth files
my %sample_name_to_position_to_read_depth = (); # key: sample name -> key: position -> value: read depth that that position
if($minimum_read_depth > 0)
{
	print STDERR "generating and reading in read depth files, using ".$cores_to_use." cores in parallel...\n" if $verbose;
	
	my $pm = Parallel::ForkManager -> new($cores_to_use);
	$pm -> run_on_finish(
	sub
	{
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
		my $q = $data_structure_reference -> {sample};
		$sample_name_to_position_to_read_depth{$q} = $data_structure_reference -> {position_to_read_depth};
	});

	foreach my $sample(keys %sample_names)
	{
		my $pid = $pm -> start and next;
		
		my $within_sample_diversity_file = $sample_name_to_within_sample_diversity_file{$sample};
		my %position_to_read_depth = (); # key: position -> value: read depth that that position
		
		# generates read depth file if needed
		my $read_depth_file = $sample_name_to_read_depth_file{$sample};
		if(!$read_depth_file)
		{
			if($within_sample_diversity_file_to_stage{$within_sample_diversity_file} eq "bam") # if this sample indeed has a bam file
			{
				# runs samtools depth
				$read_depth_file = $temp_intermediate_directory.retrieve_file_name($within_sample_diversity_file."_read_depth.txt");
				check_if_file_exists_before_writing($read_depth_file);
				print STDERR "$SAMTOOLS_EXECUTABLE_FILE_PATH depth $within_sample_diversity_file > $read_depth_file\n" if $verbose;
				`$SAMTOOLS_EXECUTABLE_FILE_PATH depth $within_sample_diversity_file > $read_depth_file`;
				print STDERR "\n" if $verbose;
			}
			else # no bam file
			{
				print STDERR "Error: expected bam file for sample ".$sample.". Could not "
					."generate read depth file.\n";
			}
		}
		
		# verifies that read depth file is generated and is not empty
		verify_input_file_exists_and_is_nonempty($read_depth_file, "read depth file", 1, 1);
		
		# reads in read depth file
		open READ_DEPTH_FILE, "<$read_depth_file" || die "Could not open $read_depth_file to read; terminating =(\n";
		while(<READ_DEPTH_FILE>) # for each line in the file
		{
			chomp;
			if($_ =~ /\S/)
			{
				# reads in mapped values
				my @items_in_row = split($DELIMITER, $_);
	
				my $position = $items_in_row[$READ_DEPTH_POSITION_COLUMN];
				my $read_depth = $items_in_row[$READ_DEPTH_COLUMN];
	
				$position_to_read_depth{$position} = $read_depth;
			}
		}
		close READ_DEPTH_FILE;
		$pm -> finish(0, {position_to_read_depth => \%position_to_read_depth, sample => $sample});
	}
	$pm -> wait_all_children;

# not parallelized; uncomment and replace parallelized version if needed
# 	print STDERR "generating read depth files (not parallelized)...\n" if $verbose;
# 
# 	foreach my $sample(keys %sample_names)
# 	{
# 		my $within_sample_diversity_file = $sample_name_to_within_sample_diversity_file{$sample};
# 		
# 		# generates read depth file if needed
# 		my $read_depth_file = $sample_name_to_read_depth_file{$sample};
# 		if(!$read_depth_file)
# 		{
# 			$read_depth_file = $temp_intermediate_directory.retrieve_file_name($within_sample_diversity_file."_read_depth.txt");
# 			if($within_sample_diversity_file_to_stage{$within_sample_diversity_file} eq "bam") # if this is indeed a bam file
# 			{
# 				# runs samtools depth
# 				check_if_file_exists_before_writing($read_depth_file);
# 				print STDERR "$SAMTOOLS_EXECUTABLE_FILE_PATH depth $within_sample_diversity_file > $read_depth_file\n" if $verbose;
# 				`$SAMTOOLS_EXECUTABLE_FILE_PATH depth $within_sample_diversity_file > $read_depth_file`;
# 				print STDERR "\n" if $verbose;
# 			}
# 			else
# 			{
# 				print STDERR "Error: expected bam file for sample ".$sample.". Could not "
# 					."generate read depth file.\n";
# 			}
# 		}
# 		
# 		# verifies that read depth file is generated and is not empty
# 		verify_input_file_exists_and_is_nonempty($read_depth_file, "read depth file", 1, 1);
# 		
# 		# reads in read depth file
# 		open READ_DEPTH_FILE, "<$read_depth_file" || die "Could not open $read_depth_file to read; terminating =(\n";
# 		while(<READ_DEPTH_FILE>) # for each line in the file
# 		{
# 			chomp;
# 			if($_ =~ /\S/)
# 			{
# 				# reads in mapped values
# 				my @items_in_row = split($DELIMITER, $_);
# 	
# 				my $position = $items_in_row[$READ_DEPTH_POSITION_COLUMN];
# 				my $read_depth = $items_in_row[$READ_DEPTH_COLUMN];
# 	
# 				$sample_name_to_position_to_read_depth{$sample}{$position} = $read_depth;
# 			}
# 		}
# 		close READ_DEPTH_FILE;
# 	}

	# verifies that each sample has at least minimum_genome_coverage * reference_sequence_length
	# positions with sufficient read depth
	print STDERR "removing samples without at least "
		.($minimum_genome_coverage * $reference_sequence_length)." unmasked bases with read depth >= "
		.$minimum_read_depth." (".(100*$minimum_genome_coverage)."% of "
		.$reference_sequence_length." total bases)...\n" if $verbose;
	my $samples_removed = remove_samples_without_minimum_genome_coverage_with_high_read_depth();
	
	# prints number of samples remaining
	print_number_samples_remaining_and_exit_if_none();
	
	if($samples_removed and scalar @plate_map_files and !$print_all_isnvs)
	{
		# removes samples that are now without plate neighbors
		print STDERR "removing samples without plate neighbors...\n" if $verbose;
		remove_samples_without_plate_neighbors();
	
		# prints number of samples remaining
		print_number_samples_remaining_and_exit_if_none();
	}
}


# aligns consensus genomes if they aren't already aligned
if(!$consensus_genomes_aligned_file)
{
	# reads in all consensus genomes
	print STDERR "retrieving included consensus genomes...\n" if $verbose;
	my %sequence_name_to_consensus = (); # key: sequence name -> value: consensus sequence, including gaps froms alignment
	foreach my $consensus_genome_fasta_file(@consensus_genome_files, $consensus_genomes_aligned_file)
	{
		if($consensus_genome_fasta_file)
		{
			open FASTA_FILE, "<$consensus_genome_fasta_file" || die "Could not open $consensus_genome_fasta_file to read; terminating =(\n";
			my $sequence = "";
			my $sample_name = "";
			while(<FASTA_FILE>) # for each line in the file
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
			close FASTA_FILE;
		}
	}
	
	# generates concatenated consensus genomes file
	print STDERR "generating concatenated consensus genomes fasta file...\n" if $verbose;
	my $all_consensus_genomes_file = $temp_intermediate_directory."all_consensus_genomes_concat.fasta";
	check_if_file_exists_before_writing($all_consensus_genomes_file);
	open CONSENSUS_GENOMES, ">$all_consensus_genomes_file" || die "Could not open $all_consensus_genomes_file to write; terminating =(\n";
	
	# prints reference sequence
	print CONSENSUS_GENOMES ">".$reference_sequence_name.$NEWLINE;
	print CONSENSUS_GENOMES $reference_sequence.$NEWLINE;
	
	# prints included consensus genomes
	foreach my $sample_name(keys %sample_names)
	{
		print CONSENSUS_GENOMES ">".$sample_name.$NEWLINE;
		print CONSENSUS_GENOMES $sequence_name_to_consensus{$sample_name}.$NEWLINE;
	}
	close CONSENSUS_GENOMES;
	
	# aligns all consensus genomes
	print STDERR "aligning consensus genomes...\n" if $verbose;
	$consensus_genomes_aligned_file = $temp_intermediate_directory."all_consensus_genomes_MAFFT_aligned.fasta";
	check_if_file_exists_before_writing($consensus_genomes_aligned_file);
	`$MAFFT_EXECUTABLE_FILE_PATH $all_consensus_genomes_file > $consensus_genomes_aligned_file`;
	
	# removes temp files
	`rm $all_consensus_genomes_file`;
}


# reads in aligned consensus genomes
print STDERR "reading in aligned consensus genomes...\n";
open ALIGNED_CONSENSUS_GENOMES, "<$consensus_genomes_aligned_file" || die "Could not open $consensus_genomes_aligned_file to read; terminating =(\n";
%sequence_name_to_consensus = (); # key: sequence name -> value: consensus sequence, including gaps froms alignment
$reference_sequence = ""; # first sequence in alignment

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
			$sequence_name_to_consensus{$sample_name} = uc($sequence);
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


# saves base indices in reference where there are gaps
print STDERR "identifying indices with gaps in reference and removing from all sequences in alignment...\n" if $verbose;
my %base_index_has_gap = (); # key: index of base in reference sequence (0-indexed) -> value: 1 if there is a gap in the reference sequence
my @reference_values = split(//, $reference_sequence);
for(my $base_index = 0; $base_index < length($reference_sequence); $base_index++)
{
	if(!is_base($reference_values[$base_index]))
	{
		$base_index_has_gap{$base_index} = 1;
	}
}

# removes bases or gaps at the corresponding positions in all other sequences in the alignment
foreach my $sample_name(keys %sequence_name_to_consensus)
{
	$sequence_name_to_consensus{$sample_name}
		= remove_bases_at_indices_with_gaps_in_reference($sequence_name_to_consensus{$sample_name});
}

# removes gaps in reference (first sequence) in alignment
$reference_sequence = remove_bases_at_indices_with_gaps_in_reference($reference_sequence);
@reference_values = split(//, $reference_sequence);


# masks any positions with read depth below minimum_read_depth
# masks any positions input by user
# removes samples that do not have sufficiently complete genomes with read depth filter applied
# removes samples without plate neighbors after read depth filter applied
my %sequence_name_to_pre_masking_consensus = ();
if($minimum_read_depth > 0 or scalar keys %position_is_masked)
{
	print STDERR "masking positions with read depth < ".$minimum_read_depth
		." and user-defined masked positions...\n" if $verbose;
	foreach my $sample_name(keys %sample_names)
	{
		# retrieves consensus genome bases
		my $consensus = $sequence_name_to_consensus{$sample_name};
		my @consensus_values = split(//, $consensus);
		
		# masks positions with low read depth
		if($minimum_read_depth > 0)
		{
			for my $position(keys %{$sample_name_to_position_to_read_depth{$sample_name}})
			{
				if($sample_name_to_position_to_read_depth{$sample_name}{$position} < $minimum_read_depth)
				{
					$consensus_values[$position - 1] = "N";
				}
			}
		}
		
		# masks any positions input by user
		if(scalar keys %position_is_masked)
		{
			foreach my $position(keys %position_is_masked)
			{
				if($position_is_masked{$position})
				{
					$consensus_values[$position - 1] = "N";
				}
			}
		}
		
		# saves masked and non-masked consensus genomes
		$sequence_name_to_pre_masking_consensus{$sample_name} = $consensus;
		$sequence_name_to_consensus{$sample_name} = join("", @consensus_values);
	}
	
	print STDERR "removing samples without at least ".($minimum_genome_coverage*100)
		."% coverage with read depth >= ".$minimum_read_depth." ("."of "
		.$reference_sequence_length." total bases)...\n" if $verbose;
	my $samples_removed = remove_samples_without_minimum_genome_coverage();
	
	# prints number of samples remaining
	print_number_samples_remaining_and_exit_if_none();
	
	if($samples_removed and scalar @plate_map_files and !$print_all_isnvs)
	{
		# removes samples that are now without plate neighbors
		print STDERR "removing samples without plate neighbors...\n" if $verbose;
		remove_samples_without_plate_neighbors();
	
		# prints number of samples remaining
		print_number_samples_remaining_and_exit_if_none();
	}
}
else
{
	# copies consensus genomes without masking
	%sequence_name_to_pre_masking_consensus = %sequence_name_to_consensus;
}


# pre-processes within-sample diversity files
# parallelization based on https://perlmaven.com/speed-up-calculation-by-running-in-parallel
print STDERR "pre-processing within-sample diversity files, using ".$cores_to_use." cores in parallel...\n" if $verbose;
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


# if plate map provided, counts number positions with heterozygosity (iSNVs) in each
# sample and generates visualization
# within-sample diversity files must be pre-processed before this step
my %plate_map_file_to_plate_iSNV_file = (); # key: plate map file -> value: file with plate iSNVs
if(scalar @plate_map_files)
{
	print STDERR "generating plate-map iSNV summaries...\n" if $verbose;
	my %sample_name_to_number_positions_with_heterozygosity = (); # key: sample name -> value: number iSNVs
	my %sample_examined = (); # key: sample name -> value: 1 if we've counted iSNVs in this sample
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
				my $plate_position = uc $items[$PLATE_MAP_POSITION_COLUMN];
				
				if($plate_position and $sample_name)
				{
					$plate_position_to_sample_name{$plate_map_file}{$plate_position} = $sample_name;
					$sample_name_to_plate_position{$plate_map_file}{$sample_name} = $plate_position;
				}
			}
		}
		close PLATE_MAP;

		# retrieves number iSNVs for each sample on plate map
		foreach my $plate_position(keys %{$plate_position_to_sample_name{$plate_map_file}})
		{
			my $sample = $plate_position_to_sample_name{$plate_map_file}{$plate_position};
			if($sample_names{$sample}) # if sample not excluded
			{
				# retrieves within-sample diversity file for potential contaminated sample
				my $within_sample_diversity_file = $sample_name_to_within_sample_diversity_file{$sample};
				if($within_sample_diversity_file_to_stage{$within_sample_diversity_file} ne "het")
				{
					print STDERR "Warning: within-sample diversity file for sample "
						.$sample." not preprocessed:\n\t".$within_sample_diversity_file."\n"
						."Processing without parallelization.\n";
		
					# if needed, processes within-sample diversity file for this sample
					$within_sample_diversity_file = process_within_sample_diversity_file_for_sample($sample);
					$sample_name_to_within_sample_diversity_file{$sample} = $within_sample_diversity_file;
					$within_sample_diversity_file_to_stage{$within_sample_diversity_file} = "het";
				}
			
				# reads in heterozygosity table generated for within-sample diversity of this sample
				if(!$sample_examined{$sample})
				{
					my %positions_with_heterozygosity = (); # key: position -> value: 1 if position has minor allele
					my $number_positions_with_heterozygosity = 0;
					open HETEROZYGOSITY_TABLE, "<$within_sample_diversity_file" || die "Could not open $within_sample_diversity_file to read; terminating =(\n";
					while(<HETEROZYGOSITY_TABLE>) # for each line in the file
					{
						chomp;
						my $line = $_;
	
						# parses this line
						my @items = split($DELIMITER, $line);
						my $position = $items[$HETEROZYGOSITY_TABLE_POSITION_COLUMN];
						my $minor_allele_readcount = $items[$HETEROZYGOSITY_TABLE_MINOR_ALLELE_READCOUNT_COLUMN];
						my $minor_allele_frequency = $items[$HETEROZYGOSITY_TABLE_MINOR_ALLELE_FREQUENCY_COLUMN];
	
						# verifies that input values all make sense
						if($position !~ /^\d+$/)
						{
							print STDERR "Warning: position is not non-zero integer"
								.$position." in heterozygosity table:\n\t"
								.$within_sample_diversity_file."\n";
						}
						if($minor_allele_readcount !~ /^\d+$/)
						{
							print STDERR "Warning: minor allele readcount is not non-zero integer"
								.$minor_allele_readcount." in heterozygosity table:\n\t"
								.$within_sample_diversity_file."\n";
						}
						if($minor_allele_frequency !~ /^[\de.-]+$/)
						{
							print STDERR "Warning: non-numerical minor allele frequency "
								.$minor_allele_frequency." in heterozygosity table:\n\t"
								.$within_sample_diversity_file."\n";
						}

						# only includes positions with minor allele readcount >= 10, minor allele frequency >= 3%
						# assumes that major allele frequency = 100% - minor allele frequency
						if($minor_allele_readcount >= $minimum_minor_allele_readcount
							and $minor_allele_frequency >= $minimum_minor_allele_frequency)
						{
							if($positions_with_heterozygosity{$position})
							{
								print STDERR "Warning: position ".$position
									." appears in more than one line in heterozygosity table:\n\t"
									.$within_sample_diversity_file."\n";
							}
							else
							{
								$positions_with_heterozygosity{$position} = 1;
								$number_positions_with_heterozygosity++;
							}
						}
					}
					close HETEROZYGOSITY_TABLE;
					
					# saves number iSNVs for this sample
					$sample_name_to_number_positions_with_heterozygosity{$sample} = $number_positions_with_heterozygosity;
					$sample_examined{$sample} = 1;
				}
			}
		}
		
		
		# prints copy of plate map with number iSNVs
		# creates file
		my $plate_iSNVs_output_file = $visualizations_directory.retrieve_file_name($plate_map_file)."_iSNVs.txt";
		check_if_file_exists_before_writing($plate_iSNVs_output_file);
		
		# prints header line
		open PLATE_ISNVS_OUT_FILE, ">$plate_iSNVs_output_file" || die "Could not open $plate_iSNVs_output_file to write; terminating =(\n";
		print PLATE_ISNVS_OUT_FILE "well".$DELIMITER;
		print PLATE_ISNVS_OUT_FILE "sample".$DELIMITER;
		print PLATE_ISNVS_OUT_FILE "iSNVs".$NEWLINE;
		
		# prints number iSNVs for each sample
		foreach my $plate_position(keys %{$plate_position_to_sample_name{$plate_map_file}})
		{
			my $sample_name = $plate_position_to_sample_name{$plate_map_file}{$plate_position};
			my $number_positions_with_heterozygosity = $sample_name_to_number_positions_with_heterozygosity{$sample_name};
			if(!$sample_names{$sample_name})
			{
				$number_positions_with_heterozygosity = $NO_DATA;
			}
			
			print PLATE_ISNVS_OUT_FILE $plate_position.$DELIMITER;
			print PLATE_ISNVS_OUT_FILE $sample_name.$DELIMITER;
			print PLATE_ISNVS_OUT_FILE $number_positions_with_heterozygosity.$NEWLINE;
		}

		# closes output file
		close PLATE_ISNVS_OUT_FILE;
		
		# generates visualization
		my $plate_iSNV_visualization_output_file = trim_off_file_extension($plate_iSNVs_output_file)."_visualization";
		check_if_file_exists_before_writing($plate_iSNV_visualization_output_file.".jpg");
		check_if_file_exists_before_writing($plate_iSNV_visualization_output_file.".pdf");
		
 		`$PLATE_VISUALIZATION_FILE_PATH $plate_iSNVs_output_file $plate_iSNV_visualization_output_file $plate_number_columns $plate_number_rows isnvs`;
	}
}


# removes samples without neighbors if we haven't already (if we delayed it to print all iSNVs passing thresholds)
if(scalar @plate_map_files and $print_all_isnvs)
{
	# removes samples without plate neighbors
	print STDERR "removing samples without plate neighbors...\n" if $verbose;
	remove_samples_without_plate_neighbors();

	# prints number of samples remaining
	print_number_samples_remaining_and_exit_if_none();
}


# prepares to process sample pairs in parallel
# parallelization based on https://perlmaven.com/speed-up-calculation-by-running-in-parallel
my %results = ();
my %plate_results = ();
$pm = Parallel::ForkManager -> new($cores_to_use);
$pm -> run_on_finish(
sub
{
	my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;

	my $sample_1 = $data_structure_reference -> {sample_1};
	my $sample_2 = $data_structure_reference -> {sample_2};
	
	$results{$sample_1}{$sample_2} = $data_structure_reference -> {result_a};
	$plate_results{$sample_1}{$sample_2} = $data_structure_reference -> {plate_result_a};

	$results{$sample_2}{$sample_1} = $data_structure_reference -> {result_b};
	$plate_results{$sample_2}{$sample_1} = $data_structure_reference -> {plate_result_b};
});

# if plate map(s) provided, compares all neighboring samples
# within-sample diversity files must be pre-processed before this step
if(scalar @plate_map_files)
{
	print STDERR "comparing all neighboring samples...\n" if $verbose;
	my %sample_pair_compared = (); # keys: sample name, sample name -> value: 1 if sample pair compared
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
				my $plate_position = uc $items[$PLATE_MAP_POSITION_COLUMN];
				
				if($plate_position and $sample_name and $sample_names{$sample_name})
				{
					$plate_position_to_sample_name{$plate_map_file}{$plate_position} = $sample_name;
					$sample_name_to_plate_position{$plate_map_file}{$sample_name} = $plate_position;
				}
			}
		}
		close PLATE_MAP;
		
		# compares all pairs of samples that are neighbors
		foreach my $plate_position(keys %{$plate_position_to_sample_name{$plate_map_file}})
		{
			my $sample_1 = $plate_position_to_sample_name{$plate_map_file}{$plate_position};
			if($sample_names{$sample_1}) # verify that this sample is included
			{
				# retrieve neighboring samples (only samples that are included)
				my @neighboring_samples = retrieve_samples_neighboring_plate_position($plate_position, $sample_1, $plate_map_file);
			
				# compares sample to all included neighboring samples
				foreach my $sample_2(@neighboring_samples)
				{
					my $neighboring_plate_position = $sample_name_to_plate_position{$plate_map_file}{$sample_2};
					if(!$sample_pair_compared{$sample_1}{$sample_2}) # this pair of samples not already compared
					{
						# records that we are comparing this pair (to avoid comparing it again)
						$sample_pair_compared{$sample_1}{$sample_2} = 1;
						$sample_pair_compared{$sample_2}{$sample_1} = 1;
				
						# tests sample 1 contaminating sample 2 and vice versa
						my $pid = $pm -> start and next;
					
						my ($result_a, $plate_result_a) = detect_potential_contamination_in_sample_pair($sample_1, $sample_2);
						my ($result_b, $plate_result_b) = detect_potential_contamination_in_sample_pair($sample_2, $sample_1);
					
						$pm -> finish(0, {sample_1 => $sample_1, sample_2 => $sample_2,
							result_a => $result_a, plate_result_a => $plate_result_a,
							result_b => $result_b, plate_result_b => $plate_result_b});
					}
				}
			}
		}
		$pm -> wait_all_children;
	}
	
	# collects outputs for each individual plate map
	my %plate_map_file_to_output = (); # key: plate map file -> value: plate-specific cross-contamination output to print
	foreach my $sample_1(sort keys %plate_results)
	{
		foreach my $sample_2(sort keys %{$plate_results{$sample_1}})
		{
			if($plate_results{$sample_1}{$sample_2})
			{
				foreach my $plate_map_file(@plate_map_files)
				{
					if($sample_name_to_plate_position{$plate_map_file}{$sample_1}
						and $sample_name_to_plate_position{$plate_map_file}{$sample_2})
					{
						$plate_map_file_to_output{$plate_map_file} .= $sample_name_to_plate_position{$plate_map_file}{$sample_1}.$DELIMITER;
						$plate_map_file_to_output{$plate_map_file} .= $sample_name_to_plate_position{$plate_map_file}{$sample_2}.$DELIMITER;
						$plate_map_file_to_output{$plate_map_file} .= $plate_results{$sample_1}{$sample_2};
						$plate_map_file_to_output{$plate_map_file} .= $NEWLINE;
					}
				}
			}
		}
	}
	
	# creates cross contamination table and visualization for each individual plate map
	foreach my $plate_map_file(@plate_map_files)
	{
		if($plate_map_file_to_output{$plate_map_file}) # if there is anything to print
		{
			my $plate_output_file = $visualizations_directory.retrieve_file_name($plate_map_file)."_potential_cross_contamination.txt";
			check_if_file_exists_before_writing($plate_output_file);
			open PLATE_OUT_FILE, ">$plate_output_file" || die "Could not open $plate_output_file to write; terminating =(\n";
	
			# prints header line
			print PLATE_OUT_FILE "well".$DELIMITER;
			print PLATE_OUT_FILE "contamination_source_well".$DELIMITER;
			print PLATE_OUT_FILE "sample".$DELIMITER;
			print PLATE_OUT_FILE "contamination_source_sample".$DELIMITER;
			print PLATE_OUT_FILE "appearance_of_potential_contamination".$DELIMITER;
			print PLATE_OUT_FILE "estimated_contamination_volume".$NEWLINE;
	
			# prints cross contamination table
			print PLATE_OUT_FILE $plate_map_file_to_output{$plate_map_file};
	
			# closes output file
			close PLATE_OUT_FILE;
	
			# generates visualization
			my $plate_visualization_output_file = trim_off_file_extension($plate_output_file)."_visualization";
			check_if_file_exists_before_writing($plate_visualization_output_file.".jpg");
			check_if_file_exists_before_writing($plate_visualization_output_file.".pdf");
			`$PLATE_VISUALIZATION_FILE_PATH $plate_output_file $plate_visualization_output_file $plate_number_columns $plate_number_rows contamination`;

			# generates visualization of consensus-level contamination only
			my $plate_visualization_output_file_consensus = trim_off_file_extension($plate_output_file)."_visualization_consensus";
			check_if_file_exists_before_writing($plate_visualization_output_file_consensus.".jpg");
			check_if_file_exists_before_writing($plate_visualization_output_file_consensus.".pdf");
			`$PLATE_VISUALIZATION_FILE_PATH $plate_output_file $plate_visualization_output_file_consensus $plate_number_columns $plate_number_rows contamination_consensus`;

			# generates visualization of minor allele-level contamination only
			my $plate_visualization_output_file_minor = trim_off_file_extension($plate_output_file)."_visualization_minor";
			check_if_file_exists_before_writing($plate_visualization_output_file_minor.".jpg");
			check_if_file_exists_before_writing($plate_visualization_output_file_minor.".pdf");
			`$PLATE_VISUALIZATION_FILE_PATH $plate_output_file $plate_visualization_output_file_minor $plate_number_columns $plate_number_rows contamination_minor`;

		}
	}
}

# if no plate map provided, compares all pairs of samples
else
{
	print STDERR "comparing all pairs of samples...\n" if $verbose;
	my @sample_names_array = keys %sample_names;
	for(my $index_1 = 0; $index_1 <= $#sample_names_array; $index_1++)
	{
		my $sample_1 = $sample_names_array[$index_1];
		for(my $index_2 = $index_1+1; $index_2 <= $#sample_names_array; $index_2++)
		{
			my $sample_2 = $sample_names_array[$index_2];
			
			# tests sample 1 contaminating sample 2 and vice versa
			my $pid = $pm -> start and next;

			my ($result_a, $plate_result_a) = detect_potential_contamination_in_sample_pair($sample_1, $sample_2);
			my ($result_b, $plate_result_b) = detect_potential_contamination_in_sample_pair($sample_2, $sample_1);
			
			$pm -> finish(0, {sample_1 => $sample_1, sample_2 => $sample_2,
				result_a => $result_a, plate_result_a => $plate_result_a,
				result_b => $result_b, plate_result_b => $plate_result_b});
		}
	}
	$pm -> wait_all_children;
}


# creates directory for output file if it does not already exist
my $output_file_directory = retrieve_file_directory($output_file_path);
if(!-d $output_file_directory and -e $output_file_directory)
{
	# directory already exists and is a file, not a directory
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

# opens output file
check_if_file_exists_before_writing($output_file_path);
open OUT_FILE, ">$output_file_path" || die "Could not open $output_file_path to write; terminating =(\n";

# prints header line
print OUT_FILE "potential_contaminated_sample".$DELIMITER;
print OUT_FILE "potential_contaminated_sample_unambiguous_bases".$DELIMITER;
print OUT_FILE "potential_contaminated_sample_genome_covered".$DELIMITER;
if($minimum_read_depth)
{
	print OUT_FILE "potential_contaminated_sample_unambiguous_bases_passing_read_depth_filter".$DELIMITER;
	print OUT_FILE "potential_contaminated_sample_genome_covered_passing_read_depth_filter".$DELIMITER;
}
print OUT_FILE "num_positions_with_heterozygosity".$DELIMITER;
print OUT_FILE "alleles_at_positions_with_heterozygosity".$DELIMITER;

print OUT_FILE "potential_contaminating_sample".$DELIMITER;
print OUT_FILE "potential_contaminating_sample_unambiguous_bases".$DELIMITER;
print OUT_FILE "potential_contaminating_sample_genome_covered".$DELIMITER;
if($minimum_read_depth)
{
	print OUT_FILE "potential_contaminating_sample_unambiguous_bases_passing_read_depth_filter".$DELIMITER;
	print OUT_FILE "potential_contaminating_sample_genome_covered_passing_read_depth_filter".$DELIMITER;
}
print OUT_FILE "minor_alleles_matched".$DELIMITER;
print OUT_FILE "major_alleles_matched".$DELIMITER;
print OUT_FILE "heterozygous_positions_matched".$DELIMITER;
print OUT_FILE "alleles_matched".$DELIMITER;
print OUT_FILE "num_mismatches".$DELIMITER;
print OUT_FILE "mismatches".$DELIMITER;

print OUT_FILE "appearance_of_potential_contamination".$DELIMITER;
print OUT_FILE "estimated_contamination_volume".$DELIMITER;
print OUT_FILE "contaminating_allele_frequency_range".$DELIMITER;
print OUT_FILE "contaminating_allele_frequencies";

if(scalar @plate_map_files)
{
	print OUT_FILE $DELIMITER;
	print OUT_FILE "potential_contaminated_sample_plate_position".$DELIMITER;
	print OUT_FILE "potential_contaminating_sample_plate_position";
	if(scalar @plate_map_files > 1)
	{
		print OUT_FILE $DELIMITER;
		print OUT_FILE "potential_contaminated_sample_plate";
		print OUT_FILE $DELIMITER;
		print OUT_FILE "potential_contaminating_sample_plate";
	}
}
print OUT_FILE $NEWLINE;

# prints output line for each pair of samples with potential contamination
foreach my $sample_1(sort keys %results)
{
    foreach my $sample_2(sort keys %{$results{$sample_1}})
    {
    	if($results{$sample_1}{$sample_2})
    	{
			print OUT_FILE $results{$sample_1}{$sample_2};
			print OUT_FILE $NEWLINE;
    	}
    }
}

# closes output file
close OUT_FILE;


# HELPER FUNCTIONS FOR DETECTING POTENTIAL CONTAMINATION

# compares these two samples to check for contamination in both directions
# sub detect_potential_contamination_in_sample_pair_both_directions
# {
# 	my $sample_1 = $_[0];
# 	my $sample_2 = $_[1];
# 	
# 	my ($result_1, $result_1_plate) = detect_potential_contamination_in_sample_pair($sample_1, $sample_2);
# 	my ($result_2, $result_2_plate) = detect_potential_contamination_in_sample_pair($sample_2, $sample_1);
# 	
# 	# returns output lines
# 	if($result_1 and $result_2)
# 	{
# 		return ($result_1.$NEWLINE.$result_2, $result_1_plate.$NEWLINE.$result_2_plate);
# 	}
# 	elsif($result_1)
# 	{
# 		return ($result_1, $result_1_plate);
# 	}
# 	elsif($result_2)
# 	{
# 		return ($result_2, $result_2_plate);
# 	}
# }

# prepares and retrieves files for the two input samples; checks if the second sample
# potentially could have contaminated the first
# assumes that within-sample diversity file has been pre-processed
sub detect_potential_contamination_in_sample_pair
{
	my $potential_contaminated_sample = $_[0];
	my $potential_contaminating_sample = $_[1];
	
	# exits if either contaminated or contaminating sample does not have a consensus genome
	if(!defined $sequence_name_to_consensus{$potential_contaminated_sample})
	{
		print STDERR "Error: no consensus genome for sample "
			.$potential_contaminated_sample.". Skipping comparison:\n\t"
			."potential contaminated:  ".$potential_contaminated_sample."\n\t"
			."potential contaminating: ".$potential_contaminating_sample."\n";
		return;
	}
	if(!defined $sequence_name_to_consensus{$potential_contaminating_sample})
	{
		print STDERR "Error: no consensus genome for sample "
			.$potential_contaminating_sample.". Skipping comparison:\n\t"
			."potential contaminated:  ".$potential_contaminated_sample."\n\t"
			."potential contaminating: ".$potential_contaminating_sample."\n";
		return;
	}
	
	# exits if contaminated sample does not have a within-sample diversity file
	if(!defined $sample_name_to_within_sample_diversity_file{$potential_contaminated_sample})
	{
		print STDERR "Error: no within-sample diversity file for sample "
			.$potential_contaminated_sample.". Skipping comparison:\n\t"
			."potential contaminated:  ".$potential_contaminated_sample."\n\t"
			."potential contaminating: ".$potential_contaminating_sample."\n";
		return;
	}
	
	# retrieves within-sample diversity file for potential contaminated sample
	my $potential_contaminated_within_sample_diversity_file = $sample_name_to_within_sample_diversity_file{$potential_contaminated_sample};
	
	# exits if within-sample diversity table is not at heterozygosity table stage
	if($within_sample_diversity_file_to_stage{$potential_contaminated_within_sample_diversity_file} ne "het")
	{
		print STDERR "Error: within-sample diversity file for sample "
			.$potential_contaminated_sample." not preprocessed:\n\t"
			.$potential_contaminated_within_sample_diversity_file."\nExiting.\n";
		die;
		
		# if needed, processes within-sample diversity file for potential contaminated sample
		# may cause race conditions and is in general a terrible idea--delete if possible
# 		my $within_sample_diversity_file = process_within_sample_diversity_file_for_sample($potential_contaminated_sample);
# 		$sample_name_to_within_sample_diversity_file{$potential_contaminated_sample} = $within_sample_diversity_file;
# 		$within_sample_diversity_file_to_stage{$within_sample_diversity_file} = "het";
# 		
# 		# updates value
# 		$potential_contaminated_within_sample_diversity_file = $sample_name_to_within_sample_diversity_file{$potential_contaminated_sample};
	}
	
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
	
		# only includes positions not masked by user,
		# with minor allele readcount >= 10, minor allele frequency >= 3%,
		# minor + consensus-level allele readcount >= minimum_read_depth,
		# and read depth in potential contaminated sample >= minimum_read_depth
		# assumes that major allele frequency = 100% - minor allele frequency
		if($minor_allele_readcount >= $minimum_minor_allele_readcount
			and $minor_allele_frequency >= $minimum_minor_allele_frequency
			and $minor_allele_readcount + $major_allele_readcount >= $minimum_read_depth
			and (!$minimum_read_depth or $sample_name_to_position_to_read_depth{$potential_contaminated_sample}{$position} >= $minimum_read_depth)
			and !$position_is_masked{$position})
		{
			if($positions_with_heterozygosity{$position})
			{
				print STDERR "Warning: position ".$position." appears in more than one line in "
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
	
	# does the same for the sequences before masking
	my $pre_masking_potential_contaminated_consensus = $sequence_name_to_pre_masking_consensus{$potential_contaminated_sample};
	my @pre_masking_potential_contaminated_consensus_values = split(//, $pre_masking_potential_contaminated_consensus);
	my $pre_masking_potential_contaminated_consensus_unambig_bases = count_unambiguous_bases_in_sequence(@pre_masking_potential_contaminated_consensus_values);
	my $pre_masking_potential_contaminated_consensus_percent_covered = $pre_masking_potential_contaminated_consensus_unambig_bases / $reference_sequence_length;
	
	my $pre_masking_potential_contaminating_consensus = $sequence_name_to_pre_masking_consensus{$potential_contaminating_sample};
	my @pre_masking_potential_contaminating_consensus_values = split(//, $pre_masking_potential_contaminating_consensus);
	my $pre_masking_potential_contaminating_consensus_unambig_bases = count_unambiguous_bases_in_sequence(@pre_masking_potential_contaminating_consensus_values);
	my $pre_masking_potential_contaminating_consensus_percent_covered = $pre_masking_potential_contaminating_consensus_unambig_bases / $reference_sequence_length;
	
	# exits if either sample does not have minimum genome coverage
	if($potential_contaminated_consensus_percent_covered < $minimum_genome_coverage
		or $potential_contaminating_consensus_percent_covered < $minimum_genome_coverage)
	{
		if(!$print_all)
		{
			return;
		}
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
		if(is_unambiguous_base($nucleotide_at_position)
			and (!$minimum_read_depth or $sample_name_to_position_to_read_depth{$potential_contaminating_sample}{$position} >= $minimum_read_depth)
			and (!$minimum_read_depth or $sample_name_to_position_to_read_depth{$potential_contaminated_sample}{$position} >= $minimum_read_depth)
			and !$position_is_masked{$position})
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
					if(!$print_all)
					{
						return;
					}
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
			
			# position not masked by user
			and !$position_is_masked{$position}
		
			 # mismatch between contaminated and contaminating consensus sequences
			and $potential_contaminated_consensus_value ne $potential_contaminating_consensus_value
			
			# bases in both contaminating and contaminated sequence are A, T, C, or G
			and is_unambiguous_base($potential_contaminating_consensus_value)
			and is_unambiguous_base($potential_contaminated_consensus_value)
			
			# both contaminating and contaminated sequence have read depth >= minimum_read_depth
			# at this position
			and (!$minimum_read_depth or $sample_name_to_position_to_read_depth{$potential_contaminating_sample}{$position}
				and $sample_name_to_position_to_read_depth{$potential_contaminating_sample}{$position} >= $minimum_read_depth)
			and (!$minimum_read_depth or $sample_name_to_position_to_read_depth{$potential_contaminated_sample}{$position}
				and $sample_name_to_position_to_read_depth{$potential_contaminated_sample}{$position} >= $minimum_read_depth))
		{
			# this base does not match consensus genome
			$number_mismatches++;
			
			# verifies that this potential contamination scenario is still worth looking at
			if($number_mismatches > $maximum_allowed_mismatches)
			{
				if(!$print_all)
				{
					return;
				}
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
		elsif(scalar @matched_allele_frequencies % 2 == 0) # even number of values
		{
			# average of the two center values
			$median_frequency = ($matched_allele_frequencies[(scalar @matched_allele_frequencies)/2 - 1]
				+ $matched_allele_frequencies[(scalar @matched_allele_frequencies)/2]) / 2;
		}
		else # odd number of values
		{
			$median_frequency = $matched_allele_frequencies[(scalar @matched_allele_frequencies - 1)/2]
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
	if($minimum_read_depth or scalar keys %position_is_masked)
	{
		$output_line .= add_comma_separators($pre_masking_potential_contaminated_consensus_unambig_bases).$DELIMITER;
		$output_line .= prepare_percentage_to_print($pre_masking_potential_contaminated_consensus_percent_covered).$DELIMITER;
	}
	$output_line .= add_comma_separators($potential_contaminated_consensus_unambig_bases).$DELIMITER;
	$output_line .= prepare_percentage_to_print($potential_contaminated_consensus_percent_covered).$DELIMITER;
	$output_line .= $number_positions_with_heterozygosity.$DELIMITER;
	$output_line .= $list_of_alleles_to_print.$DELIMITER; # TODO
	
	# adds columns about contaminating sample
	$output_line .= $potential_contaminating_sample.$DELIMITER;
	if($minimum_read_depth or scalar keys %position_is_masked)
	{
		$output_line .= add_comma_separators($pre_masking_potential_contaminating_consensus_unambig_bases).$DELIMITER;
		$output_line .= prepare_percentage_to_print($pre_masking_potential_contaminating_consensus_percent_covered).$DELIMITER;
	}
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
		$output_line_plate .= $potential_contaminated_sample.$DELIMITER;
		$output_line_plate .= $potential_contaminating_sample.$DELIMITER;
		$output_line_plate .= $contamination_type.$DELIMITER;
		$output_line_plate .= $median_frequency;
	}

	return ($output_line, $output_line_plate);
}

# finishes processing within-sample diversity file for sample, starting at whatever stage
# we have input file for:
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
# 		print STDERR "$LOFREQ_EXECUTABLE_FILE_PATH call -f $reference_genome_file -o $output_vcf_file $within_sample_diversity_file\n" if $verbose;
		`$LOFREQ_EXECUTABLE_FILE_PATH call -f $reference_genome_file -o $output_vcf_file $within_sample_diversity_file`;
		
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
# 		print STDERR "$VCF_TO_HETEROZYGOSITY_TABLE_SCRIPT_FILE_PATH $within_sample_diversity_file > $output_heterozygosity_table\n" if $verbose;
		`$VCF_TO_HETEROZYGOSITY_TABLE_SCRIPT_FILE_PATH $within_sample_diversity_file > $output_heterozygosity_table`;
		print STDERR "\n" if $verbose;
		
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

# removes bases at indices that have gaps in reference and prints modified sequence
# returns updated sequence
# reference sequence must be read into $reference_sequence
# indices with gaps in reference must be marked in %base_index_has_gap
sub remove_bases_at_indices_with_gaps_in_reference
{
	my $sequence = $_[0]; # must be capitalized

	# exit if no sequence input
	if(!$sequence)
	{
		return;
	}
	
	# exit if reference sequence not read in
	if(!$reference_sequence)
	{
		return;
	}
	
	# removes bases at indices that have gaps in reference
	my @sequence_values = split(//, $sequence);
	my @updated_sequence_values = ();
	for(my $base_index = 0; $base_index < length($sequence); $base_index++)
	{
		if(!$base_index_has_gap{$base_index})
		{
			push(@updated_sequence_values, $sequence_values[$base_index]);
		}
	}
	
	# generates and returns updated sequence from array
	return join("", @updated_sequence_values);
}

# returns 1 if base is not gap, 0 if base is a gap (or whitespace or empty)
sub is_base
{
	my $base = $_[0];
	
	# empty value
	if(!$base)
	{
		return 0;
	}
	
	# only whitespace
	if($base !~ /\S/)
	{
		return 0;
	}
	
	# gap
	if($base eq "-")
	{
		return 0;
	}
	
	# base
	return 1;
}


# HELPER FUNCTIONS FOR FILTERING SAMPLES

# prints number of samples remaining and exits if no samples remain
sub print_number_samples_remaining_and_exit_if_none
{
	# prints number of samples remaining
	print STDERR (keys %sample_names)." samples remain...\n" if $verbose;

	# verifies that we still have samples to compare
	if(scalar keys %sample_names < 2)
	{
		print STDERR "Error: no pairs of samples to compare. Exiting.\n";
		die;
	}
}

# removes samples that do not have genome coverage >= minimum_genome_coverage
# assumes that consensus genomes have been read into %sequence_name_to_consensus
# returns number of samples removed
sub remove_samples_without_minimum_genome_coverage
{
	my $number_samples_removed = 0;
	foreach my $sample_name(keys %sample_names)
	{
		# retrieves consensus genome bases
		my $consensus = $sequence_name_to_consensus{$sample_name};
		my @consensus_values = split(//, $consensus);

		# counts unambiguous bases in consensus genome
		my $consensus_unambig_bases = count_unambiguous_bases_in_sequence(@consensus_values);
		my $consensus_percent_covered = $consensus_unambig_bases / $reference_sequence_length;

		# removes sample name if sample does not have minimum genome coverage
		if($consensus_percent_covered < $minimum_genome_coverage)
		{
			delete $sample_names{$sample_name};
			$number_samples_removed++;
		}
	}
	
	return $number_samples_removed;
}

# removes samples without enough positions with sufficient read depth to reach
# minimum_genome_coverage
# assumes that sample_name_to_position_to_read_depth has been read in
# returns number of samples removed
sub remove_samples_without_minimum_genome_coverage_with_high_read_depth
{
	my $number_samples_removed = 0;
	for my $sample_name(keys %sample_name_to_position_to_read_depth)
	{
		my $number_positions_with_high_read_depth = 0;
		for my $position(keys %{$sample_name_to_position_to_read_depth{$sample_name}})
		{
			if($sample_name_to_position_to_read_depth{$sample_name}{$position} >= $minimum_read_depth
				and !$position_is_masked{$position})
			{
				$number_positions_with_high_read_depth++;
			}
		}
		
		if($number_positions_with_high_read_depth < $minimum_genome_coverage * $reference_sequence_length)
		{
			delete $sample_names{$sample_name};
			$number_samples_removed++;
		}
	}
	
	return $number_samples_removed;
}

# removes samples that don't have at least one plate neighbor
# assumes that %plate_position_to_sample_name and %sample_name_to_plate_position have been read in
# returns number of samples removed
sub remove_samples_without_plate_neighbors
{
	# checks if each sample has at least one neighbor
	my %sample_has_plate_neighbors = (); # key: sample name -> value: 1 if sample has plate neighbors
	foreach my $plate_map_file(@plate_map_files)
	{
		foreach my $plate_position(keys %{$plate_position_to_sample_name{$plate_map_file}})
		{
			my $sample_name = $plate_position_to_sample_name{$plate_map_file}{$plate_position};
			my @neighboring_samples = retrieve_samples_neighboring_plate_position($plate_position, $sample_name, $plate_map_file);
			
			if(scalar @neighboring_samples)
			{
				$sample_has_plate_neighbors{$sample_name} = 1;
			}
		}
	}
	
	# removes samples that don't have at least one plate neighbor
	my $number_samples_removed = 0;
	foreach my $sample_name(keys %sample_names)
	{
		if(!$sample_has_plate_neighbors{$sample_name})
		{
			delete $sample_names{$sample_name};
			$number_samples_removed++;
		}
	}
	
	return $number_samples_removed;
}


# HELPER FUNCTIONS FOR HANDLING PLATE NEIGHBORS

# retrieves samples neighboring given plate position
# assumes that %plate_position_to_sample_name and %sample_name_to_plate_position have been read in
# for this particular plate
# input: plate position (example input: H9)
# output: list of samples at neighboring plate positions
# (example output: sample at H8, sample at H10, sample at G9, sample at I9)
sub retrieve_samples_neighboring_plate_position
{
	my $plate_position = $_[0];
	my $plate_position_sample_name = $_[1];
	my $plate_map_file = $_[2];
	
	my @neighbors = ();
	if($plate_position =~ /^([A-Z]+)(\s*)(\d+)$/)
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
			foreach my $sample_name(keys %{$sample_name_to_plate_position{$plate_map_file}})
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
					if($plate_position_to_sample_name{$plate_map_file}{$neighbor_position})
					{
						push(@neighbors, $plate_position_to_sample_name{$plate_map_file}{$neighbor_position});
					}
				}
			}
			
			if($compare_row)
			{
				# retrieves all samples on same row (e.g., row A) on plate
				# compares this plate position's letter to letters of all plate positions
				foreach my $sample_name(keys %{$sample_name_to_plate_position{$plate_map_file}})
				{
					my $plate_position_option = $sample_name_to_plate_position{$plate_map_file}{$sample_name};
					if($plate_position_option =~ /^([A-Z+])\s*\d+$/)
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
				foreach my $sample_name(keys %{$sample_name_to_plate_position{$plate_map_file}})
				{
					my $plate_position_option = $sample_name_to_plate_position{$plate_map_file}{$sample_name};
					if($plate_position_option =~ /^[A-Z]+\s*(\d+)$/)
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
	
	# remove duplicate neighbors, empty neighbors, samples we're not looking at,
	# and neighbors identical to sample we're getting neighbors of
	my %neighbors_hash = ();
	foreach my $neighbor(@neighbors)
	{
		if($neighbor and $sample_names{$neighbor} and $neighbor ne $plate_position_sample_name)
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
	
	if($plate_position =~ /^[A-Z]+\s*\d+$/)
	{
		return 1;
	}
	return 0;
}

# returns next letter(s)
# A -> B, B -> C, Z -> AA, AB -> AC, AZ -> BA, ZZ -> AAA, etc.
sub get_next_letter
{
	my $letters_string = $_[0];
	$letters_string++; # thank you perl!
	return $letters_string;
}

# returns previous letter
# returns empty string if A
# B -> A, C -> B, AC -> AB, BA -> AZ, AA -> Z, AAA -> ZZ
sub get_previous_letter
{
	my $letters_string = $_[0];
	my @letters = split(//, $letters_string);
	for(my $letter_index = $#letters; $letter_index >= 0; $letter_index--) # starts at rightmost letter
	{
		my $letter = $letters[$letter_index];
		if($LETTER_TO_PREVIOUS_LETTER{$letter}) # not A
		{
			# we can just decrement this letter and return what we have
			$letters[$letter_index] = $LETTER_TO_PREVIOUS_LETTER{$letter};
			return join("", @letters);
		}
		else # A
		{
			if($letter_index == 0) # this is the leftmost letter
			{
				# this letter disappears
				$letters[$letter_index] = "";
			}
			
			else # there are more letters to the left
			{
				# this letter turns to Z
				$letters[$letter_index] = "Z";
			
				# and we decrement the letter to the left (next iteration through the loop!)
			}
		}
	}
	return join("", @letters);
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
			print STDERR "To allow file overwriting, use option --overwrite TRUE. Exiting.\n";
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

# prepares directory: adds / to end of path, verifies that directory is not a file,
# creates directory if it does not already exist
# returns directory path
sub prepare_directory
{
	my $directory = $_[0];

	# adds / to end of directory path if it isn't already there
	if($directory !~ /\/$/) # if doesn't end in /
	{
		$directory .= "/"; # adds / to the end
	}

	# creates directory for temporary and intermediate files if it doesn't already exist
	if(-e $directory and -d $directory)
	{
		# directory already exists
		print STDERR "Warning: adding files to already existing directory:\n\t"
			.$directory."\n";
	}
	elsif(-e $directory)
	{
		# directory already exists and is a file
		print STDERR "Error: directory already exists and is a file:\n\t"
			.$directory."\nExiting.\n";
		die;
	}
	else
	{
		# directory doesn't exist
		# create directory and all necessary parent directories
		`mkdir -p $directory`;
	}
	
	return $directory;
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
			print STDERR "Warning: ignoring ".$argument_option_1." | ".$argument_option_2
				." argument with no input file. Setting to empty string.\n";
			return "";
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
			print STDERR "Warning: ignoring ".$argument_option_1." | ".$argument_option_2
				." argument with no int. Setting to 0.\n";
			return 0;
		}
		else
		{
			$argument_index++;
			return $next_item;
		}
	}
	return -1; # this argument was not entered
}

# checks if parameter argument was supplied and, if so, reads in string that follows
# returns string supplied after argument
# returns -1 if this argument was not entered
sub read_in_string_argument
{
	my $argument_option_1 = $_[0]; # for example: -p
	my $argument_option_2 = $_[1]; # for example: --cores
	
	if($argument eq $argument_option_1 or $argument eq $argument_option_2)
	{
		my $next_item = $ARGV[$argument_index+1];
		if($argument_index + 1 > $#ARGV or $next_item !~ /^.+$/)
		{
			print STDERR "Warning: ignoring ".$argument_option_1." | ".$argument_option_2
				." argument with no string. Setting to empty string.\n";
			return "";
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
			print STDERR "Warning: ignoring ".$argument_option_1." | ".$argument_option_2
				." argument with no float. Setting to 0.\n";
			return 0;
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
	my $argument_option_1 = $_[0]; # short argument, case sensitive, for example: -d
	my $argument_option_2 = $_[1]; # long argument, not case sensitive, for example: --diagonal
	
	if($argument eq $argument_option_1 or uc($argument) eq uc($argument_option_2))
	{
		my $next_item = $ARGV[$argument_index+1];
		if($argument_index + 1 > $#ARGV or $next_item =~ /^-\w$/ or $next_item =~ /^--[\w-]+$/)
		{
			print STDERR "Warning: ignoring ".$argument_option_1." | ".$argument_option_2
				." argument with no boolean. Setting to false.\n";
			return 0;
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
			print STDERR "Warning: ignoring ".$argument_option_1." | ".$argument_option_2
				." argument with no input file. Setting to empty array.\n";
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

# prints error or warning if input file does not exist or is empty
# if exit_if_nonexistent is true (1), exits if input file does not exist
# if exit_if_empty is true (1), exits if input file is empty
sub verify_input_file_exists_and_is_nonempty
{
	my $input_file_path = $_[0]; # must be non-empty string
	my $input_file_name = $_[1]; # for printing error to console
	my $exit_if_nonexistent = $_[2]; # if 1, prints error and exits if file does not exist; if 0, prints warning and does not exit
	my $exit_if_empty = $_[3]; # if 1, prints error and exits if file is empty; if 0, prints warning and does not exit
	
	# verifies that input file exists
	if(!-e $input_file_path)
	{
		if($exit_if_nonexistent)
		{
			print STDERR "Error: ";
		}
		else
		{
			print STDERR "Warning: ";
		}
		print STDERR $input_file_name." does not exist:\n\t".$input_file_path."\n";
		if($exit_if_nonexistent)
		{
			print STDERR "Exiting.\n";
			die;
		}
	}
	
	# verifies that input file is non-empty
	if(-z $input_file_path)
	{
		if($exit_if_empty)
		{
			print STDERR "Error: ";
		}
		else
		{
			print STDERR "Warning: ";
		}
		print STDERR $input_file_name." is empty:\n\t".$input_file_path."\n";
		if($exit_if_empty)
		{
			print STDERR "Exiting.\n";
			die;
		}
	}
}

# June 1, 2021
