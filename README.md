# polyphonia
Detects cross-contamination in sequence data.

Usage: `perl detect_potential_cross_contamination.pl [options]`

```
OPTIONS:
- Reference (required):
	-f | --ref FILE			Reference fasta file [null]

- Consensus genomes (aligned or not aligned, not both; at least one file required):
	-c | --consensus FILE(S)	Unaligned consensus genome or genomes [null]
	-a | --consensus-aligned FILE	Consensus genomes pre-aligned to reference as fasta alignment; reference provided by --ref must be first [null]
	-g | --min-covered FLOAT	Minimum proportion genome covered for a sample to be included [0.98]

- Within-sample diversity (any combination; at least one file required):
	-b | --bam FILE(S)		Aligned and trimmed reads as bam file(s); must use reference provided by --ref [null]
	-v | --vcf FILE(S)		VCF file(s) output by LoFreq or GATK; must use reference provided by --ref [null]
	-h | --het FILE(S)		Tab-separated heterozygosity summary tables; see documentation for format [null]
	-e | --min-readcount INT	Minimum minor allele readcount for position to be considered heterozygous [10]
	-i | --min-maf FLOAT		Minimum minor allele frequency for position to be considered heterozygous [0.03]

- Plate map and neighbors (any combination, all optional):
	-m | --map FILE(S)		Optional plate map (tab-separated, no header: sample name, plate position (e.g., A8)); provides substantial speed-up [null]
	-n | --direct BOOL		Compare direct plate neighbors (left, right, top, bottom) [TRUE]
	-d | --diagonal BOOL		Compare diagonal plate neighbors (top-right, bottom-right, top-left, bottom-left) [FALSE]
	-w | --row BOOL			Compare samples in the same row (e.g., row A) [FALSE]
	-l | --column BOOL		Compare samples in the same column (e.g., column 8) [FALSE]
	-t | --plate BOOL		Compare all samples in the same plate map [FALSE]

- Misc:
	-a | --max-mismatches INT	Maximum allowed bases in contaminating sample consensus not matching contaminated sample alleles [1]
	-p | --cores INT		Optional number of cores to use for preprocessing in parallel [1]
	-u | --verbose BOOL		Print progress to STDOUT [TRUE]
	-i | --directory DIRECTORY	Path of directory to store intermediate and temp files [current working directory]
	-o | --output FILE		Output file path [current working directory/potential_cross-contamination.txt]
	-j | --overwrite FILE		Overwrite output, intermediate, and temp files at input paths [FALSE]
```
