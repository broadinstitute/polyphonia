# polyphonia
Detects potential cross-contamination in sequence data.

Usage: `detect_potential_cross_contamination.pl [options]`

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

## Getting Started

We use Docker to distribute and run polyphonia pre-packaged in an isolated, friendly environment with the software it depends on. If you are new to Docker, you can learn more about it [here](https://docs.docker.com/get-started/).


1. Install Docker. [Download Docker for your system](https://docs.docker.com/get-docker/) and follow instructions to install. When prompted, grant permissions to access your file system. Open the Docker application on your computer and leave it open in the background.

2. Use Docker to create a new isolated, ephemeral pre-built file system containing polyphonia installed alongside the software it requires, with access to your computer's file system. In your terminal application, enter:

   `docker run -v $(pwd):/mnt/data -it --rm quay.io/broadinstitute/polyphonia`

   - [`docker run`](https://docs.docker.com/engine/reference/commandline/run/)` quay.io/broadinstitute/polyphonia` creates a new, isolated file system (a container) copied from the snapshot (docker image) we created, which is stored at quay.io/broadinstitute/polyphonia
   - `-v $(pwd):/mnt/data` connects your current directory on your computer (`pwd`) to the `/mnt/data` directory in the new file system. (If you navigate to `/mnt/data` in the new file system, you will see your own familiar files.)
   - `-it` connects and provides a terminal for you to communicate with the new file system
   - `--rm` automatically destroys the file system once it is exited

   You should see your new container appear in the Docker application interface.

3. Explore the space. Regardless of your computer's file system, your new, isolated file system is an Ubuntu file system. Type `ls` to look in a directory and `cd` to navigate to it (or `cd ../` to navigate to the parent directory). Polyphonia lives in the directory `opt/polyphonia`. Your own files are in the directory `/mnt/data/`.

4. Run polyphonia. Enter `/opt/polyphonia/detect_potential_cross_contamination.pl`. You should see a description of the software and a list of options. To detect potential cross-contamination in a set of example files included with polyphonia, enter:

   `TODO`

   You can detect potential cross-contamination in your own files by adding `/mnt/data` to the start of your current directory when you created your container. On my computer, my files live in a directory called `/Users/lakras/myfiles`. When I created my container I was in `/Users/lakras`—from the perspective of the container, my files are in `/mnt/data/myfiles`. To detect potential cross-contamination in files on my computer, I would enter:

   `TODO`

   Follow the guide below to learn more about polyphonia's input and output options.

5. Exit. When you are done, destroy your ephemeral file system by typing `exit`. You should see your container disappear from the Docker interface.

## Required Input Files

TODO

## Optional Plate Map Inputs

TODO

## Other Options

TODO

## Output Files

TODO

