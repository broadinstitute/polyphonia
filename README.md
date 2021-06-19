# polyphonia
Detects potential cross-contamination in sequence data.

Usage: `polyphonia cross_contamination [options]`

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
	-m | --plate-map FILE(S)	Optional plate map (tab-separated, no header: sample name, plate position (e.g., A8)); provides substantial speed-up [null]
	-z | --plate-size INT		Standard plate size (6-well, 12-well, 24, 48, 96, 384, 1536, 3456) [96]
	-q | --plate-columns INT	Number columns in plate (e.g., 1, 2, 3, 4) [12]
	-k | --plate-rows INT		Number rows in plate (e.g., A, B, C, D) [8]
	-n | --compare-direct BOOL	Compare direct plate neighbors (left, right, top, bottom) [TRUE]
	-d | --compare-diagonal BOOL	Compare diagonal plate neighbors (top-right, bottom-right, top-left, bottom-left) [FALSE]
	-w | --compare-row BOOL		Compare samples in the same row (e.g., row A) [FALSE]
	-l | --compare-column BOOL	Compare samples in the same column (e.g., column 8) [FALSE]
	-t | --compare-plate BOOL	Compare all samples in the same plate map [FALSE]

- Output:
	-o | --output FILE		Output file path [current working directory/potential_cross-contamination.txt]
	-s | --out-figures DIRECTORY	Path of directory to store plate visualization files [current working directory]
	-x | --out-temp DIRECTORY	Path of directory to store intermediate and temporary files [current working directory]

- Misc:
	-y | --max-mismatches INT	Maximum allowed bases in contaminating sample consensus not matching contaminated sample alleles [1]
	-p | --cores INT		Optional number of cores to use for preprocessing in parallel [1]
	-u | --verbose BOOL		Print progress to STDOUT [TRUE]
	-j | --overwrite BOOL		Overwrite output, intermediate, and temp files at input paths [FALSE]
```

## Contents
- [Getting Started](#getting-started)
- [Under the Hood](#under-the-hood)
- [Required Input Files](#required-input-files)
   - [Reference Genome](#reference-genome)
   - [Consensus Genomes](#consensus-genomes)
   - [Within-Sample Diversity Files](#within-sample-diversity-files)
- [Optional Plate Map Inputs](#optional-plate-map-inputs)
   - [Plate Map File(s)](#plate-map-files)
   - [Plate Map Size](#plate-map-size)
   - [Well Comparison Options](#well-comparison-options)
- [Other Options](#other-options)
   - [Sample Inclusion Thresholds](#sample-inclusion-thresholds)
   - [Allele Filtering Thresholds](#allele-filtering-thresholds)
   - [Cross-Contamination Detection Thresholds](#cross-contamination-detection-thresholds)
   - [Parallelization](#parallelization)
   - [Output File Paths](#output-file-paths)
   - [Verbose](#verbose)
- [Output Files](#output-files)
   - [Potential Cross-Contamination Table](#potential-cross-contamination-table)
   - [Plate Map Visualization of Potential Cross-Contamination](#plate-map-visualization-of-potential-cross-contamination)
   - [Plate Map Visualization of iSNVs](#plate-map-visualization-of-isnvs)
- [Example Run-Throughs](#example-run-throughs)
   - [With VCF Files and Unaligned Consensus Genomes](#with-vcf-files-and-unaligned-consensus-genomes)
   - [With Heterozygosity Tables and Aligned Consensus Genomes](#with-heterozygosity-tables-and-aligned-consensus-genomes)
- [FAQ](#faq)
   - [Dependencies](#dependencies)
   - [Name](#name)


## Getting Started

We use [Docker](https://docs.docker.com/get-started/overview/) to distribute and run polyphonia pre-packaged in an isolated environment with the software it depends on. Using Docker will allow you to run polyphonia on your computer, on the cloud, or elsewhere without worrying about dependencies. If you are new to Docker, you can learn more about it [here](https://docs.docker.com/get-started/).


1. Install Docker. [Download Docker for your system](https://docs.docker.com/get-docker/) and follow instructions to install. When prompted, grant permission to access your file system. Open the Docker application on your computer and leave it open in the background.

2. Update your version of polyphonia. If there have been any updates to polyphonia since you last ran it, your version will be out of date. Update it by entering:

   `docker pull quay.io/broadinstitute/polyphonia:latest`

3. Use Docker to create a new, isolated, ephemeral, pre-built file system containing polyphonia installed alongside the software it requires, with access to your computer's file system. In your terminal application, enter:

   `docker run -v $(pwd):/mnt/data -it --rm quay.io/broadinstitute/polyphonia`

   - [`docker run`](https://docs.docker.com/engine/reference/commandline/run/)` quay.io/broadinstitute/polyphonia` creates a new file system (a container) copied from the snapshot (docker image) we created. (Our snapshot is stored at `quay.io/broadinstitute/polyphonia`.) This new file system is isolated from the rest of your computer.
   - `-v $(pwd):/mnt/data` connects your current directory on your computer (`pwd`) to the `/mnt/data` directory in the new file system. (If you navigate to `/mnt/data` in the new file system, you will see your own familiar files.)
   - `-it` connects and provides a terminal for you to communicate with the new file system.
   - `--rm` automatically destroys the new file system once it is exited.

   You should see your new container appear in the Docker application interface.

4. Explore the space. Regardless of your computer's file system, your new, isolated file system is an Ubuntu file system. Type `ls` to look in a directory and `cd` to navigate to it (or `cd ../` to navigate to the parent directory). Polyphonia lives in the directory `opt/polyphonia`. Your own files are in the directory `/mnt/data`.

5. Run polyphonia. Enter `polyphonia cross_contamination`. You should see a description of the software and a list of options.

   To detect potential cross-contamination in a set of example files included with polyphonia, enter:

   ```
   polyphonia cross_contamination \
   --ref /opt/polyphonia/test/input/NC_045512.2.fasta \
   --vcf /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf \
   --vcf /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf \
   --consensus /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01315-2021.fasta \
   --consensus /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01323-2021.fasta
   ```

   You can read from or write to files outside the isolated file system by replacing the directory you were in when you created your container with `/mnt/data`. On my computer, my files live in a directory called `/Users/lakras/myfiles`. When I created my container I was in `/Users/lakras`—from the perspective of the container, my files are in `/mnt/data/myfiles`. To detect potential cross-contamination in files on my computer, I would enter:

    ```
   polyphonia cross_contamination \
   --ref /mnt/data/myfiles/NC_045512.2.fasta \
   --vcf /mnt/data/myfiles//USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf \
   --vcf /mnt/data/myfiles//USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf \
   --consensus /mnt/data/myfiles/USA-MA-Broad_CRSP-01315-2021.fasta \
   --consensus /mnt/data/myfiles/USA-MA-Broad_CRSP-01323-2021.fasta \
   --output /mnt/data/myfiles/potential_cross_contamination.txt
   ```

   Follow the guide below to learn more about polyphonia's input and output options and to run through this example in more detail.

6. Exit. When you are done, destroy your ephemeral file system by typing `exit`. You should see your container disappear from the Docker interface.

## Under the Hood

TODO

## Required Input Files

### Reference Genome

`--ref`

You must include a reference genome in fasta format. For SARS-CoV-2 samples, we use [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2).

The reference genome input using `--ref` must match that used in other input files. If you provide aligned [consensus genomes](#consensus-genomes) using `--consensus-aligned`, then the reference genome indicated by `--ref` must match the first sequence in the fasta file. If you provide [aligned reads](#within-sample-diversity-files) through `--bam`, then the reads must be aligned to the same reference genome as that indicated by `--ref`. If you provide a pre-processed [within-sample diversity file](#within-sample-diversity-files) through `--vcf` or `--het`, then the positions of loci must be relative to the same reference as that indicated by `--ref`.

### Consensus Genomes

`--consensus`
`--consensus-aligned`

You must include a consensus genome for every sample you want to compare. You can either enter unaligned consensus genomes in one or more fasta files using `--consensus` (in which case they and the [reference](#reference-genome) provided using `--ref` will be aligned using MAFFT), or you can input all consensus genomes aligned to the reference in a single aligned fasta file using `--consensus-aligned`. If you provide aligned consensus genomes using `--consensus-aligned`, then the first sequence in the fasta file must match the [reference genome](#reference-genome) input using `--ref`.

The sample name must occupy the full header line after the `>`. The provided sample name must match the *filename* of a [within-sample diversity file](#within-sample-diversity-files) up to a `.` (if your within-sample diversity file name contains multiple `.`s, polyphonia will try to match all possible names starting with the longest); any sample without a matching within-sample diversity file will be excluded. If you provide a [plate map](#optional-plate-map-inputs) using `--plate-map`, then sample names not appearing in a plate map will be excluded.

Any sample without a consensus genome will be excluded.

### Within-Sample Diversity Files

`--bam`
`--vcf`
`--het`

You must include a within-sample diversity file for every sample you want to compare. A within-sample diversity file can be a bam file with aligned reads provided through `--bam`, a vcf file provided through `--vcf`, or a heterozygosity table provided through `--het`. It is best for all samples to be processed the same way, but polyphonia does not require all within-sample diversity files to be at the same stage.

If you provide aligned reads through `--bam`, then within-sample diversity (base substitutions only) will be called using `LoFreq call` and the output will be processed into a heterozygosity table. These intermediate files are saved in the [directory](#output-file-paths) provided by `out-temp`.

Processing large bam files can be very slow. `--vcf` can be helpful if you have already processed your bam files using LoFreq, or if you prefer to use GATK or another tool. Polyphonia is set up to process vcf files output by LoFreq or GATK. If you provide a within-sample diversity file through `--vcf`, the file will be processed into a heterozygosity table saved in the [directory](#output-file-paths) provided by `out-temp`. You can view example vcf files here: [USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf](/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf) and [USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf](/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf).

If polyphonia cannot read your vcf files, or if you have catalogued within-sample diversity using a tool that produces an output in a different format, you can use `--het` to input a heterozygosity table directly. Each line a heterozygosity table summarizes the minor and major alleles at a locus with heteroyzgosity (base substitutions only). A heterozysity table contains the following columns, tab separated, without a header line:
1. reference genome name (e.g., NC_045512.2)
2. position of locus relative to reference genome, 1-indexed (e.g., 28928)
3. major allele at that position (e.g., C)
4. major allele readcount (e.g., 1026)
5. major allele frequency (e.g., 0.934426)
6. minor allele at that position (e.g., T)
7. minor allele readcount (e.g., 72)
8. minor allele frequency (e.g., 0.065574)

Note that the loci listed in the heterozygosity table will be filtered according to the [allele filtering thresholds](#allele-filtering-thresholds) provided by `--min-readcount` and `--min-maf`. Not all rows will be included in sample comparisons.

You can view example heterozygosity tables here: [USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf_heterozygosity.txt](/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf_heterozygosity.txt) and [USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf_heterozygosity.txt](/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf_heterozygosity.txt).

Any sample without a within-sample diversity file will be excluded.

The sample name associated with a within-sample diversity file is retrieved from the *filename*. Polyphonia will try to match the filename of each within-sample diversity file up to a `.` to a sample name (from the [plate map file(s)](#optional-plate-map-inputs) if provided, or from the [consensus genome files](#consensus-genomes) if not). If a within-sample diversity file name contains multiple `.`s, polyphonia will try to match all possible names starting with the longest. Any within-sample diversity file without an associated [consensus genome](#consensus-genomes) is excluded. If you provide a [plate map](#optional-plate-map-inputs) using `--plate-map`, then any sample whose name does not appear in a plate map is also excluded.

If you provide aligned reads through `--bam`, then the reads must be aligned to the same [reference genome](#reference-genome) as that indicated by `--ref`. If you provide a pre-processed within-sample diversity file through `--vcf` or `--het`, then the positions of loci must be relative to the same [reference](#reference-genome) as that indicated by `--ref`.

## Optional Plate Map Inputs

### Plate Map File(s)

`--plate-map`

If you would like, you can include a plate map file, or multiple plate map files, using `--plate-map`. Your plate map must contain two columns, tab separated, without a header:
1. The name of the sample. Each sample name in a plate map must match a full header line (after the `>`) in the [consensus genome fasta file(s)](#consensus-genomes) and the *name* of a [within-sample diversity file](#within-sample-diversity-files) up to a `.`. (If your within-sample diversity file name contains multiple `.`s, polyphonia will try to match all possible names starting with the longest.) Any sample names that do not have a corresponding consensus genome and a corresponding within-sample diversity file will not be included.
3. The well the sample is in. Each well must be a letter, denoting the row of the well, followed by a number, denoting the column of the well. Column numbers are 1-indexed: the first column is column 1, the second column is column 2, and so on. The letters are A-Z for the first 26 rows, followed by AA-AZ for the next 26 rows, followed by BA-BZ, and so on.

You can view an example of a plate map here: [USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt](/test/input/USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt).

If you include at least one plate map, the plate map(s) will be used to determine what samples should be compared to what other samples based on the [well comparison options](#well-comparison-options) you enter. Any samples not on a plate map will be excluded from comparison, as will any samples that do not have neighbors based on the [well comparison options](#well-comparison-options) entered.

There are two benefits to including a plate map. First, including a plate map allows polyphonia to generate visualizations of [iSNVs](#plate-map-visualization-of-isnvs) and [potenial cross-contamination](#plate-map-visualization-of-potential-cross-contamination) on the plate. Second, including a plate map can provide substantial speed-up, since it limits the number of comparisons that are needed. If you include a plate map, only neighbors are compared, based on the [well comparison options](#well-comparison-options) entered. If you do not include a plate map, all samples are compared to all other samples. This takes [*O(n²)*](https://rob-bell.net/2009/06/a-beginners-guide-to-big-o-notation) time, scaling with *n* samples entered, which gets very slow very fast.

### Plate Map Size

`--plate-size`
`--plate-columns`
`--plate-rows`

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/5/5e/Microplates.jpg/2560px-Microplates.jpg" alt="plate maps" width="500" align="right" />

By default, polyphonia will assume you have a 96-well plate. You can use `--plate-size` to indicate that your plate is any of the following standard plate layouts:
- 6-well plate: 3 columns (A, B, C) x 2 rows (1, 2)
- 12-well plate: 4 columns x 3 rows
- 24-well plate: 6 columns x 4 rows
- 48-well plate: 8 columns x 6 rows
- 96-well plate: 12 columns x 8 rows
- 384-well plate: 24 columns x 16 rows
- 1536-well plate: 48 columns x 32 rows
- 3456-well plate: 72 columns x 48 rows

If you are using a different plate layout, you can indicate it using `--plate-columns` and `--plate-rows`.

For example, if you are using a 384-well plate you can describe its size by entering either `--plate-size 384` or `--plate-columns 24 --plate-rows 16`.

*Image from [Wikipedia](https://en.wikipedia.org/wiki/Microplate#/media/File:Microplates.jpg).*

### Well Comparison Options

`--compare-direct`
`--compare-diagonal`
`--compare-row`
`--compare-column`
`--compare-plate`

If you enter `--compare-direct TRUE`, then each sample will be compared to its direct plate neighbors in the wells to the left, right, top, and bottom.

If you enter `--compare-diagonal TRUE`, then each sample will be compared to its diagonal plate neighbors in the wells to the top-right, bottom-right, top-left, and bottom-left.

If you enter `--compare-row TRUE`, then each sample will be compared to all other samples in the same row.

If you enter `--compare-column TRUE`, then each sample will be compared to all other samples in the same column.

If you enter `--compare-plate TRUE`, then each sample will be compared to all other samples on the plate.

By default, samples are compared only to their direct plate neighbors to the left, right, top, and bottom: `--compare-direct` is set to `TRUE` and all other options are set to `FALSE`.

You can mix and match well comparison options. For example, `--compare-direct TRUE --compare-diagonal TRUE` will compare each sample to the samples in all the wells that surround it: left, right, top, bottom, top-right, bottom-right, top-left, and bottom-left.

## Other Options

### Sample Inclusion Thresholds

`--min-covered`

Use `--min-covered` to set the minimum proportion of the genome that must be covered in order for a sample to be included. The proportion of the genome covered is calculated by counting the number of unambiguous (`A`, `T`, `C`, or `G`)  bases in the sample's consensus genome provided using `--consensus` or `--consensus-aligned`, then dividing by the total number of unambiguous bases in the reference provided using `--ref`.

By default, ≥98% of the genome must be unambigously covered for a sample to be included. If `--min-covered` is set too low, polyphonia may erroneously call potential cross-contamination in or by samples with low genome coverage.

### Allele Filtering Thresholds

`--min-readcount`
`--min-maf`

By default, polyphonia does not include all alleles called by `LoFreq call`. We have found that requiring a minimum minor allele readcount of 10 and a minimum minor allele frequency of 3% after running  `LoFreq call` with default parameters does the best job of separating true within-sample diversity from sequencing errors and other noise. If you would like to specify a different minor allele readcount, you can do that using `--min-readcount`. You can specify a different minimum minor allele frequency using `--min-maf`. These options work regardless of the type(s) of [within-sample diversity files](#within-sample-diversity-files) provided.

### Cross-Contamination Detection Thresholds

`--max-mismatches`

By default, a sample is considered potentially contaminated if at most 1 contaminating consensus-level allele does not appear as a minor or consensus-level allele in the contaminated sample. You can allow for more sequencing errors or missed minor allele calls by specifying a larger threshold using `--max-mismatches`. A too-large `--max-mismatches` will result in many erroneous cross-contamination calls.

(A note that insertions and deletions are ignored: we only include base substitutions.)

### Parallelization

`--cores`

Parallelization is used to substantially speed up processing of [within-sample diversity files](#within-sample-diversity-files) and comparison of neighboring samples.

Use `--cores` to specify how many cores are available. (For example, I use `--cores 4`.) By default, polyphonia will use one core.

### Output File Paths

`--output`
`--out-figures`
`--out-temp`
`--overwrite`

By default, the main output, the [potential cross-contamination table](#potential-cross-contamination-table), are printed to a file named `potential_cross-contamination.txt` in the current working directory within the container, and all output figures, intermediate files, and temporary files are printed to the current working directory within the container.

You can set the output file path for the [potential cross-contamination table](#potential-cross-contamination-table) using `--output`.

If you include a [plate map](#optional-plate-map-inputs), you can set the directory for plate map visualizations of [potential cross-contamination](#plate-map-visualization-of-potential-cross-contamination) and [iSNVs](#plate-map-visualization-of-isnvs) using `--out-figures`.

You can set the output directory for temporary and intermediate files using `--out-temp`.

If an output directory does not already exist, polyphonia will create it.

By default, polyphonia will not overwrite an existing file. You can allow overwriting using `--overwrite TRUE`.

### Verbose

`--verbose`

By default, polyphonia will print updates on its progress to the console. You can silence most updates using `--verbose FALSE`.

## Output Files

### Potential Cross-Contamination Table

`--output`

TODO

You can set the [output file path](#output-file-paths) for the potential cross-contamination table using `--output`.

### Plate Map Visualization of Potential Cross-Contamination

`--out-figures`

TODO

You can set the [output directory](#output-file-paths) for the plate map visualization of potential cross-contamination using `--out-figures`.

### Plate Map Visualization of iSNVs

`--out-figures`

TODO

You can set the [output directory](#output-file-paths) for the plate map visualization of iSNVs using `--out-figures`.

## Example Run-Throughs

You can practice using polyphonia by detecting cross-contamination between two of our SARS-CoV-2 clinical samples: plate neighbors [USA-MA-Broad_CRSP-01315-2021](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13971573) and [USA-MA-Broad_CRSP-01323-2021](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13971565).

We will run polyphonia twice: once with vcf files and unaligned consensus genomes as input, once with heterozygosity tables and aligned consensus genomes as input. (After you are done following along, you can run polyphonia with vcf files and *aligned* consensus genomes and with heterozygosity tables and *unaligned* consensus genomes.)

Before we start, follow [Getting Started](#getting-started) Steps 1-4 to set up your container.

### With VCF Files and Unaligned Consensus Genomes

#### Inputs

Inside your container, the input files we will use for this example are in the directory `/opt/polyphonia/test/input`:

- [**NC_045512.2.fasta**](/test/input/NC_045512.2.fasta): the SARS-CoV-2 reference genome, retrieved from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2).
- [**USA-MA-Broad_CRSP-01315-2021.fasta**](/test/input/USA-MA-Broad_CRSP-01315-2021.fasta) and [**USA-MA-Broad_CRSP-01323-2021.fasta**](/test/input/USA-MA-Broad_CRSP-01323-2021.fasta): the SARS-CoV-2 consensus genome assembled from samples USA-MA-Broad_CRSP-01315-2021 and USA-MA-Broad_CRSP-01323-2021.
- [**USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf**](/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf) and [**USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf**](/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf): within-sample heteroyzgosity in USA-MA-Broad_CRSP-01315-2021 and USA-MA-Broad_CRSP-01323-2021. These vcf files are the result of running `LoFreq call` on USA-MA-Broad_CRSP-01315-2021 and USA-MA-Broad_CRSP-01323-2021 reads aligned to [NC_045512.2.fasta](/test/input/NC_045512.2.fasta) and trimmed. If we used the aligned reads themselves (the very large and therefore not included USA-MA-Broad_CRSP-01315-2021.bam and USA-MA-Broad_CRSP-01323-2021.bam) as input to polyphonia, these vcf files would be generated as intermediate files.
- [**USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt**](/test/input/USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt): a plate map mapping USA-MA-Broad_CRSP-01315-2021 and USA-MA-Broad_CRSP-01323-2021 to their wells in the 96-well plate in which they were processed.

#### Running Polyphonia

To check for potential cross-contamination between USA-MA-Broad_CRSP-01315-2021 and USA-MA-Broad_CRSP-01323-2021, enter:

```
polyphonia cross_contamination \
--ref /opt/polyphonia/test/input/NC_045512.2.fasta \
--vcf /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf \
--vcf /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf \
--consensus /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01315-2021.fasta \
--consensus /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01323-2021.fasta \
--plate-map /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt \
--plate-size 96 \
--output /opt/polyphonia/test/outputs/test1_potential_cross_contamination.txt \
--out-figures /opt/polyphonia/test/test1_plate_visualizations \
--out-temp /opt/polyphonia/test/test1_intermediate_files
```

If it is more convenient, we can alternatively enter vcf and consensus input files as one option each:

```
polyphonia cross_contamination \
--ref /opt/polyphonia/test/input/NC_045512.2.fasta \
--vcf /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf \
--consensus /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01315-2021.fasta /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01323-2021.fasta \
--plate-map /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt \
--plate-size 96 \
--output /opt/polyphonia/test/outputs/test1_potential_cross_contamination.txt \
--out-figures /opt/polyphonia/test/test1_plate_visualizations \
--out-temp /opt/polyphonia/test/test1_intermediate_files
```

(If you run both of these commands in succession, you will run into an error: `Error: file already exists at file path to write to`. Use option `--overwrite TRUE` to allow file overwriting, or `exit` and [regenerate the container](#getting-started) to start over.)

#### Output Files

By default, all intermediate and output files are stored in the current working directory. In this example, we specified output paths and directories for our output files.

Our output table in is in `/opt/polyphonia/test/outputs/test1_potential_cross_contamination.txt`:

<img width="1435" alt="output_table_screenshot" src="https://user-images.githubusercontent.com/6245320/122626233-3becea80-d077-11eb-97af-dccbb3007b1b.png">

The output table shows that USA-MA-Broad_CRSP-01323-2021, in well H10, is potentially contaminated by USA-MA-Broad_CRSP-01315-2021, in well H9. USA-MA-Broad_CRSP-01323-2021 has 21 positions with intrahost variation passing our thresholds. Alleles at all of these positions match the USA-MA-Broad_CRSP-01315-2021 consensus genome—in this case, all matched alleles in positions of heterozygosity are minor alleles, rather than consensus-level alleles. There are no mismatches, or USA-MA-Broad_CRSP-01315-2021 consensus-level alleles not appearing as minor or consensus-level alleles in USA-MA-Broad_CRSP-01323-2021. (A note that insertions and deletions are ignored: we only examine base substitutions.)

The USA-MA-Broad_CRSP-01323-2021 alleles appearing at consensus level in USA-MA-Broad_CRSP-01315-2021 range in frequency from 3.2% to 12.3%, with a median allele frequency of 4.9%. In this case, all potential contaminating alleles appear as minor, rather than consensus-level, alleles in USA-MA-Broad_CRSP-01323-2021.

Our plate visualizations are in the directory `/opt/polyphonia/test/test1_plate_visualizations`:

potential cross-contamination visualization | iSNVs visualization
:-------------------------:|:-------------------------:
![USA-MA-Broad_CRSP-01315_23-2021_plate_map txt_potential_cross_contamination_visualization](https://user-images.githubusercontent.com/6245320/122625873-663da880-d075-11eb-975a-0f7c9ddd19f0.jpg)   |   ![USA-MA-Broad_CRSP-01315_23-2021_plate_map txt_iSNVs_visualization](https://user-images.githubusercontent.com/6245320/122625879-6b9af300-d075-11eb-9b5e-dfdf1308b2aa.jpg)

In the visualization of potential cross-contamination, we see that USA-MA-Broad_CRSP-01323-2021, in well H10, is potentially contaminated by USA-MA-Broad_CRSP-01315-2021, in well H9, and that the median allele frequency of potentially contaminating alleles is 4.9%.

In the iSNVs visualization, we see that USA-MA-Broad_CRSP-01323-2021, in well H10, has 21 positions with intrahost variation (base substitutions only) passing our thresholds, and that USA-MA-Broad_CRSP-01315-2021, in well H9, has 3. We did not examine intrahost variation in any of the other wells.

Our intermediate and temporary files are in the directory `/opt/polyphonia/test/test1_intermediate_files`.

To access these output files on your own system, copy them to `/mnt/data`:

```
cp /opt/polyphonia/test/outputs/test1_potential_cross_contamination.txt /mnt/data/test1_potential_cross_contamination.txt
cp /opt/polyphonia/test/test1_plate_visualizations/* /mnt/data/
```

or run polyphonia with output options set to directories in `/mnt/data`:

```
--output /mnt/data/test1_potential_cross_contamination.txt
--out-figures /mnt/data/test1_plate_visualizations \
--out-temp /mnt/data/test1_intermediate_files
```

### With Heterozygosity Tables and Aligned Consensus Genomes

#### Inputs

As before, the input files we will use for this example are in your directory `/opt/polyphonia/test/input`:

- [**NC_045512.2.fasta**](/test/input/NC_045512.2.fasta): the SARS-CoV-2 reference genome, retrieved from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2).
- [**USA-MA-Broad_CRSP-01315_23-2021_MAFFT_aligned.fasta**](/test/input/USA-MA-Broad_CRSP-01315_23-2021_MAFFT_aligned.fasta): an alignment generated by MAFFT of [NC_045512.2.fasta](/test/input/NC_045512.2.fasta) and the SARS-CoV-2 consensus genomes assembled from USA-MA-Broad_CRSP-01315-2021 and USA-MA-Broad_CRSP-01323-2021 ([USA-MA-Broad_CRSP-01315-2021.fasta](/test/input/USA-MA-Broad_CRSP-01315-2021.fasta) and [USA-MA-Broad_CRSP-01323-2021.fasta](/test/input/USA-MA-Broad_CRSP-01323-2021.fasta)). In our previous example, this alignment was generated as an intermediate file.
- [**USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf_heterozygosity.txt**](/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf_heterozygosity.txt) and [**USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf_heterozygosity.txt**](/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf_heterozygosity.txt): tables summarizing within-sample heteroyzgosity detected by LoFreq in USA-MA-Broad_CRSP-01315-2021 and USA-MA-Broad_CRSP-01323-2021. In our previous example, these tables were generated as intermediate files after processing [USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf](/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf) and [USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf](/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf). (Note that not all rows will pass our filters—most will not be included in the output.)
- [**USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt**](/test/input/USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt): a plate map mapping USA-MA-Broad_CRSP-01315-2021 and USA-MA-Broad_CRSP-01323-2021 to their wells in the 96-well plate in which they were processed.

#### Running Polyphonia

To check for potential cross-contamination between USA-MA-Broad_CRSP-01315-2021 and USA-MA-Broad_CRSP-01323-2021 using these more processed input files, enter:

```
polyphonia cross_contamination \
--ref /opt/polyphonia/test/input/NC_045512.2.fasta \
--het /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf_heterozygosity.txt \
--het /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf_heterozygosity.txt \
--consensus-aligned USA-MA-Broad_CRSP-01315_23-2021_MAFFT_aligned.fasta \
--plate-map /opt/polyphonia/test/input/USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt \
--plate-size 96 \
--output /opt/polyphonia/test/outputs/test2_potential_cross_contamination.txt \
--out-figures /opt/polyphonia/test/test2_plate_visualizations \
--out-temp /opt/polyphonia/test/test2_intermediate_files
```

Like the vcf files in the previous example, we can instead choose to list multiple heterozygosity tables following `--het`. In contrast, we cannot enter multiple aligned fasta files (the single alignment file listed after `--consensus-aligned` must include all of our consensus genomes).

#### Output Files

You should see the same output files as those generated in the previous example.

## FAQ

### Dependencies

If you are using polyphonia through Docker, you do not need to worry about dependencies.

You can view the software polyphonia uses in `requirements-conda.txt`:
- Scripts are written in [Perl 5](https://www.perl.org/get.html).
   - Output float values are rounded using [Math::Round](https://metacpan.org/pod/Math::Round) ≥0.07.
   - Parallelization is performed using [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager) ≥2.02.
- Allele frequencies are calculated by [LoFreq](https://csb5.github.io/lofreq/installation/) ≥2.1.5.
- Sequence alignments are generated by [MAFFT](https://mafft.cbrc.jp/alignment/software/) ≥7.480.
- Visualizations are created in [R](https://www.r-project.org/) ≥4.0.5 using [tidyverse](https://www.tidyverse.org/) ≥1.2.1.

### Name

[**Polyphony**](https://www.youtube.com/watch?v=teh22szdnRQ) describes music containing multiple simultaneous voices with their own melodies, like parts of ["We Will Fall Together" by Streetlight Manifesto](https://www.youtube.com/watch?v=SOqenYis1iQ), ["Somebody That I Used to Know" by Walk off the Earth](https://www.youtube.com/watch?v=d9NF2edxy-M), ["exile" by Taylor Swift with Bon Iver](https://www.youtube.com/watch?v=osdoLjUNFnA), ["Daddy Cool" by Boney M](https://www.youtube.com/watch?v=tHyHOUShOg8), ["Come As You Are" by Nirvana](https://www.youtube.com/watch?v=vabnZ9-ex7o), ["Severed" by The Decemberists](https://www.youtube.com/watch?v=ksTFj6L0mao), and to at least some extent most music created in recent memory.

[**Diafonia**](https://translate.google.com/?sl=it&tl=en&text=diafonia&op=translate) means crosstalk in Italian, Portuguese, and Haitian Creole, and is very close to the words for crosstalk in Spanish (diafonía) and French (diaphonie). In Italian, diafonia also refers to the simultaneous presence of multiple voices, in music and otherwise.

In addition to being the name of our software package, **Polyphonia** is also [a particularly beautiful abstract ballet](https://www.nycballet.com/discover/ballet-repertory/polyphonia/) created by Christopher Wheeldon for New York City Ballet. I think you will find it captures the essence of cross-contamination quite nicely. You can watch an excerpt by clicking the preview image below:

[<img width="1157" alt="ballet_video_preview" src="https://user-images.githubusercontent.com/6245320/121761826-3421c900-cb00-11eb-8712-544ab801fe3c.png">](http://www.youtube.com/watch?v=MOhvllhQo8A)
