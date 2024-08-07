# polyphonia
Detects potential cross-contamination in sequence data.

Usage: `polyphonia cross_contamination [options]`

```
OPTIONS:
- Reference (required):
	-f | --ref FILE			Reference fasta file [null]

- Consensus genomes (aligned or not aligned, not both; at least one file required):
	-c | --consensus FILE(S)	Unaligned consensus genome or genomes [null]
	-a | --consensus-aligned FILE	Fasta alignment with consensus genomes pre-aligned to reference; reference provided by --ref must appear first [null]

- Within-sample diversity (any combination; at least one file required):
	-b | --bam FILE(S)		Aligned and trimmed reads as bam file(s); must use reference provided by --ref [null]
	-v | --vcf FILE(S)		VCF file(s) output by LoFreq; must use reference provided by --ref [null]
	-h | --het FILE(S)		Tab-separated heterozygosity summary tables; see documentation for format [null]

- Filtering options:
	-i | --min-maf FLOAT		Minimum minor allele frequency for a position to be considered heterozygous [0]
	-e | --min-readcount INT	Minimum minor allele readcount for a position to be considered heterozygous [0]
	-r | --min-depth INT		Minimum read depth for a position to be used for comparison [100]
	-1 | --read-depths FILE(S)	Read depth tables; provide alongside vcf files or heterozygosity tables if min-depth>0; see documentation for format [null]
	-g | --min-covered FLOAT	Minimum proportion genome that must be covered at minimum read depth for a sample to be included [0.95]
	-3 | --masked-positions STRING	1-indexed positions to mask (e.g., 1-10,50,55-70) [null]
	-4 | --masked-positions-file FILE	1-indexed positions to mask, one per line [null]
	-y | --max-mismatches INT	In flagged potential cross-contamination, maximum allowed unambiguous bases in contaminating sample consensus not matching contaminated sample alleles [0]
	-5 | --min-matches INT		Of positions at which the two consensus genomes differ, the minimum number of positions at which contamination is detected as a minor allele [3]
	-6 | --min-matches-proportion FLOAT	Of positions at which the two consensus genomes differ, the minimum proportion of positions at which contamination is detected as a minor allele [1]

- Plate map and neighbors (any combination, all optional):
	-m | --plate-map FILE(S)	Optional plate map(s) (tab-separated, no header: sample name, plate position (e.g., A8)); provides substantial speed-up [null]
	-z | --plate-size INT		Standard plate size (6-well, 12-well, 24, 48, 96, 384, 1536, or 3456) [96]
	-q | --plate-columns INT	Number columns in plate (e.g., 1, 2, 3, 4) [12]
	-k | --plate-rows INT		Number rows in plate (e.g., A, B, C, D) [8]
	-n | --compare-direct BOOL	Compare direct plate neighbors (left, right, top, bottom) [TRUE]
	-d | --compare-diagonal BOOL	Compare diagonal plate neighbors (top-right, bottom-right, top-left, bottom-left) [FALSE]
	-w | --compare-row BOOL		Compare samples in the same row (e.g., row A) [FALSE]
	-l | --compare-column BOOL	Compare samples in the same column (e.g., column 8) [FALSE]
	-t | --compare-plate BOOL	Compare all samples in each plate map [FALSE]

- Output:
	-o | --output FILE		Output file path [current working directory/potential_cross-contamination.txt]
	-s | --out-figures DIRECTORY	Path of directory to store plate visualization files [current working directory]
	-x | --out-temp DIRECTORY	Path of directory to store intermediate and temporary files [current working directory]

- Misc:
	-p | --cores INT		Optional number of cores to use for preprocessing in parallel [1]
	-u | --verbose BOOL		Print progress updates to STDOUT [TRUE]
	-j | --overwrite BOOL		Overwrite files that already exist at output, intermediate, and temp file paths [FALSE]
	-2 | --print-all-iSNVs BOOL	Include all threshold-passing samples in iSNVs visualizations, including samples without plate neighbors [FALSE]
	-0 | --print-all BOOL		Output outcomes of all comparisons (all comparisons are marked as potential cross-contamination) [FALSE]
```

## Contents
- [Getting Started](#getting-started)
   - [Docker](#docker)
   - [Bioconda](#bioconda)
   - [Dockstore](#dockstore)
- [Under the Hood](#under-the-hood)
- [Important Caveats](#important-caveats)
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
   - [Position Inclusion Thresholds](#position-inclusion-thresholds)
   - [Allele Filtering Thresholds](#allele-filtering-thresholds)
   - [Cross-Contamination Detection Thresholds](#cross-contamination-detection-thresholds)
   - [Parallelization](#parallelization)
   - [Output File Paths](#output-file-paths)
   - [Verbose](#verbose)
   - [Print All](#print-all)
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

### Docker

We use [Docker](https://docs.docker.com/get-started/overview/) to distribute and run polyphonia pre-packaged in an isolated environment with the software it depends on. Using Docker will allow you to run polyphonia on your computer, on the cloud, or elsewhere without worrying about dependencies. If you are new to Docker, you can learn more about it [here](https://docs.docker.com/get-started/).


1. Install Docker. [Download Docker for your system](https://docs.docker.com/get-docker/) and follow instructions to install. When prompted, grant permission to access your file system. Open the Docker application on your computer and leave it open in the background.

2. Update your version of polyphonia. If there have been any updates to polyphonia since you last ran it, your version will be out of date. Update it by entering:

   `docker pull quay.io/broadinstitute/polyphonia:latest`

3. Use Docker to create a new, isolated, ephemeral, pre-built file system containing polyphonia installed alongside the software it requires and with access to your computer's file system. In your terminal application, enter:

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

6. Exit. When you are done, destroy your ephemeral file system by typing `exit`. You should see your container disappear from the Docker interface. Make sure to also quit the Docker application when you are done using it.

### Bioconda

While Docker runs polyphonia together with dependencies in an isolated environment, bioconda solves the same problem by instead installing all [dependencies](#dependencies) onto your machine.

Polyphonia will be available as a [bioconda](https://bioconda.github.io/) package soon.

### Dockstore

Polyphonia is available through the following [WDL](https://github.com/openwdl/wdl) workflows and tasks included in [broadinstitute/viral-pipelines](https://github.com/broadinstitute/viral-pipelines):
- Workflows:
    - [`detect_cross_contamination`](https://github.com/broadinstitute/viral-pipelines/blob/master/pipes/WDL/workflows/detect_cross_contamination.wdl) is available [on Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/viral-pipelines/detect_cross_contamination:master?tab=info) for execution on [Terra](https://terra.bio/), [DNAnexus](https://www.dnanexus.com/), and other WDL-compatible platforms. The `detect_cross_contamination` workflow runs polyphonia with [unaligned consensus genomes](#consensus-genomes) and with within-sample diversity captured in [aligned reads in bam files](#--bam). In contrast to polyphonia run on [Docker](#docker), this workflow runs [`LoFreq call`](https://csb5.github.io/lofreq/commands/#call) on all samples, including samples that do not pass [sample inclusion criteria](#sample-inclusion-thresholds) or do not have neighbors according to the [well comparison options](#well-comparison-options): these filters are applied after rather than before running `LoFreq call`.
    - [`detect_cross_contamination_precalled_vcfs`](https://github.com/broadinstitute/viral-pipelines/blob/master/pipes/WDL/workflows/detect_cross_contamination_precalled_vcfs.wdl) is also available [on Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/viral-pipelines/detect_cross_contamination_precalled_vcfs:master?tab=info). The `detect_cross_contamination_precalled_vcfs` workflow runs polyphonia with [unaligned consensus genomes](#consensus-genomes) and with within-sample diversity captured in [vcf files](#--vcf) generated by LoFreq.
- Tasks:
    [`detect_cross_contamination`](https://github.com/broadinstitute/viral-pipelines/blob/master/pipes/WDL/tasks/tasks_intrahost.wdl).

## Under the Hood
![for_github](https://github.com/broadinstitute/polyphonia/assets/6245320/e04a8bd2-0eda-4111-b5fd-059c1b82dd9a)

Polyphonia starts off by verifying and printing input options and preparing a list of samples to analyze. Samples without a [consensus genome](#consensus-genomes) or [within-sample diversity file](#within-sample-diversity-files) are excluded from analysis. If at least one [optional plate map](#optional-plate-map-inputs) is provided, samples not appearing in any plate map are excluded from analysis. A [read depth filter](#position-inclusion-thresholds) is applied. Samples not passing [sample inclusion thresholds](#sample-inclusion-thresholds) are excluded from analysis. To save time, samples without any plate neighbors as specified by provided [well comparison options](#well-comparison-options) are excluded.

If [consensus genomes](#consensus-genomes) are not already aligned, they are aligned using [`MAFFT`](https://mafft.cbrc.jp/alignment/software/).

[Within-sample diversity files](#within-sample-diversity-files) are pre-processed depending on their stage in processing, [in parallel](#parallelization) if possible. If aligned reads are provided in a bam file, they are processed using [`LoFreq call`](https://csb5.github.io/lofreq/commands/#call) into a vcf file, which is in turn processed into a heterozygosity table cataloguing within-sample diversity (base substitutions only). If the aligned reads bam file is large, processing it can take a long time. If a vcf file is provided, it is processed into a heterozygosity table. (If a heterozygosity table is provided, it does not need to be processed.) [Allele filtering thresholds](#allele-filtering-thresholds) are applied to the heterozygosity tables as they are read in.

Pairs of samples are then compared to detect potential cross-contamination. A sample B is marked as potentially contaminated by another sample A if at positions where the consensus genomes of A and B differ, A's consensus alleles appear as minor alleles in B, with a minimum number and proportion of genome-defining positions and allowing for a number of mismatches specified by input [cross-contamination detection thresholds](#cross-contamination-detection-thresholds). Only positions with unambiguous bases (`A`, `T`, `C`, or `G`) in both samples are compared.

If at least one plate map is provided, each sample is compared to the samples in neighboring wells as determined by provided [well comparison options](#well-comparison-options). If no plate map is provided, all samples are compared to all other samples. (As the number of samples increases, comparing all samples to all other samples very quickly becomes a large and intractably slow task.) If possible, comparisons are made [in parallel](#parallelization).

The primary output is [a table describing potential instances of cross-contamination](#potential-cross-contamination-table). If at least one plate map is provided, visualizations of [potential cross-contamination](#plate-map-visualization-of-potential-cross-contamination) and of [intrahost variation](#plate-map-visualization-of-isnvs) are also produced for each plate.

The rest is up to you. Depending on availability of material and the stage at which contamination occurred, contaminated samples could be resequenced or subjected to other follow-up. It is important to critically examine flagged potential cross-contamination and consider [alternative explanations](#important-caveats) for shared alleles: not all flagged sample pairs will be true instances of cross-contamination or require resequencing.

## Important Caveats

Not all potential cross-contamination flagged by polyphonia will be true cross-contamination. You should examine and verify any potential cross-contamination that polyphonia flags. Potential cross-contamination at the level of minor alleles can also be explained by the "contaminated" sample containing multiple infections. In the example we run through from the Broad, we suspect that USA-MA-Broad_CRSP-01323-2021 is contaminated by USA-MA-Broad_CRSP-01315-2021 because they are plate neighbors; however, it is also possible that USA-MA-Broad_CRSP-01323-2021 contains a second, minor infection identical to the consensus-level infection in USA-MA-Broad_CRSP-01315-2021. There is no way for us to know for sure.

On the other hand, there are situations in which polyphonia will not detect cross-contamination. Polyphonia does not detect consensus-level or near-consensus-level potential cross-contamination (with contaminating allele frequency at, above, or near 50%). Polyphonia will not detect cross-contamination if the contaminated and contaminating samples have identical genomes, if the contaminating material does not cover the entire genome (specifically if it does not cover the loci at which the contaminating and consensus-level genomes differ), if either the contaminated or contaminating consensus sequence covers less of the genome than specified by [`--min-covered`](#sample-inclusion-thresholds), or if the contaminating material appears at very low readcount or allele frequency below the [allele filtering thresholds](#allele-filtering-thresholds) set by `--min-readcount` and `--min-maf`.

## Required Input Files

### Reference Genome

`--ref`

You must include a reference genome in fasta format. For SARS-CoV-2 samples, we use [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2).

The reference genome that is input using `--ref` must match that used in other input files. If you provide aligned [consensus genomes](#consensus-genomes) using `--consensus-aligned`, then the first sequence in the alignment must match the reference genome indicated by `--ref`. If you provide [aligned reads](#within-sample-diversity-files) through `--bam`, then the reads must be aligned to the same reference genome as that indicated by `--ref`. If you provide a pre-processed [within-sample diversity file](#within-sample-diversity-files) through `--vcf` or `--het` or a [read depth table](#position-inclusion-thresholds) through `--read-depths`, then the positions of loci must be relative to the same reference as that indicated by `--ref`.

### Consensus Genomes

`--consensus`
`--consensus-aligned`

You must include a consensus genome for every sample you want to compare. Any sample without a consensus genome will be excluded.

You can either enter unaligned consensus genomes in one or more fasta files using `--consensus` (in which case they and the [reference](#reference-genome) provided using `--ref` will be aligned using [`MAFFT`](https://mafft.cbrc.jp/alignment/software/)), or you can input all consensus genomes aligned to the reference in a single aligned fasta file using `--consensus-aligned`. If you provide aligned consensus genomes using `--consensus-aligned`, then the first sequence in the alignment must match the [reference genome](#reference-genome) that is input using `--ref`.

The sample name must occupy the full header line after the `>`. The provided sample name must match the *filename* of the corresponding [within-sample diversity file](#within-sample-diversity-files) up to a `.` (if your within-sample diversity file name contains multiple `.`s, polyphonia will try to match all possible names starting with the longest); any sample without a matching within-sample diversity file will be excluded. If you provide at least one [plate map](#optional-plate-map-inputs) using `--plate-map`, then sample names not appearing in any plate map will be excluded.

### Within-Sample Diversity Files

`--bam`
`--vcf`
`--het`

You must include a within-sample diversity file for every sample you want to compare. Any sample without a within-sample diversity file will be excluded.

A within-sample diversity file can be a bam file with aligned reads provided through `--bam`, a vcf file provided through `--vcf`, or a heterozygosity table provided through `--het`. It is best for all samples to be processed the same way, but polyphonia does not require all within-sample diversity files to be at the same stage.

The sample name associated with a within-sample diversity file is retrieved from the *filename*. Polyphonia will try to match the filename of each within-sample diversity file up to a `.` to an existing sample name (from the [plate map file(s)](#optional-plate-map-inputs) if provided, or from the [consensus genome files](#consensus-genomes) if not). If the name of a within-sample diversity file contains multiple `.`s, polyphonia will try to match all possible names, starting with the longest, until it finds a matching sample name. Any within-sample diversity file without an associated [consensus genome](#consensus-genomes) is excluded. If you provide a [plate map](#optional-plate-map-inputs) using `--plate-map`, any within-sample diversity file without an associated plate map entry is excluded.

#### `--bam`
If you provide aligned reads through `--bam`, then within-sample diversity (base substitutions only) will be called using [`LoFreq call`](https://csb5.github.io/lofreq/commands/#call) and the output vcf file will be processed into a heterozygosity table. These intermediate files are saved in the [directory](#output-file-paths) provided by `out-temp`. Within the bam file, reads must be aligned to the same [reference genome](#reference-genome) as that provided by `--ref`.

#### `--vcf`
Processing large bam files can be very slow. `--vcf` can be helpful if you have already processed your bam files using [`LoFreq call`](https://csb5.github.io/lofreq/commands/#call). If you provide a within-sample diversity file through `--vcf`, the file will be processed into a heterozygosity table saved in the [directory](#output-file-paths) provided by `out-temp`. You can view example vcf files here: [USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf](/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf) and [USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf](/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf).

Within the vcf file, locus positions must be relative to the same [reference](#reference-genome) as that provided by `--ref`.

If a non-zero minimum readcount is set by [`--min-depth`](#position-inclusion-thresholds), then each vcf file must be accompanied by a read depth table provided by [`--read-depths`](#position-inclusion-thresholds). If within-sample diversity files are provided as vcf files or heterozygosity tables without any accompanying read depth tables, then [`--min-depth`](#position-inclusion-thresholds) is set to 0. If [`--min-depth`](#position-inclusion-thresholds) is non-zero and at least one read depth table is provided, then any vcf files without corresponding read depth tables are excluded.

Polyphonia is set up to process vcf files that were output by LoFreq. If you use a different tool, you can preprocess the output yourself and input it using `--het`.

#### `--het`
If polyphonia cannot read your vcf files, or if you have catalogued within-sample diversity using a tool that produces an output in a different format, you can pre-process the output files yourself and use `--het` to input a heterozygosity table directly. Each line in a heterozygosity table summarizes the major and minor alleles at a locus with heteroyzgosity (base substitutions only). In a heterozysity table, every row lists a position in the genome relative to the reference genome. A heterozysity table contains the following tab-separated columns without a header line:
1. name of reference genome (e.g., NC_045512.2)
2. position of locus relative to reference genome, 1-indexed (e.g., 28928)
3. major allele at that position (e.g., C)
4. major allele readcount (e.g., 1026)
5. major allele frequency (e.g., 0.934426)
6. minor allele at that position (e.g., T)
7. minor allele readcount (e.g., 72)
8. minor allele frequency (e.g., 0.065574)

Locus positions must be relative to the same [reference](#reference-genome) as that provided by `--ref`.

Note that the loci listed in the heterozygosity table will be filtered according to the [allele filtering thresholds](#allele-filtering-thresholds) provided by `--min-readcount` and `--min-maf`. Not all alleles will be included in sample comparisons.

If a non-zero minimum readcount is set by [`--min-depth`](#position-inclusion-thresholds), then each heterozygosity table must be accompanied by a read depth table provided by [`--read-depths`](#position-inclusion-thresholds). If within-sample diversity files are provided as heterozygosity tables or vcf files without any accompanying read depth tables, then [`--min-depth`](#position-inclusion-thresholds) is set to 0. If [`--min-depth`](#position-inclusion-thresholds) is non-zero and at least one read depth table is provided, then any heterozygosity tables without corresponding read depth tables are excluded.

You can view example heterozygosity tables here: [USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf_heterozygosity.txt](/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf_heterozygosity.txt) and [USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf_heterozygosity.txt](/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf_heterozygosity.txt).

## Optional Plate Map Inputs

### Plate Map File(s)

`--plate-map`

If you would like, you can include a plate map file, or multiple plate map files, using `--plate-map`. Your plate map must contain two tab-separated columns without a header line:
1. The name of the sample. Each sample name in a plate map must match a full header line (after the `>`) in the [consensus genome fasta file(s)](#consensus-genomes) and the *name* of a [within-sample diversity file](#within-sample-diversity-files) up to a `.`. (If your within-sample diversity file name contains multiple `.`s, polyphonia will try to match all possible names starting with the longest.) Any sample names that do not have a corresponding consensus genome and a corresponding within-sample diversity file will not be included.
3. The well the sample is in. Each well must be a letter, denoting the row of the well, followed by a number, denoting the column of the well. Column numbers are 1-indexed: the first column is column 1, the second column is column 2, and so on. Row letters are A-Z for the first 26 rows, followed by AA-AZ for the next 26 rows, followed by BA-BZ, and so on.

You can view an example of a plate map here: [USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt](/test/input/USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt).

If you include at least one plate map, the plate map(s) will be used to determine what samples should be compared to what other samples based on the [well comparison options](#well-comparison-options) you enter. If you include at least one plate map, any samples not on a plate map will be excluded from comparison, as will any samples that do not have neighbors based on the [well comparison options](#well-comparison-options) entered.

There are two benefits to including a plate map. First, including a plate map allows polyphonia to generate visualizations of [iSNVs](#plate-map-visualization-of-isnvs) and [potenial cross-contamination](#plate-map-visualization-of-potential-cross-contamination) on the plate. Second, including a plate map can provide substantial speed-up, since it limits the number of comparisons that are needed. If you include a plate map, only neighbors are compared, based on the [well comparison options](#well-comparison-options) entered. If you do not include a plate map, all samples are compared to all other samples. This takes [*O(n²)*](https://rob-bell.net/2009/06/a-beginners-guide-to-big-o-notation) time, scaling with *n* samples entered, which gets very slow very fast.

### Plate Map Size

`--plate-size`
`--plate-columns`
`--plate-rows`

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/5/5e/Microplates.jpg/2560px-Microplates.jpg" alt="plate maps" width="500" align="right" />

By default, polyphonia will assume you have a 96-well plate. You can use `--plate-size` to indicate that your plate is any of the following standard plate layouts:
- 6-well plate: 3 columns (1, 2, 3) x 2 rows (A, B)
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

If you enter `--compare-direct TRUE`, then each sample is compared to its direct plate neighbors in the wells to the left, right, top, and bottom. This option is on (`TRUE`) by default and unless otherwise specified using `--compare-direct FALSE`.

<p align="center"><img src="https://user-images.githubusercontent.com/6245320/122707970-17387480-d229-11eb-9da0-4f84909341ff.png" alt="compare-direct" width=400></p>

If you enter `--compare-diagonal TRUE`, then each sample is compared to its diagonal plate neighbors in the wells to the top-right, bottom-right, top-left, and bottom-left. This option is off (`FALSE`) by default.

<p align="center"><img src="https://user-images.githubusercontent.com/6245320/122708126-70080d00-d229-11eb-838f-635b8cde4609.png" alt="compare-diagonal" width=400></p>

If you enter `--compare-row TRUE`, then each sample is compared to all other samples in the same row. This option is off (`FALSE`) by default.

<p align="center"><img src="https://user-images.githubusercontent.com/6245320/122708147-79917500-d229-11eb-98db-8d52c91d33ea.png" alt="compare-row " width=400></p>

If you enter `--compare-column TRUE`, then each sample is compared to all other samples in the same column. This option is off (`FALSE`) by default.

<p align="center"><img src="https://user-images.githubusercontent.com/6245320/122708153-7dbd9280-d229-11eb-8845-c4ccd703dbe9.png" alt="compare-column" width=400></p>

If you enter `--compare-plate TRUE`, then each sample is compared to all other samples on the plate. This option is off (`FALSE`) by default.

<p align="center"><img src="https://user-images.githubusercontent.com/6245320/122708159-81511980-d229-11eb-83ad-144cb6191441.png" alt="compare-plate" width=400></p>

You can mix and match well comparison options. For example, `--compare-direct TRUE --compare-diagonal TRUE` compares each sample to the samples in all the wells that surround it: left, right, top, bottom, top-right, bottom-right, top-left, and bottom-left.

<p align="center"><img src="https://user-images.githubusercontent.com/6245320/122708120-6a122c00-d229-11eb-8b98-e4b845a48b73.png" alt="compare-direct and compare-diagonal" width=400></p>

By default, samples are compared only to their direct plate neighbors to the left, right, top, and bottom: `--compare-direct` is set to `TRUE` and all other options are set to `FALSE`. Note that even if another option is set to `TRUE`, `--compare-direct` is only set to `FALSE` if `--compare-direct FALSE` is explicitly specified.

## Other Options

### Sample Inclusion Thresholds

`--min-covered`

Use `--min-covered` to set the minimum proportion of the genome that must be covered in order for a sample to be included. The proportion of the genome covered is calculated by counting the number of unambiguous (`A`, `T`, `C`, or `G`)  bases in the sample's consensus genome provided using `--consensus` or `--consensus-aligned`, then dividing by the total number of unambiguous bases in the reference provided using `--ref`. If [`--min-depth`](#position-inclusion-thresholds) is non-zero, then only positions passing the read depth filter are included.

By default, ≥95% of the genome must be unambigously covered for a sample to be included. If `--min-covered` is set too low, polyphonia may erroneously call potential cross-contamination in or by samples with low genome coverage.

### Position Inclusion Thresholds

`--min-depth`
`--read-depths`
`--masked-positions`
`--masked-positions-file`

Use `--masked-positions` or `--masked-positions-file` to specify positions that should be masked, for example positions known to have high rates of sequencing error. Masked positions are not included as input consensus or minor alleles.

Use `--min-depth` to set the minimum number of reads that must overlap a position in a sample in order for that position to be included. If a position does not pass the read depth filter, it is not included in the numerator of [genome coverage](#sample-inclusion-thresholds) and is not included as a heterozygous position.

If [aligned reads](#--bam) are provided as input [within-sample diversity files](#within-sample-diversity-files), read depths are calculated from aligned reads using [`SAMtools depth`](http://www.htslib.org/doc/samtools-depth.html).

If pre-processed [within-sample diversity file](#within-sample-diversity-files) are provided through `--vcf` or `--het`, a read depth table corresponding to each pre-processed within-sample diversity file can be input using the `--read-depths` option.

The read depth table should be identical to that generated by running [`SAMtools depth`](http://www.htslib.org/doc/samtools-depth.html) on the [aligned reads](#--bam). It must contain three tab-separated columns without a header line:
1. name of reference genome (e.g., NC_045512.2)
2. position of locus relative to reference genome, 1-indexed (e.g., 28928)
3. read depth at that position (e.g., 4710)

Within the read depth table, locus positions must be relative to the same [reference](#reference-genome) as that provided by `--ref`.

By default, minimum read depth is set to 100 reads, or 0 reads if no read depth tables and at least one non-bam [within-sample diversity file](#within-sample-diversity-files) are provided as input. If a minimum read depth is non-zero, any sample with a non-bam [within-sample diversity file](#within-sample-diversity-files) but no read depth file is excluded from analysis.

If `--min-depth` is set too low, polyphonia may erroneously call mismatches at positions where read depth is too low for a minor allele to appear, and non-consensus-level cross-contamination may be missed.

### Allele Filtering Thresholds

`--min-readcount`
`--min-maf`

By default, polyphonia includes all alleles called by `LoFreq call`. You can set a minimum minor allele readcount and a minimum minor allele frequency after running `LoFreq call` to help separate true within-sample diversity from sequencing errors and other noise. If you would like to specify a minimum minor allele readcount, you can do that using `--min-readcount`. You can specify a minimum minor allele frequency using `--min-maf`. These options work regardless of the type(s) of [within-sample diversity files](#within-sample-diversity-files) provided.

If [`--min-depth`](#position-inclusion-thresholds) is non-zero, then an allele will be excluded if it appears at a position that does not pass the read depth threshold—even if the minor allele frequency and readcount pass the thresholds set by `--min-readcount` and `--min-maf`.

### Cross-Contamination Detection Thresholds

`--min-matches`
`--min-matches-proportion`
`--max-mismatches`

A sample is considered potentially contaminated if a contaminating sample's consensus alleles appear as minor alleles in positions where the two samples' consensus genomes differ. Of positions where the consensus genomes of the two samples differ, `--min-matches` specifies the minimum number (by default 3) and `--min-matches-proportion` specifies the minimum proportion (by default 100%) of positions where contaminating alleles must appear as minor alleles for a cross-contamination event to be called.

You can allow for more sequencing errors or missed minor allele calls by specifying a larger threshold using `--max-mismatches`, by default 0. A too-large `--max-mismatches` will result in many erroneous cross-contamination calls.

Alleles are only compared or, therefore, counted as matches or mismatches if they pass the [allele filtering thresholds](#allele-filtering-thresholds) and appear at positions that pass the [read depth filter](#position-inclusion-thresholds).

(A note that insertions and deletions are ignored: we only include base substitutions.)

### Parallelization

`--cores`

Parallelization is used to substantially speed up processing of [within-sample diversity files](#within-sample-diversity-files) and comparison of neighboring samples.

Use `--cores` to specify how many cores are available for polyphonia to use. (For example, I use `--cores 4`.) By default, polyphonia will use one core.

### Output File Paths

`--output`
`--out-figures`
`--out-temp`
`--overwrite`

By default, the main output file, the [potential cross-contamination table](#potential-cross-contamination-table), is printed to a file named `potential_cross-contamination.txt` in the current working directory within the container, and all output figures, intermediate files, and temporary files are printed to the current working directory within the container.

You can set the output file path for the [potential cross-contamination table](#potential-cross-contamination-table) using `--output`.

If you include a [plate map](#optional-plate-map-inputs), you can set the directory for plate map visualizations of [potential cross-contamination](#plate-map-visualization-of-potential-cross-contamination) and [iSNVs](#plate-map-visualization-of-isnvs) using `--out-figures`.

You can set the output directory for temporary and intermediate files using `--out-temp`.

If an output directory does not already exist, polyphonia will create it.

By default, polyphonia will not overwrite an existing file. You can allow overwriting using `--overwrite TRUE`.

### Verbose

`--verbose`

By default, polyphonia will print updates on its progress to the console. You can silence most updates using `--verbose FALSE`.

### Print All

`--print-all`

Use `--print-all TRUE` for all comparisons to be output as potential cross-contamination. In other words, the output potential [cross-contamination table](#potential-cross-contamination-table) will instead list all comparisons. Similarly, the [plate map visualization of potential cross-contamination](#plate-map-visualization-of-potential-cross-contamination) will indicate that every sample is potentially contaminated by every sample it was compared to. `--print-all` is useful for debugging or otherwise peeking under the hood; it is not meant for everyday use.

## Output Files

### Potential Cross-Contamination Table

`--output`

The primary output file generated by polyphonia is a table detailing evidence for potential cross-contamination events. You can set the [output file path](#output-file-paths) for the potential cross-contamination table using `--output`.

You can view an example of a potential cross-contamination table at the [end](#output-files-1) of the [example run-through](#example-run-throughs):

<img width="1749" alt="output_table" src="https://github.com/broadinstitute/polyphonia/assets/6245320/cb2180f2-dc54-45d3-ab1d-d07bbef5249f">

Each row in the output table corresponds to a single potential cross-contamination event. The output table has the following columns—

Columns summarizing the potentially contaminated sample:
- `potential_contaminated_sample`: the sample proposed to potentially be contaminated. (E.g., `USA-MA-Broad_CRSP-01323-2021`.)
- `potential_contaminated_sample_unambiguous_bases`: the number of unambiguous (`A`, `T`, `C`, or `G`) bases in the potentially contaminated samples's consensus genome. (E.g., `29,830`.)
- `potential_contaminated_sample_genome_covered`: the number of unambiguous (`A`, `T`, `C`, or `G`) bases in the potentially contaminated samples's consensus genome divided by the number of unambiguous bases in the reference genome. The proportion is presented as a percentage rounded to one decimal place. (E.g., `99.8%`.)
- `potential_contaminated_sample_unambiguous_bases_passing_read_depth_filter`: the number of unambiguous (`A`, `T`, `C`, or `G`) bases in the potentially contaminated samples's consensus genome that pass the [read depth filter](#position-inclusion-thresholds). This column only appears if `--min-depth` is non-zero. (E.g., `29,728`.)
- `potential_contaminated_sample_genome_covered_passing_read_depth_filter`: the number of unambiguous (`A`, `T`, `C`, or `G`) bases in the potentially contaminated samples's consensus genome that pass the [read depth filter](#position-inclusion-thresholds) divided by the number of unambiguous bases in the reference genome. The proportion is presented as a percentage rounded to one decimal place. This column only appears if `--min-depth` is non-zero. (E.g., `99.4%`.)
- `num_positions_with_heterozygosity`: the number of loci in the potentially contaminated sample at which there is heterozygosity (base substitutions only) passing [allele filtering thresholds](#allele-filtering-thresholds). (E.g., `21`.)
- `alleles_at_positions_with_heterozygosity`: the alleles at each heterozygous locus (base substitutions only) in the potentially contaminated sample. Each locus is presented as its position, the consensus-level allele, and the minor allele. Loci are separated by a semicolon and a space (`; `). (E.g., `2,162 AG; 2,813 AG; 3,044 CA; 5,452 GA; 6,429 TC; 6,762 CT; 7,348 TA; 9,391 TC; 10,702 TC; 13,119 CT; 14,484 CT; 18,131 CT; 19,072 GT; 19,868 CT; 21,203 AG; 24,703 TC; 27,630 TC; 29,095 TC; 29,272 TC; 29,360 TC; 29,367 TC`.)

Columns summarizing the potentially contaminating sample:
- `potential_contaminating_sample`: the sample proposed to be the source of this contamination event. (E.g., `USA-MA-Broad_CRSP-01315-2021`.)
- `potential_contaminating_sample_unambiguous_bases`: the number of unambiguous (`A`, `T`, `C`, or `G`) bases in the potentially contaminating samples's consensus genome. (E.g., `29,782`.)
- `potential_contaminating_sample_genome_covered`: the number of unambiguous (`A`, `T`, `C`, or `G`) bases in the potentially contaminating samples's consensus genome divided by the number of unambiguous bases in the reference genome. The propotion is presented as a percentage rounded to one decimal place. (E.g., `99.6%`.)
- `potential_contaminating_sample_unambiguous_bases_passing_read_depth_filter`: the number of unambiguous (`A`, `T`, `C`, or `G`) bases in the potentially contaminating samples's consensus genome that pass the [read depth filter](#position-inclusion-thresholds). This column only appears if `--min-depth` is non-zero. (E.g., `29,509`.)
- `potential_contaminating_sample_genome_covered_passing_read_depth_filter`: the number of unambiguous (`A`, `T`, `C`, or `G`) bases in the potentially contaminating samples's consensus genome that pass the [read depth filter](#position-inclusion-thresholds) divided by the number of unambiguous bases in the reference genome. The propotion is presented as a percentage rounded to one decimal place. This column only appears if `--min-depth` is non-zero. (E.g., `98.7%`.)

Columns summarizing potentially contaminating alleles:
- `number_consensus_differences`: the number of positions at which the contaminating sample and contaminated sample consensus genomes differ. Only substitutions are included; insertions and deletions are not included. (E.g., `21`.)
- `consensus_differences`: a list of positions at which the contaminating sample and contaminated sample consensus genomes differ and the alleles at those positions. The putatively contaminated sample's consensus allele is listed first, followed by the putatively contaminating sample's consensus allele. Only substitutions are included; insertions and deletions are not included. (E.g., `2,162 AG; 2,813 AG; 3,044 CA; 5,452 GA; 6,429 TC; 6,762 CT; 7,348 TA; 9,391 TC; 10,702 TC; 13,119 CT; 14,484 CT; 18,131 CT; 19,072 GT; 19,868 CT; 21,203 AG; 24,703 TC; 27,630 TC; 29,095 TC; 29,272 TC; 29,360 TC; 29,367 TC`.)
- `number_consensus_differences_matched_as_minor_alleles`: the number of minor alleles in the contaminated sample that appear in the contaminating sample's consensus genome. (E.g., `21`.)
- `proportion_consensus_differences_matched_as_minor_alleles` the proportion of positions where the two samples' consensus genomes differ at which the contaminating sample's consensus allele appears as a minor allele in the contaminated sample. (E.g., `100%`)
- `minor_alleles_matched`: a list of the alleles at heterozygous positions in the potentially contaminated sample that appear in the potentially contaminating sample's consensus genome. Each allele is presented as its position followed by its base. Alleles are separated by a semicolon and a space (`; `). (E.g., `2,162 G; 2,813 G; 3,044 A; 5,452 A; 6,429 C; 6,762 T; 7,348 A; 9,391 C; 10,702 C; 13,119 T; 14,484 T; 18,131 T; 19,072 T; 19,868 T; 21,203 G; 24,703 C; 27,630 C; 29,095 C; 29,272 C; 29,360 C; 29,367 C`.)

Columns summarizing mismatches between the potential contaminating and contaminated samples:
- `num_mismatches`: the number of positions at which the two consensus genomes differ where the contaminating sample's consensus allele does not appear as a minor allele in the contaminated sample; in other words, the number of potentially contaminating consensus-level alleles that do not appear as a minor or consensus-level allele in the potentially contaminated sample. (E.g., `0`.)
- `mismatches`: a list of potentially contaminating consensus-level alleles that do not appear as a minor or consensus-level allele in the potentially contaminated sample. Each allele is presented as its position followed by its base in the potentially contaminating sample's consensus genome. Alleles are separated by a semicolon and a space (`; `).

Columns summarizing the frequencies of potential contaminating alleles:
- `estimated_contamination_volume`: the median frequency of contaminating alleles in the potentially contaminated sample, including 0%s at positions where the contaminating allele does not appear. The frequency is presented as a percentage rounded to one decimal place. (E.g., `4.9%`.)
- `contaminating_allele_frequency_range`: the minimum and maximum frequencies of contaminating alleles in the potentially contaminated sample, including 0%s at positions where the contaminating allele does not appear. The minimum and maximum are presented as percentages rounded to one decimal place, separated by a ` - `. (E.g., `3.2% - 12.3%`.)
- `contaminating_allele_frequencies`: a list of the frequencies in the potentially contaminated sample of matched potentially contaminating alleles (the alleles at heterozygous positions in the potentially contaminated sample that appear in the potentially contaminating sample's consensus genome). Frequencies are separated by a comma and a space (`, `). (E.g., `3.2%, 3.6%, 3.6%, 3.8%, 3.8%, 3.9%, 4.3%, 4.5%, 4.5%, 4.7%, 4.8%, 4.9%, 5.0%, 5.2%, 5.2%, 5.4%, 5.6%, 6.2%, 7.2%, 9.2%, 12.3%`.)

Columns summarizing plate positions:
- `potential_contaminated_sample_plate_position`: the well occupied by the sample listed in `potential_contaminated_sample`, if a plate map is provided. If multiple plate maps are provided, and the potentially contaminated sample appears in multiple plate maps, this field will contain a list of all wells the potentially contaminated sample appears in, separated by a comma and a space (`, `). (E.g., `H10`.)
- `potential_contaminating_sample_plate_position`: the well occupied by the sample listed in `potential_contaminating_sample`, if a plate map is provided. If multiple plate maps are provided and the potential contaminating sample appears in multiple plate maps, this field will contain a list of all wells the potential contaminating sample appears in, separated by a comma and a space (`, `). (E.g., `H9`.)
- `potential_contaminated_sample_plate`: a list of all plate maps the potentially contaminated sample appears in, in order corresponding to the list of wells in `potential_contaminated_sample_plate_position`. This column is only present if multiple plate maps are provided.
- `potential_contaminating_sample_plate`: a list of all plate maps the potential contaminating sample appears in, in order corresponding to the list of wells in `potential_contaminating_sample_plate_position`. This column is only present if multiple plate maps are provided.

### Plate Map Visualization of Potential Cross-Contamination

`--out-figures`

If you enter at least one plate map file, polyphonia will generate a visualization of potential cross-contamination detected on each plate map entered. Each potentially contaminated well is colored by the median contaminating allele frequency. In cases where a sample is proposed to be contaminated by multiple sources of contamination, each potentially contaminated well is colored by the sum of the median contaminating allele frequency from each proposed source of contamination. This metric is labelled the total estimated contamination volume. Proposed contamination events are shown by arrows to each potentially contaminated well from the proposed source(s) of contamination. Arrow thickness is determined by the number of positions at which the contaminating consensus allele appears as a minor allele in the contaminated sample. Wells proposed to be either contaminated or sources of contamination have a black border, while all other wells have a grey border. Note that sample comparisons are made only according to the entered [well comparison options](#well-comparison-options)—other contamination events would not be detected.

You can view an example of a plate map visualization of potential cross-contamination at the [end](#output-files-1) of the [example run-through](#example-run-throughs):

![USA-MA-Broad_CRSP-01315_23-2021_plate_map txt_potential_cross_contamination_visualization](https://github.com/broadinstitute/polyphonia/assets/6245320/b423ad07-365d-408b-b320-d4a92ef56a00)

You can set the [output directory](#output-file-paths) for the plate map visualization of potential cross-contamination using `--out-figures`.

### Plate Map Visualization of iSNVs

`--out-figures`
`--print-all-iSNVs`

If you enter at least one plate map file, polyphonia will generate a visualization of the plate with each well colored by and labelled with the number of iSNVs, or positions with intrahost variation (base substitutions only) passing our [allele filtering thresholds](#allele-filtering-thresholds) and [read depth filter](#position-inclusion-thresholds).

You might notice that some or even many wells are grey and labelled `NA` ("no data") even if the well had been mapped to a sample in a plate map. This is because processing [within-sample diversity files](#within-sample-diversity-files) can be very slow, and we only do it if a sample has passed [sample inclusion thresholds](#sample-inclusion-thresholds) and has at least one neighbor (which also passes [sample inclusion thresholds](#sample-inclusion-thresholds)) according to the provided [well comparison options](#well-comparison-options). You can include all wells mapped to samples passing [sample inclusion thresholds](#sample-inclusion-thresholds), including samples without neighbors, by adding `--print-all-iSNVs TRUE`.

You can view an example of a plate map visualization of iSNVs at the [end](#output-files-1) of the [example run-through](#example-run-throughs):

![USA-MA-Broad_CRSP-01315_23-2021_plate_map txt_iSNVs_visualization](https://github.com/broadinstitute/polyphonia/assets/6245320/90064057-3ac3-46b4-90fc-ac8d741b2754)

You can set the [output directory](#output-file-paths) for the plate map visualization of iSNVs using `--out-figures`.

## Example Run-Throughs

You can practice using polyphonia by detecting cross-contamination between two of our SARS-CoV-2 clinical samples: plate neighbors [USA-MA-Broad_CRSP-01315-2021](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13971573) and [USA-MA-Broad_CRSP-01323-2021](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR13971565).

We will run polyphonia twice: once with vcf files and unaligned consensus genomes as input, once with heterozygosity tables and aligned consensus genomes as input. (After you are done following along, you can run polyphonia with vcf files and *aligned* consensus genomes and with heterozygosity tables and *unaligned* consensus genomes.)

Before we start, follow [Getting Started](#getting-started) Steps 1-4 to set up your container.

### With VCF Files and Unaligned Consensus Genomes

#### Inputs

Inside your container, the input files we will use for this example are in the directory `/opt/polyphonia/test/input`:

- [**NC_045512.2.fasta**](/test/input/NC_045512.2.fasta): the SARS-CoV-2 reference genome, retrieved from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2).
- [**USA-MA-Broad_CRSP-01315-2021.fasta**](/test/input/USA-MA-Broad_CRSP-01315-2021.fasta) and [**USA-MA-Broad_CRSP-01323-2021.fasta**](/test/input/USA-MA-Broad_CRSP-01323-2021.fasta): the SARS-CoV-2 consensus genomes assembled from samples USA-MA-Broad_CRSP-01315-2021 and USA-MA-Broad_CRSP-01323-2021.
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

Our output table is in `/opt/polyphonia/test/outputs/test1_potential_cross_contamination.txt`:

<img width="1749" alt="output_table" src="https://github.com/broadinstitute/polyphonia/assets/6245320/62faf13a-3f88-448a-bd16-c23a685a1ace">

The output table shows that USA-MA-Broad_CRSP-01323-2021, in well H10, is potentially contaminated by USA-MA-Broad_CRSP-01315-2021, in well H9. USA-MA-Broad_CRSP-01323-2021 has 21 positions with intrahost variation passing our thresholds. Alleles at all of these positions match the USA-MA-Broad_CRSP-01315-2021 consensus genome. There are no mismatches, or USA-MA-Broad_CRSP-01315-2021 consensus-level alleles not appearing as minor or consensus-level alleles in USA-MA-Broad_CRSP-01323-2021. (A note that insertions and deletions are ignored: we only examine base substitutions.)

The USA-MA-Broad_CRSP-01323-2021 alleles appearing at consensus level in USA-MA-Broad_CRSP-01315-2021 range in frequency from 3.2% to 12.3%, with a median allele frequency of 4.9%. In this case, all potential contaminating alleles appear as minor, rather than consensus-level, alleles in USA-MA-Broad_CRSP-01323-2021.

Our plate visualizations are in the directory `/opt/polyphonia/test/test1_plate_visualizations`, `USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt_potential_cross_contamination_visualization.jpg` and `USA-MA-Broad_CRSP-01315_23-2021_plate_map.txt_iSNVs_visualization.jpg`:

![USA-MA-Broad_CRSP-01315_23-2021_plate_map txt_potential_cross_contamination_visualization](https://github.com/broadinstitute/polyphonia/assets/6245320/c6af09dc-0eb3-48ea-89cc-846049800fbb)
![USA-MA-Broad_CRSP-01315_23-2021_plate_map txt_iSNVs_visualization](https://github.com/broadinstitute/polyphonia/assets/6245320/64dcc158-0656-40b6-a1df-6f4bafca0b22)

In the visualization of potential cross-contamination, we see that USA-MA-Broad_CRSP-01323-2021, in well H10, is potentially contaminated by USA-MA-Broad_CRSP-01315-2021, in well H9, and that the median allele frequency of potentially contaminating alleles is 4.9%.

In the iSNVs visualization, we see that USA-MA-Broad_CRSP-01323-2021, in well H10, has 21 positions with intrahost variation (base substitutions only) passing our [allele filtering thresholds](#allele-filtering-thresholds), and that USA-MA-Broad_CRSP-01315-2021, in well H9, has 3. We did not examine intrahost variation in any of the other wells.

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
- [**USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf_heterozygosity.txt**](/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf_heterozygosity.txt) and [**USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf_heterozygosity.txt**](/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf_heterozygosity.txt): tables summarizing within-sample heteroyzgosity detected by [`LoFreq call`](https://csb5.github.io/lofreq/commands/#call) in USA-MA-Broad_CRSP-01315-2021 and USA-MA-Broad_CRSP-01323-2021. In our previous example, these tables were generated as intermediate files after processing [USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf](/test/input/USA-MA-Broad_CRSP-01315-2021.bam_LoFreq.vcf) and [USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf](/test/input/USA-MA-Broad_CRSP-01323-2021.bam_LoFreq.vcf). (Note that not all rows will pass our filters—most will not be included in the output.)
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

If you are using polyphonia through [Docker](#docker), you do not need to worry about dependencies.

You can view the software polyphonia uses in [`requirements-conda.txt`](https://github.com/broadinstitute/polyphonia/blob/main/requirements-conda.txt):
- Scripts are written in [Perl 5](https://www.perl.org/get.html).
   - Output float values are rounded using [Math::Round](https://metacpan.org/pod/Math::Round) ≥0.07.
   - Parallelization is performed using [Parallel::ForkManager](https://metacpan.org/pod/Parallel::ForkManager) ≥2.02.
- Allele frequencies are calculated by [LoFreq](https://csb5.github.io/lofreq/installation/) ≥2.1.5.
- Sequence alignments are generated by [MAFFT](https://mafft.cbrc.jp/alignment/software/) ≥7.480.
- Read depths are calculated by [SAMtools depth](http://www.htslib.org/doc/samtools-depth.html).
- Visualizations are created in [R](https://www.r-project.org/) ≥4.0.5 using [tidyverse](https://www.tidyverse.org/) ≥1.2.1.

### Name

[**Polyphony**](https://www.youtube.com/watch?v=teh22szdnRQ) describes music containing multiple simultaneous voices with their own melodies, like parts of ["We Will Fall Together" by Streetlight Manifesto](https://www.youtube.com/watch?v=SOqenYis1iQ), ["Somebody That I Used to Know" by Walk off the Earth](https://www.youtube.com/watch?v=d9NF2edxy-M), ["exile" by Taylor Swift with Bon Iver](https://www.youtube.com/watch?v=osdoLjUNFnA), ["Daddy Cool" by Boney M](https://www.youtube.com/watch?v=tHyHOUShOg8), ["Come As You Are" by Nirvana](https://www.youtube.com/watch?v=vabnZ9-ex7o), ["Severed" by The Decemberists](https://www.youtube.com/watch?v=ksTFj6L0mao), ["That's Not My Name" by The Ting Tings](https://www.youtube.com/watch?v=v1c2OfAzDTI), and to at least some extent most music created in recent memory.

[**Diafonia**](https://translate.google.com/?sl=it&tl=en&text=diafonia&op=translate) means crosstalk in Italian, Portuguese, and Haitian Creole, and is very close to the words for crosstalk in Spanish (diafonía) and French (diaphonie). In Italian, diafonia also refers to the simultaneous presence of multiple voices, in music and otherwise.

In addition to being the name of our software package, **Polyphonia** is also [a particularly beautiful abstract ballet](https://www.nycballet.com/discover/ballet-repertory/polyphonia/) created by Christopher Wheeldon for New York City Ballet. I think you will find it captures the essence of cross-contamination quite nicely. You can watch an excerpt by clicking the preview image below:

[<img width="1157" alt="ballet_video_preview" src="https://user-images.githubusercontent.com/6245320/121761826-3421c900-cb00-11eb-8712-544ab801fe3c.png">](http://www.youtube.com/watch?v=MOhvllhQo8A)
