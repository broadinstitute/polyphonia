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
- [Required Input Files](#getting-started)
- [Optional Plate Map Inputs](#optional-plate-map-inputs)
- [Other Options](#other-options)
- [Output Files](#output-files)
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

## Required Input Files

TODO

## Optional Plate Map Inputs

TODO

## Other Options

TODO

## Output Files

TODO

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

In the cross-contamination visualization, we see that USA-MA-Broad_CRSP-01323-2021, in well H10, is potentially contaminated by USA-MA-Broad_CRSP-01315-2021, in well H9, and that the median allele frequency of potentially contaminating alleles is 4.9%.

In the iSNVs visualization, we see that USA-MA-Broad_CRSP-01323-2021, in well H10, has 21 positions with intrahost variation (base substitutions only) passing our thresholds, and that USA-MA-Broad_CRSP-01315-2021, in well H9, has 3. We did not examine intrahost variation in the other wells.

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
