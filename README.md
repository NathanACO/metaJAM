<img src="https://github.com/NathanACO/metaJAM/blob/main/metaJAM_logo.png" width="200" /> 

# metaJAM v1.1.1

Metagenomic Pipeline for ancient DNA analysis performed at the Centre for Palaeogenetics, Stockholm University.

This pipeline performs processing and analysis of metagenomic data, starting from paired-end fastq files or any intermediate files used in the pipeline if not the entire workflow is desired to run. 

> [!NOTE]
> Further developments will be included in metaJAM v1.2.1 (Release 03.02.2026):
> - Conversion of the pipeline into Nextflow to ensure compatability across different HPC clusters
> - Extra module to mask databases based on microbial content
> And in metaJAM v1.2.2 in early March:
> - Addition of leeHom as an alternative choice to fastp
> - Choice of the mapping strategy by iteration or not


## Overview of the pipeline
![alt text](https://github.com/NathanACO/metaJAM/blob/main/metaJAM_diagram.png)

## How to launch it:
`nextflow run main.nf -profile conda`\
if you want to resume, add also `-resume`, and `--with-trace` for output a trace*.txt reporting memory and time for each process.
It will creates different folders for the different steps of the pipeline, where the files will be stored:
- 00_Samples_prefix -> Containing files with list of processed samples at each step 
- 01_fastp
- 02_sga (or 02_prinseq)
- 03_kraken2_filter
- 04_mapping
- 05_merged_bam
- 06_masked_bam 
- 07_ngslca
- 08_bamdam
- 09_kronaTools
- 99_metrics

> [!TIP]
> You can run it with the Test samples present in the test folder of this github.

A few requirements are needed to run this pipeline, and only the config `nextflow.config` file need to be modify in order to run it.

## Required program
nextflow (developed in v25.10.3)\
conda 

## Input files
1. If running the pipeline from scratch
Path of the raw sequencing samples to be processed

2. If running a specific tool or the pipeline from any step after preprocessing
File containing a list of samples with absolute path to be processed

## Path to scripts and databases
So far metaCPG is configured to be run through different databases iteratively.
Please specify any databases that you want to use in the MAP_DB_LIST variable of the config file.
If the same samples are planning to be run against different databases, we strongly advise to define a name for each combination of databases in the MAP_LAST_DB_TAG variable of the config.

## Tools activation with "enable" or "disable"

## Parameters to precise for specific tools
1. *fastp*\
-overlap_len_require        (default=20)\
-l                          (default=30)
2. *SGA*\
--dust-threshold            (default=4)
3. *PRINSEQ*\
-lc_method                  (default=dust)\
-lc_threshold               (default=4)\
-min_len                    (default=35)\
-derep                      (default=1)
4. *bamdam*\
--stranded                  (default=ds)\
--minreads                  (default=5)\
--maxdamage                 (default=0.5)\
TOP_GENUS                   (default=10)           # Number of the most abundant genus to plot for damage\
5.1 *MMSeqs2*\
MMSEQS2_THREADS             (default=60)\
MMSEQS2_MAX_SEQS            (default=300)\
MMSEQS2_MIN_LENGTH          (default=30)\
MIN_SEQID                   (default=0.93)\
MIN_BITS                    (default=50)\
MIN_QUERY_COV               (default=0.95)\
MAX_EVALUE                  (default="1e-5")\
MMSEQS2_S                   (default=7.5)\
MMSEQS2_SPACED_KMER_MODE    (default=1)\
MMSEQS2_SPLIT_MEM_LIMIT     (default=220G)\
5.2 *MMSeqs2 evaluation*\
MMSEQS2_TOP_GENERA          (default=10)\
MMSEQS2_MIN_DMG             (default=3.5)\
MMSEQS2_MAX_READS           (default=100)\
MMSEQS2_MIN_READS           (default=30\
MMSEQS2_SEED                (default=42)            # Leave empty to have different reads each run, othewise set an integer (e.g. 42)\
MMSEQS2_AMBIG_FRAC          (default=0.05)\         # Percent of e-value difference for the second best hit for a different genera than the mmseqs2 best hit 

6. *Plots*\
PLOTS_BAMDAM_MIN_READS (default=50)    - Minimum reads per sample to include in bamdam plots\
PLOTS_BAMDAM_PLOT_MODE (default=both)  - Chose which plots to produce: heatmap, bubble or both\
PLOTS_DAMAGE_THRESHOLD (default=3.5)   - Minimum percentage of damage for plotting

## SBATCH parameters
To be refine based on samples size and database using
> [!TIP]
> We advise to set up at least 650G and 256CPUs for the Kraken step, and to allow at least 20 hours of running time for massive fastq.gz (>60M reads)\
> For smaller datasets (around 10-20M), 15h for 3-4 samples seems a coherent value, depending on your cluster queueing system.

## Plot output
Different plots are created by metaJAM and required a metadata file (see in test for metadata format and content).\
The first plot created is giving an overview of the metrics of the samples processed:
![alt text](https://github.com/NathanACO/metaJAM/blob/main/reads_per_step_dots.png)
&copy; Martin et al. 2025 plot, peatbog paper - in revision

The second plot is giving the taxa composition as a bubble plot or a heatmap (can be chosen in the config file).\
The user can filter the minimum number of reads required to plot a taxa in the config file.\
It will produce two different plots each time, one for the genus and one for the family, plotting in alphabetical order the Viridiplantae taxa first and then the Opisthokonta taxa.

![alt text](https://github.com/NathanACO/metaJAM/blob/main/bamdam_family_bubbleplot_reads.png)
![alt text](https://github.com/NathanACO/metaJAM/blob/main/bamdam_family_heatmap.png)
&copy; Inspired by Tyler Murchie and Scott Cocker plot - in revision

> [!NOTE]
> An extra plot containing information about the damage is also produced for supplementary information.\
> The damage is coded by 3 colours (red,orange,green).\
> Green if for this taxa the damage are superior to 5% and if the damage percent is within the interval (+-5%) of the mean of the damage of the 3 main taxa that have more than 5% damage.\
> If only one condition is met, the associated damage column will be orange and the user will be invited to investigate deeper this taxa.
> 
> ![alt text](https://github.com/NathanACO/metaJAM/blob/main/bamdam_family_bubbleplot_damage.png)


_Why metaJAM?_\
We like to see ancient sediment metagenomic data as a jar of jam:
- You never really know what is inside
- It is most of the time quite old and degraded
- We are always seeking for some help to open the jar and define which fruits were used to make it!

&copy; Ernst Johnson; Jamie Alumbaugh; Nathan Martin | Centre for Palaeogenetics - Stockholm University
