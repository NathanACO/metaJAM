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

## Input files specified in the config `nextflow.config` file
### 1. If running the pipeline from scratch
Path of the raw sequencing samples to be processed: 
`FASTQ_direct_path='/path/to/*_{R1,R2,1,2,R1_001,R2_001}*.{fastq,fq}.gz'` OR/AND `FASTQ_list_path="/path/to/List_fastq.txt"`. Please leave the other as "" when you only use one format. You can also use both formats, and the two sets of fastqs will all be processed.

Give absolute path to the overall output directory: `OUTPUT_Dir="\path\to"`. It will creates different folders for the different steps of the pipeline inside your given `OUTPUT_Dir`, where the outputs of each step will be stored correspondingly:
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

1.2 Or If running a specific tool or the pipeline from any step after preprocessing
File containing a list of samples with absolute path to be processed (OVERRIDE_LIST_*) [In development]

### 2. Path to databases
### 3. Tools activation with "enable" or "disable"
### 4. Parameters to precise for specific tools
1) *fastp*\
-overlap_len_require        (default=20)\    FASTP_OVERLAP_LEN_REQUIRE
-l                          (default=30)    FASTP_MIN_LENGTH
2) *SGA*\
--dust-threshold            (default=4)    SGA_DUST_THRESHOLD
3) *PRINSEQ*\
-lc_method                  (default=dust)\  PRINSEQ_COMPLEXITY_METHOD
-lc_threshold               (default=4)\   PRINSEQ_COMPLEXITY_THRESHOLD
-min_len                    (default=35)\  PRINSEQ_MIN_LE
-derep                      (default=1)    PRINSEQ_DEREP
4)*Bowtie2*
-k                          (default=1000)    BOWTIE2_N_ALLOW_MULTIMAPPER
5) *bamdam*\
--stranded                  (default=ds)\        BAMDAM_STRANDED
--minreads                  (default=5)\        BAMDAM_MINREADS
--maxdamage                 (default=0.5)\    BAMDAM_MAXDAMAGE
TOP_GENUS                   (default=10)        BAMDAM_TOP_GENUS   # Number of the most abundant genus to plot for damage\
5.1) *MMSeqs2*\
MMSEQS2_THREADS             (default=60)\
MMSEQS2_MAX_SEQS            (default=300)\
MMSEQS2_MIN_LENGTH          (default=30)\
MMSEQS2_MIN_SEQID                   (default=0.93)\
MMSEQS2_MIN_BITS                    (default=50)\
MMSEQS2_MIN_QUERY_COV               (default=0.95)\
MMSEQS2_MAX_EVALUE                  (default="1e-5")\
MMSEQS2_S                   (default=7.5)\
MMSEQS2_SPACED_KMER_MODE    (default=1)\
MMSEQS2_SPLIT_MEM_LIMIT     (default=220G)\
5.2) *MMSeqs2 evaluation*\
MMSEQS2_TOP_GENERA          (default=10)\
MMSEQS2_GENERA_FILE         # path to a file containing genera of interst at each line (?); otherwise use 'assets/NO_FILE' 
MMSEQS2_MIN_DMG             (default=3.5)\
MMSEQS2_MAX_READS           (default=100)\
MMSEQS2_MIN_READS           (default=30\
MMSEQS2_SEED                (default=42)            # Leave empty to have different reads each run, othewise set an integer (e.g. 42)\
MMSEQS2_AMBIG_FRAC          (default=0.05)\         # Percent of e-value difference for the second best hit for a different genera than the mmseqs2 best hit 

7) *Plots*\
PLOTS_BAMDAM_MIN_READS (default=50)    - Minimum reads per sample to include in bamdam plots\
PLOTS_BAMDAM_PLOT_MODE (default=both)  - Chose which plots to produce: heatmap, bubble or both\
PLOTS_DAMAGE_THRESHOLD (default=3.5)   - Minimum percentage of damage for plotting
PLOTS_PLOT_LOW_DAMAGE_TAXA (default=0)    # 1 = keep low-damage taxa (default); 0 = drop taxa whose max damage across samples is < PLOTS_DAMAGE_THRESHOLD
PLOTS_EXCLUDE_TAXA # e.g., "Homo;Zea;Canis;Veronica" in the format of comma / semicolon / space separated (exact taxon name matches)
PLOTS_KRONA=1                   # Plots the bamdam results for all samples in metadata if enable. 1=enable, 0=disable
PLOTS_LIST_TAXA_EVOLUTION_FILE # Give a list of taxa with one taxa per line, to produce plot of abundance of each taxa as a line representation

### 5. SBATCH task parameter for memory and time assigned to the processing using tools
To be refine based on samples size and database using
> [!TIP]
> We advise to set up at least 650G and 256CPUs for the Kraken step, and to allow at least 20 hours of running time for massive fastq.gz (>60M reads)\
> For smaller datasets (around 10-20M), 15h for 3-4 samples seems a coherent value, depending on your cluster queueing system.

Apart from the memory and time for each process, you also need to change the project number if you are using slurm task system in HPC

```
params {
    project                    = "naiss2025-xx-xxx" // for example for Naiss project allocation
```
## Required program
nextflow (developed in v25.10.3)\
conda\
taxadb #use to download a file for setup

run the script once before running anything
```
mkdir assets
touch assets/NO_FILE
cd assets
python -m venv taxadb
source taxadb/bin/activate
pip3 install taxadb
taxadb download -o taxadb 
taxadb create -i taxadb --dbname taxadb.sqlite 
cd ../ 
```

## How to launch it:
`nextflow run main.nf -profile conda`\
if you want to resume, add also `-resume`, and `--with-trace` for output a trace*.txt reporting memory and time for each process.

> [!TIP]
> You can run it with the Test samples present in the test folder of this github.

## Plot output
Different plots are created by metaJAM and required a metadata file (see in test for metadata format and content)[In development].\
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

&copy; Ernst Johnson; Chenyu Jin; Jamie Alumbaugh; Nathan Martin | Centre for Palaeogenetics - Stockholm University
