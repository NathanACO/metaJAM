
## Input files specified in the config `nextflow.config` file

### 1. Tools activation with "enable" or "disable"
Indentation of the process suggest that the process is included in the process at the line above. For example, to enable FASTP you need to enable preprocess. 
```
ENABLE_PREPROCESS="enable"     // If disabled, neither fastp or SGA will run\
    ENABLE_FASTP="enable"\
    ENABLE_SGA="enable"\
    ENABLE_PRINSEQ="enable"        // enable to use PRINSEQ instead of SGA         \
    ENABLE_LOW_COMPLEXITY_FILTER="PRINSEQ" // use 'PRINSEQ' prinseq result for downstream analysis; Use 'SGA' to use SGA result for downstream analysis; Use '' if BOTH are disabled

ENABLE_KRAKEN_GTDB="enable"\
ENABLE_MAPPING="enable"\
    ENABLE_PARALLEL_MAPPING='enable'\
    ENABLE_SEQUENTIAL_MAPPING='enable'\
    USE_MAPPING='PARALLEL_MAPPING' // use the output bam files from 'PARALLEL_MAPPING' or 'SEQUENTIAL_MAPPING' for downstream analysis, it cannot be '' when mapping is enabled.\

ENABL_MERGE_BAM="enable" // merge the bam of the same sample mapped against different databases

ENABLE_MASK_REGIONS="enable"\
    ENABLE_GENERATE_BEDFILE_TO_MASK="enable" // can be enabled only if 'ENABLE_MASK_REGIONS' is enabled\
        ENABLE_MASK_FASTA="enable" // can be enabled only if the subworkflow 'ENABLE_GENERATE_BEDFILE_TO_MASK' is enabled

ENABLE_NGSLCA="enable"\
ENABLE_BAMDAM="enable"\
ENABLE_KRONATOOLS="enable"\
ENABLE_MMSEQS2="enable" // enable this requires ENABLE_BAMDAM\
ENABLE_METRICS="enable"\
ENABLE_PLOTS="enable"
```

### 2. Input files
2.1 Path of the raw sequencing samples to be processed: 
`FASTQ_direct_path='/path/to/*_{R1,R2,1,2,R1_001,R2_001}*.{fastq,fq}.gz'` OR/AND `FASTQ_list_path="/path/to/List_fastq.txt" (each row: ID[tab]fastq1[tab]fastq2)`.\ 
Please leave the other as "" when you only use one format. You can also use both `FASTQ_direct_path` and `FASTQ_list_path`, and both sets of fastqs will all be processed.

Provide `metadata` for plotting, you can see ./test/metadata.txt as an example:\
`metadata="./test/metadata.txt"`\
    #the columns are spaced with tab:\
    sample	age_ka	depth_cm	sample_type	layer	notes	site\
    X	1	25	    sediment	NA	Good DNA	LakeA

2.2 If you are (1) running after any process (all the previous processes are not 'enable') or (2) supplements data which already run through the process to the new batch of samples you are running, you can input a file containing a list of samples with absolute path to be processed (OVERRIDE_LIST_*). The format of each row in the input file is specified at the end of each line.

- specify a file if you are skipping fastq, otherwise use ""\
`OVERRIDE_LIST_FASTP="/path/to/file"` // sample_ID[tab]merged_fastq

- specify a file if you are skipping the preprocessing meant for merging and low-complexity reads removal, otherwise use ""\
`OVERRIDE_PREPROCESSED="/path/to/file"` // sample_ID[tab]preprocessed_fastq

- specify a file if you are skipping filtering out microbial reads mapped against contamination(or microbial) database, otherwise use ""\
`OVERRIDE_LIST_KRAKEN="/path/to/file"` // sample_ID[tab]unclassified_fastq

- specify a file if you are skipping mapping with bowtie2, otherwise use ""\
`OVERRIDE_LIST_BAM="/path/to/file"` // give a file with each line: sample_ID[tab]absolute_path_to_sample_bam
    - if you would not like to merge among the mapping result against different databases, make sure sample_ID is sample_ID+database_name for the workflow to differentiate them. 
    - Otherwise in default all bam of the same sample are always merged and analyzed together later.

- specify a file if you are skipping ngsLCA, otherwise use ""\
`OVERRIDE_LIST_NGSLCA="/path/to/file"` // sample_ID[tab]sample_lca_file

- specify a file if you are skipping bamdam, otherwise use ""\
`OVERRIDE_LIST_BAMDAM="/path/to/file"` // sample_ID[tab]bam[tab]lca[tab]tsv[tab]xml_file

- specify a file if you are skipping mmseq2, otherwise use ""\
`OVERRIDE_LIST_MMSEQ2=""` // ID[tab]mmseq2_output_evaluation_file
    - here the mmseq2_output_evaluation_file is output by a summary script after running mmseq2 and has suffix *.bamdam_mmseqs.evaluation.summary.tsv

- specify a file if you are skipping running the script to generate metrics, otherwise use ""\
`OVERRIDE_LIST_METRICS="/path/to/file"`
  - the header is "sample[tab]count_raw_fq1[tab]count_raw_fq2[tab]count_merged_fq[tab]count_rm_low_complex_fq[tab]count_k2_mic_unclas_fq[tab]bam_header[tab]bamdam_bam" and the values in the following rows locate correspondingly
    
2.3 Give absolute path to the overall output directory: `OUTPUT_Dir="\path\to"`. It will create different folders for the different steps of the pipeline inside your given `OUTPUT_Dir`, where the outputs of each step will be stored correspondingly:
01_fastp -> Containing trimmed and merged pair-end read files
02_sga/02_prinseq -> Containing read files where low-complexity sequences have been removed
03_kraken2_filter -> Containing the unclassified and classified microbial read files
04_parallel_mapping/04_sequential_mapping -> Containing the mapped and unmapped read files
acc2taxid (optional) -> all acc2taxid files generated
05_merged_bam -> bam of sample sample mapped against different databases are merged and output here
06_masked_bam -> bam with microbial/contamination-like regions masked
07_ngslca -> ngsLCA output after assigning multimappers to the lowest common ancestor taxonomical node
08_bamdam -> bamdam output for quality filtering and plot kronaplot for individual sample
09_mmseq2 -> MMseq2 output after blasting the results
10_metrics -> Containing metrics that are the number of reads after each process
11_plots -> kronaplot of samples of each site, heatmap/bubbleplot summarizing the MMseq2 results or DNA damage along with the number of reads
note that krona plot of each sample generated by bamdam will be in `08_bamdam/`, while kronaplots of the each site and other summary plots of MMseq2 result or DNA damage are in `11_plots/`.

2.4 The path to where you download the metajam:
`METAJAM_DIR`="/path/to/metaJAM" 


### 3. Path to databases
for DB[1-10], set as "" if not used
`BOWTIE2_MAPPING_DB1="/path/to/index_name" `: each line as the bowtie2 mapping index header (e.g., '/path/to/header' where the header refers to header.*.bt2*), the first database reads are aligned to if you specified `USE_MAPPING='SEQUENTIAL_MAPPING`, otherwise the order of database does not matter.
`BOWTIE2_MAPPING_DB2`...`BOWTIE2_MAPPING_DB10`: same format as DB1, maximum 10 databases so far.

`KRAKEN2_FILTER_DATABASE="/path/to/kraken2_db" `: path to Kraken2 database (e.g., GTDB database), the directory containing kraken2 indexes (*.k2d)\


```
NAMES="assets/names.dmp"
NODES="assets/nodes.dmp"
ACC2TAXID="/path/to/acc2taxid.txt" #each line in the format: contig[tab]contig[tab]NCBI_taxonomy_ID
```
NCBI taxonomy files (names.dmp and nodes.dmp downloaded in the `assets/` dir in the preparation script given in this tutorial)

input a bed file with all regions to be masked\
`REGIONS_TO_MASK="/path/to/databaseA.bed"`\
- if aiming to mask microbial-like regions without supplying the bed file for masking, we could also generate and use bed file using by specify `ENABLE_GENERATE_BEDFILE_TO_MASK="enable"` and fill in requirement of GENEX workflow

GENEX workflow to generate bedfile for microbial/contamination-like region for reference genomes\
`MCWORKFLOW_input_dir= ""`  // absolute path to the dir containing the reference genome fasta\
`MCWORKFLOW_input_list=""`  // a list with each line as the absolute path to the reference genome fasta\
`MCWORKFLOW_pseudo_reads_file_dir=""` // a directory to store the pseudo reads generated of GTDB/other suspected containmation source\
`MCWORKFLOW_type_of_pseudo_reads="GTDB"` // label of pseudo reads\
`MCWORKFLOW_n_allowed_multimappers=1000`

MMSeqs2 database\
`MMSEQS2_DB="/sw/data/MMseqs2_data/latest/rackham/NT"`\
`MMSEQS2_TAXADB_SQLITE="assets/taxadb_nucl.sqlite"` \\ taxadb_nucl.sqlite is in the `assets/` dir in the preparation script given in this tutorial


### 4. Parameters to precise for specific tools
1) *fastp*\
`FASTP_OVERLAP_LEN_REQUIRE`    -overlap_len_require        (default=20)    \
`FASTP_MIN_LENGTH`    -l                          (default=30)    
2) *SGA*\
`SGA_DUST_THRESHOLD`    --dust-threshold            (default=4)    
3) *PRINSEQ*\
`PRINSEQ_COMPLEXITY_METHOD`    -lc_method                  (default=dust)  \
`PRINSEQ_COMPLEXITY_THRESHOLD`    -lc_threshold               (default=4)   \
`PRINSEQ_MIN_LEN`    -min_len                    (default=35)  \
`PRINSEQ_DEREP`    -derep                      (default=1)    \
4) *Bowtie2*
`BOWTIE2_N_ALLOW_MULTIMAPPER`    -k                          (default=1000)    
5) *bamdam*\
`BAMDAM_STRANDED`    --stranded                  (default=ds)        
`BAMDAM_MINREADS`    --minreads                  (default=5)        
`BAMDAM_MAXDAMAGE`    --maxdamage                (default=0.5)    
`BAMDAM_TOP_GENUS`    TOP_GENUS                   (default=10)         Number of the most abundant genus to plot for damage\
5.1) *MMSeqs2*\
`MMSEQS2_THREADS`             (default=60)\
`MMSEQS2_MAX_SEQS`            (default=300)\
`MMSEQS2_MIN_LENGTH`          (default=30)\
`MMSEQS2_MIN_SEQID`                   (default=0.93)\
`MMSEQS2_MIN_BITS`                    (default=50)\
`MMSEQS2_MIN_QUERY_COV`               (default=0.95)\
`MMSEQS2_MAX_EVALUE`                  (default="1e-5")\
`MMSEQS2_S`                   (default=7.5)\
`MMSEQS2_SPACED_KMER_MODE`    (default=1)\
`MMSEQS2_SPLIT_MEM_LIMIT`     (default=220G)\
5.2) *MMSeqs2 evaluation*\
`MMSEQS2_TOP_GENERA`          (default=10)\
`MMSEQS2_GENERA_FILE`         # path to a file containing genera of interst at each line (?); otherwise use 'assets/NO_FILE' 
`MMSEQS2_MIN_DMG`             (default=3.5)\
`MMSEQS2_MAX_READS`           (default=100)\
`MMSEQS2_MIN_READS`           (default=30\
`MMSEQS2_SEED`                (default=42)            # Leave empty to have different reads each run, othewise set an integer (e.g. 42)\
`MMSEQS2_AMBIG_FRAC`          (default=0.05)\         # Percent of e-value difference for the second best hit for a different genera than the mmseqs2 best hit 

7) *Plots*\
`PLOTS_BAMDAM_MIN_READS` (default=50)    - Minimum reads per sample to include in bamdam plots\
`PLOTS_BAMDAM_PLOT_MODE` (default=both)  - Chose which plots to produce: heatmap, bubble or both\
`PLOTS_DAMAGE_THRESHOLD` (default=3.5)   - Minimum percentage of damage for plotting\
`PLOTS_PLOT_LOW_DAMAGE_TAXA` (default=0)    # 1 = keep low-damage taxa (default); 0 = drop taxa whose max damage across samples is < PLOTS_DAMAGE_THRESHOLD\
`PLOTS_EXCLUDE_TAXA` # e.g., "Homo;Zea;Canis;Veronica" in the format of comma / semicolon / space separated (exact taxon name matches)\
`PLOTS_KRONA`=1                   # Plots the bamdam results for all samples in metadata if enable. 1=enable, 0=disable\
`PLOTS_LIST_TAXA_EVOLUTION_FILE` # Give a list of taxa with one taxa per line, to produce plot of abundance of each taxa as a line representation\
`PLOTS_TAXA_PER_PLOT`=50       # Number of taxa passing the filters to plot per graph
`MAP_LAST_DB_TAG`="mam-bird-fish_v2" #specify database name for plotting
`PLOTS_SAMPLES_FOR_PLOTS`=params.METAJAM_DIR+"assets/NO_FILE8" // if you want to plot only a subset of samples, then supply a subset of metadata containing the samples that you want to plot. If no subsetting is needed, then use the placeholder "assets/NO_FILE8"


### 5. SBATCH task parameter for memory and time assigned to the processing using tools

The process will automatically retry for 6 times and each time with more cpu and memory. However, if it still fails, you can modify the cpus, memory accordingly and resume with `nextflow run ... -resume`

```
      withName: PRINSEQ {
         cpus={check_max(62 * task.attempt, 'cpus' ) }
         memory={ check_max(50.GB * task.attempt, 'memory' ) }
         time={ check_max(8.h * task.attempt, 'time' ) }
     }
```

To be refine based on samples size and database using
> [!TIP]
> We advise to set up at least 650G and 256CPUs for the Kraken step, and to allow at least 20 hours of running time for massive fastq.gz (>60M reads)\
> For smaller datasets (around 10-20M), 15h for 3-4 samples seems a coherent value, depending on your cluster queueing system.

Apart from the memory and time for each process, you also need to change the project number if you are using slurm task system in HPC

```
params {
    project                    = "naiss2025-xx-xxx" // for example for Naiss project allocation
```
