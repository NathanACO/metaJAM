# metaJAM
Metagenomic Pipeline for ancient DNA analysis performed at the Centre for Palaeogenetics - Stockholm

A few requirements are needed to run this pipeline, and only the config file need to be modify in order to run it.

# Required modules
1. Modules from Darde:\
- Fastp -v. 0.24+\
- Prinseq lite - v.0.20.4+\
- Kraken2 - v.2.1.2+\
- Bowtie2 - v.2.5.4+\
- Samtools - v1.20+\
- Kronatools - v2.8.1+

2. Conda environments:\
- SGA\
- FilterBAM\
- ngsLCA\

# Input files
1. If running the pipeline from scratch\
Path of the raw sequencing samples to be processed

2. If running a specific tool or the pipeline from any step after preprocessing\
File containing a list of samples with absolute path to be processed

# Path to scripts and databases
So far metaCPG is configured to be run through 4 different databases:\
1. PhyloNorway\
2. RefSeq Plastid\
3. RefSeq Mito\
4. Custom Database - containing any specific taxa of interest, build by the user

# Tools activation
Precise 1 or 0 for each step, 1=enable, 0=disable

# Parameters to precise for specific tools
* fastp *
FASTP_OVERLAP_LEN_REQUIRE=20   # --overlap_len_require\
FASTP_MIN_LENGTH=30            # -l\
* SGA *
SGA_DUST_THRESHOLD=4           # --dust-threshold\
* PRINSEQ params *
PRINSEQ_COMPLEXITY_METHOD="dust"\
PRINSEQ_COMPLEXITY_THRESHOLD=4\
PRINSEQ_MIN_LEN=35\
PRINSEQ_DEREP="14"                      # exact fwd+rev duplicates\
* Mapping run name *
#MAP_RUN="default_run"\
* bamdam parameters *
BAMDAM_STRANDED="ds"          # "ds" or "ss"\
BAMDAM_MINREADS=5             # for bamdam krona --minreads\
BAMDAM_MAXDAMAGE=0.5          # for bamdam krona --maxdamage\

# SBATCH parameters
To be refine based on samples size and database using
