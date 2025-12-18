<img src="https://github.com/NathanACO/metaJAM/blob/main/metaJAM_logo.png" width="200" /> 

# metaJAM v1.0.1

Metagenomic Pipeline for ancient DNA analysis performed at the Centre for Palaeogenetics - Stockholm

> [!CAUTION]
> Ongoing fixes for the MMSeqs2 part - It won't work if you run it

> [!NOTE]
> Further developments will be included in metaJAM v1.0.2 soon:
> - Addition of leeHom as an alternative choice to fastp
> - Addition of parameters to define in fastp, sga, prinseq and bamdam
> - Different features to chose to concatenate diverse plant database or run specific databases iteratively
> - Choice of the mapping strategy by iteration or not
> - Extra pipeline to mask databases will be added in the future

## How to launch it:
sbatch -x Run_pipeline_metage.sh config_pipeline_metage.sh\
It will creates different folders for the different steps of the pipeline, where the files will be stored:
- 00_logs -> Containing files with list of processed samples at each step
- 01_fastp
- 02_sga
- 03_kraken_gtdb
- 04_mapping
- 05_filtering
- 99_metrics
- log -> Containing one error and one out folder where the sbatch information are stored

> [!TIP]
> You can run it with the Test samples present in the test folder of this github.

A few requirements are needed to run this pipeline, and only the config file need to be modify in order to run it.

## Required modules
1. Modules from Dardel:
- Fastp -v. 0.24+
- Prinseq lite - v.0.20.4+
- Kraken2 - v.2.1.2+
- Bowtie2 - v.2.5.4+
- Samtools - v1.20+
- Kronatools - v2.8.1+

2. Conda environments:
- SGA
- FilterBAM
- ngsLCA

## Input files
1. If running the pipeline from scratch
Path of the raw sequencing samples to be processed

2. If running a specific tool or the pipeline from any step after preprocessing
File containing a list of samples with absolute path to be processed

## Path to scripts and databases
So far metaCPG is configured to be run through 4 different databases:
1. PhyloNorway
2. RefSeq Plastid
3. RefSeq Mito
4. Custom Database - containing any specific taxa of interest, build by the user

> [!CAUTION]
> So far metaJAM is designed to run with 4 databases present together, you can play around only with the Custom Database. \
> New version of metaJAM 1.0.3 will include the choice of using 1 to x databases.

## Tools activation
Precise 1 or 0 for each step, 1=enable, 0=disable

## Parameters to precise for specific tools
1. *fastp*\
-overlap_len_require   (default=20)\
-l                     (default=30)
2. *SGA*\
--dust-threshold       (default=4)
3. *PRINSEQ*\
-lc_method             (default=dust)\
-lc_threshold          (default=4)\
-min_len               (default=35)\
-derep                 (default=1)
4. *bamdam*\
--stranded             (default=ds)\
--minreads             (default=5)\
--maxdamage            (default=0.5)\
TOP_GENUS              (default=10)    - Number of the most abundant genus to plot for damage
6. *Plots*\
PLOTS_BAMDAM_MIN_READS (default=50)    - Minimum reads per sample to include in bamdam plots\
PLOTS_BAMDAM_PLOT_MODE (default=both)  - Chose which plots to produce: heatmap, bubble or both\
PLOTS_DAMAGE_THRESHOLD (default=3.5)   - Minimum percentage of damage for plotting

## SBATCH parameters
To be refine based on samples size and database using
> [!TIP]
> We advise to set up at least 650G and 256CPUs for the Kraken step, and to allow at least 20 hours of running time for massive fastq.gz (>60M reads)\
> For smaller datasets (around 10-20M), 15h for 3-4 samples seems a coherent value, depending on your cluster queueing system.

## Overview of the pipeline
![alt text](https://github.com/NathanACO/metaJAM/blob/main/metaJAM_diagram.png)

## Plot output
Two different plots are created by metaJAM and required a metadata file (see in test for metadata format and content).\
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
> The damage will be coded by 3 colours (red,orange,green).\
> Green if for this taxa the damage are superior to 5% and if the damage percent is within the interval (+-5%) of the mean of the damage of the 3 main taxa that have more than 5% damage.\
> If only one condition is met, the associated damage column will be orange and the user will be invited to investigate deeper this taxa.
> 
> ![alt text](https://github.com/NathanACO/metaJAM/blob/main/bamdam_family_bubbleplot_damage.png)
> 
> Also add a feature of the maximum number of taxa to be plotted together on the graph and the maximum number of sites will be implemented soon.


_Why metaJAM?_\
We like to see ancient sediment metagenomic data as a jar of jam:
- You never really know what is inside
- It is most of the time quite old and degraded
- We are always seeking for some help to open the jar and define which fruits were used to make it!

&copy; Ernst Johnson; Jamie Alumbaugh; Nathan Martin | Centre for Palaeogenetics - Stockholm University
