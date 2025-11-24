# metaJAM v1.0.1
Metagenomic Pipeline for ancient DNA analysis performed at the Centre for Palaeogenetics - Stockholm

_Notes:_
Further developments will include:
- Addition of leeHom as an alternative choice to fastp
- Addition of parameters to define in fastp, sga, prinseq and bamdam

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

## Tools activation
Precise 1 or 0 for each step, 1=enable, 0=disable

## Parameters to precise for specific tools
1. *fastp*\
-overlap_len_require (default=20)\
-l (default=30)
2. *SGA*\
--dust-threshold (default=4)
3. *PRINSEQ*\
-lc_method (default=dust)\
-lc_threshold (default=4)\
-min_len (default=35)\
-derep (default=1)
4. *bamdam*\
--stranded (default=ds)\
--minreads (default=5)\
--maxdamage (default=0.5)

## SBATCH parameters
To be refine based on samples size and database using

## Overview of the pipeline
![alt text](https://github.com/NathanACO/metaJAM/blob/main/metaJAM_diagram.png)

