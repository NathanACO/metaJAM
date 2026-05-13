<img src="https://github.com/NathanACO/metaJAM/blob/main/metaJAM_logo.png" width="200" /> 

# metaJAM v1.1.1 nf version

Metagenomic Pipeline for ancient DNA analysis performed at the Centre for Palaeogenetics, Stockholm University.

This pipeline performs processing and analysis of metagenomic data, starting from paired-end fastq files or any intermediate files used in the pipeline if not the entire workflow is desired to run.


> [!NOTE]
> Further developments will be included in metaJAM v1.2.2 in mid 2026:
> - Microbial version of the pipeline to investigate the microbial signal, notably using a parallel competitive mapping against the GTDB database



## Overview of the pipeline
![alt text](https://github.com/NathanACO/metaJAM/blob/metaJAM-nf/metaJAM_diagram_v1.1.1.png)

## Required program
nextflow (developed in v25.10.3)\
conda\
python3\
pip3\
git\ #optional, as for downloading the repository from github with git. No need for git if you download zip of this repository directly
taxadb #this can be downloaded with pip3 with the tutorial script below

## Download this pipeline
run in your Linux terminal
```
git clone https://github.com/NathanACO/metaJAM.git
cd metaJAM
```

## setup script to run before running Nextflow

Run the tutorial script once before running anything to set up the taxonomy files (taxdump.tar.gz names.dmp nodes.dmp), taxadb.sqlite (for MMseqs2) and placeholders for the flexibility of the pipeline (NO_FILE*) in the assets/ directory.
```
cd assets
#prepare NCBI taxonomy information: names.dmp and nodes.dmp
curl -O https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz names.dmp nodes.dmp
#prepare MMseq2 taxonomy database
python -m venv taxadb
source taxadb/bin/activate
pip3 install taxadb
taxadb download -o taxadb -t nucl
taxadb create -i taxadb --dbname taxadb.sqlite -d nucl --fast
cd ../ 
```

## A thorough explanation of the rationale of each process of the workflow is in [wiki](https://github.com/NathanACO/metaJAM/wiki)
## To specify input and parameter setting in `nextflow.config`, please refer to [input readme](https://github.com/NathanACO/metaJAM/blob/ae8b5cfbb4d81ccd059c7ae82d2ea0aeb371e042/README.input_nextflow.config.md)


## How to launch it:
After running the setup script and having your nextflow.config (or other name) ready, you could run\
`nextflow run main.nf -profile conda -c other_name_nextflow.config`\
If you want to resume, add also `-resume`, and `--with-trace` for output a trace*.txt reporting memory and time for each process. You can omit the  `-c nextflow.config`, if you are running exactly nextflow.config. Else if your nextflow.config is named in another way, you need to specify with the config name `-c other_name_nextflow.config`.

> [!TIP]
> You can run it with the Test samples present in the test folder of this github.

## Plot output
Different plots are created by metaJAM and require a metadata file (see in test for metadata format and content)[In development].\
The first plot created gives an overview of the metrics of the samples processed:
![alt text](https://github.com/NathanACO/metaJAM/blob/main/reads_per_step_dots.png)
&copy; Martin et al. 2025 plot, peatbog paper - in revision

The second plot gives the taxa composition as a bubble plot or a heatmap (can be chosen in the config file).\
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
