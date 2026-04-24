#!/bin/bash -l
#SBATCH -A project
#SBATCH -p partition
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -J launch_pipeline
#SBATCH -o launch_pipeline.%j.out
#SBATCH -e launch_pipeline.%j.err

bash -x Run_pipeline_metage.sh config_pipeline_metage.sh
