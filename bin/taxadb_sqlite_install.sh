#!/bin/bash -l
#SBATCH -A naiss2025-5-185
#SBATCH --mem=880GB
#SBATCH -p memory
#SBATCH --cpus-per-task=256
#SBATCH -t 36:00:00
#SBATCH -J taxadb_sqlite_install
#SBATCH -e /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.err
#SBATCH -o /cfs/klemming/projects/snic/sediment_paleogenomics/nobackup/JAMIE_ALUMBAUGH/shotgun/scripts/slurm_output/slurm_%A-%a.out   

OUTPUT_DIR="${OUT_ROOT}"

#download NCBI taxonomy data
#-t full is too large to make into a database for Dardel
taxadb download -t nucl -o $OUTPUT_DIR
wait

#make database, include --fast flag otherwise it will take forever (>36 hours), need to specify database type with -d
taxadb create -i $OUTPUT_DIR --dbname $OUTPUT_DIR/taxadb_nucl.sqlite -d nucl --fast
wait 
