#!/bin/bash -l
#SBATCH -A naiss2025-5-78
#SBATCH -p shared
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH -J launch_pipeline
#SBATCH -o /cfs/klemming/projects/supr/sllstore2017093/sediment/nathan/log/launch_pipeline.%j.out
#SBATCH -e /cfs/klemming/projects/supr/sllstore2017093/sediment/nathan/log/launch_pipeline.%j.err

bash -x /cfs/klemming/projects/supr/sllstore2017093/sediment/nathan/Run_pipeline_metage.sh /cfs/klemming/projects/supr/sllstore2017093/sediment/nathan/config_pipeline_metage.sh
#sbatch -x /cfs/klemming/projects/supr/sllstore2017093/sediment/nathan/Run_pipeline_metage.sh /cfs/klemming/projects/supr/sllstore2017093/sediment/nathan/config_pipeline_metage.sh
#Use bash -x to follow debuging with xtrace
