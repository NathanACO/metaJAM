#!/bin/bash -l
#SBATCH -A naiss2025-5-172
#SBATCH -p shared
#SBATCH --cpus-per-task=16
#SBATCH --time 08:00:00
#SBATCH --job-name=mmseqs2
#SBATCH -e /cfs/klemming/projects/supr/naiss2025-23-301/output/error/%x-%j.error
#SBATCH -o /cfs/klemming/projects/supr/naiss2025-23-301/output/out/%x-%j.out

set -euxo pipefail
echo "MMSeqs2 stub — insert your exact command lines here (defaults or your params)."