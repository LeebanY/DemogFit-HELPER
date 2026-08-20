#!/bin/bash
#SBATCH --job-name=demog_pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=<your_partition>

# --- Before ever submitting this script, run ONCE on the login node -------
#   snakemake --use-conda --conda-frontend mamba --conda-create-envs-only --cores 1
# so all tool environments are already cached and no internet access is
# needed on the compute node.

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate demography_pipeline

snakemake \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --use-conda \
    --conda-frontend mamba \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds
