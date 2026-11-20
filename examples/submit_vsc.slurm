#!/bin/bash
# Head job: Nextflow itself runs here and submits every step as its own Slurm job.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --account=lp_big_mindwell_cpu
#SBATCH --partition=dedicated_big_bigmem
#SBATCH --clusters=mindwell
#SBATCH --job-name=hybridphase_nf
#SBATCH --output=%x_%j.out
set -euo pipefail

source /vsc-hard-mounts/leuven-data/359/vsc35908/micromamba/etc/profile.d/micromamba.sh
micromamba activate nextflow_env          # any env with nextflow >= 23.10 and java 17
micromamba activate --stack singularity-env

cd "${SLURM_SUBMIT_DIR}"
# all intermediate files go here; override with WORKDIR=/some/big/place sbatch ...
export NXF_WORK="${WORKDIR:-${VSC_SCRATCH}/nf_work/hybridphase}"
mkdir -p "$NXF_WORK"
# sbatch AND squeue must look at the mindwell cluster, or Nextflow loses track of its jobs
export SLURM_CLUSTERS=mindwell
export NXF_OPTS="-Xms1g -Xmx6g"

PARAMS=${PARAMS:-examples/colo829bl.params.yaml}
nextflow run main.nf -profile vsc_mindwell -params-file "$PARAMS" -resume "$@"
