#!/bin/bash
#SBATCH --job-name=bbh-population
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Submit from the checkout after activating the Python environment:
# POSTERIOR=/path/to/official.h5 OUTPUT_DIR=outputs/my-paper-run sbatch examples/generate.slurm
# Add a partition/account to sbatch as required by your cluster.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the repository root}"
: "${POSTERIOR:?set POSTERIOR to the official BP2P posterior file}"
: "${OUTPUT_DIR:?set OUTPUT_DIR to a new output directory}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg
python -m bbh_population generate \
  --config "${CONFIG:-configs/paper_append.json}" \
  --posterior "$POSTERIOR" --output-dir "$OUTPUT_DIR" \
  --workers "${SLURM_CPUS_PER_TASK:-1}"
python -m bbh_population validate --run-dir "$OUTPUT_DIR" \
  --workers "${SLURM_CPUS_PER_TASK:-1}"
