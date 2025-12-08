#!/bin/bash
#SBATCH --job-name=sim_array
#SBATCH --output=/scratch/richards/guillaume.butler-laporte/non_linear_mr/simulations/logs/sim_%A_%a.out
#SBATCH --error=/scratch/richards/guillaume.butler-laporte/non_linear_mr/simulations/logs/sim_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-23

source ~/.bashrc
conda activate r4-env

cd /scratch/richards/guillaume.butler-laporte/non_linear_mr/simulations/

# Create logs directory if it doesn't exist
mkdir -p results

# Define parameter arrays
INNER_INTERVALS=0.90
N_REPLICATES=50
PROP_INTS=(0 0.1 0.25 0.5)
INST_INTS=(30 50)
SS_INTS=(10000 25000 50000)

# Calculate total number of combinations
N_PROP=${#PROP_INTS[@]}
N_INST=${#INST_INTS[@]}
N_SS=${#SS_INTS[@]}

TOTAL_JOBS=$((N_PROP * N_INST * N_SS))

echo "Total job combinations:" $TOTAL_JOBS
echo "Current array task ID:" $SLURM_ARRAY_TASK_ID

# Make new arrays:

new_prop=()
new_inst=()
new_ss=()

for v1 in "${PROP_INTS[@]}"; do
  for v2 in "${INST_INTS[@]}"; do
    for v3 in "${SS_INTS[@]}"; do
      new_prop+=("$v1")
      new_inst+=("$v2")
      new_ss+=("$v3")
    done
  done
done


# Create output filename
OUTPUT_NAME="sim_inner${INNER_INTERVALS}_nrep${N_REPLICATES}_prop_pleio_${new_prop[$SLURM_ARRAY_TASK_ID]}_inst${new_inst[$SLURM_ARRAY_TASK_ID]}_ss${new_ss[$SLURM_ARRAY_TASK_ID]}.rds"

echo "=========================================="
echo "Job Configuration:"
echo "  Inner interval: $INNER_INTERVALS"
echo "  N replicates: $N_REPLICATES"
echo "  Prop pleiotropy: ${new_prop[$SLURM_ARRAY_TASK_ID]}"
echo "  Instruments: ${new_inst[$SLURM_ARRAY_TASK_ID]}"
echo "  Sample size: ${new_ss[$SLURM_ARRAY_TASK_ID]}"
echo "  Output: $OUTPUT_NAME"
echo "=========================================="

# Run the simulation
Rscript /scratch/richards/guillaume.butler-laporte/non_linear_mr/simulations/rscript_simulations.R \
  --inner_interval $INNER_INTERVALS \
  --n_replicates $N_REPLICATES \
  --prop_int "${new_prop[$SLURM_ARRAY_TASK_ID]}" \
  --inst_int "${new_inst[$SLURM_ARRAY_TASK_ID]}" \
  --ss_int "${new_ss[$SLURM_ARRAY_TASK_ID]}" \
  --output_dir ./results \
  --output_name "$OUTPUT_NAME"

echo "Job completed successfully!"