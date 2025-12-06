#!/bin/bash
#SBATCH --job-name=sim_array
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-89

source activate myother_env

# Create logs directory if it doesn't exist
mkdir -p logs
mkdir -p results

# Define parameter arrays
INNER_INTERVALS=(0.90)
N_REPLICATES=(50)
PROP_INTS=("0,0.1,0.25,0.5")
INST_INTS=("30" "50")
SS_INTS=("10000" "25000" "50000")

# Calculate total number of combinations
N_INNER=${#INNER_INTERVALS[@]}
N_REPS=${#N_REPLICATES[@]}
N_PROP=${#PROP_INTS[@]}
N_INST=${#INST_INTS[@]}
N_SS=${#SS_INTS[@]}

TOTAL_JOBS=$((N_INNER * N_REPS * N_PROP * N_INST * N_SS))

echo "Total job combinations: $TOTAL_JOBS"
echo "Current array task ID: $SLURM_ARRAY_TASK_ID"

# Calculate indices for each parameter
idx=$SLURM_ARRAY_TASK_ID

ss_idx=$((idx % N_SS))
idx=$((idx / N_SS))

inst_idx=$((idx % N_INST))
idx=$((idx / N_INST))

prop_idx=$((idx % N_PROP))
idx=$((idx / N_PROP))

rep_idx=$((idx % N_REPS))
idx=$((idx / N_REPS))

inner_idx=$((idx % N_INNER))

# Get parameter values
INNER_INTERVAL=${INNER_INTERVALS[$inner_idx]}
N_REP=${N_REPLICATES[$rep_idx]}
PROP_INT=${PROP_INTS[$prop_idx]}
INST_INT=${INST_INTS[$inst_idx]}
SS_INT=${SS_INTS[$ss_idx]}

# Create output filename
OUTPUT_NAME="sim_inner${INNER_INTERVAL}_nrep${N_REP}_inst${INST_INT}_ss${SS_INT}.rds"

echo "=========================================="
echo "Job Configuration:"
echo "  Inner interval: $INNER_INTERVAL"
echo "  N replicates: $N_REP"
echo "  Prop int: $PROP_INT"
echo "  Instruments: $INST_INT"
echo "  Sample size: $SS_INT"
echo "  Output: $OUTPUT_NAME"
echo "=========================================="

# Load R module (adjust based on your cluster)
module load r/4.3.0  # Adjust this to match your cluster's R module

# Run the simulation
Rscript run_simulation_cli.R \
  --inner_interval $INNER_INTERVAL \
  --n_replicates $N_REP \
  --prop_int "$PROP_INT" \
  --inst_int "$INST_INT" \
  --ss_int "$SS_INT" \
  --output_dir ./results \
  --output_name "$OUTPUT_NAME"

echo "Job completed successfully!"
