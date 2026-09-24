#!/bin/bash

# ==============================================================================
# 1. APPTAINER CACHE SETUP (Pointing to your existing RNA-seq cache!)
# ==============================================================================
SHARED_CACHE="/hpc/umc_kaaij/gkarsters/Snakemake/Pipelines/RNA-seq/containers/cache"
mkdir -p "$SHARED_CACHE"

export APPTAINER_CACHEDIR="$SHARED_CACHE"
export SINGULARITY_CACHEDIR="$SHARED_CACHE"

# Local tmp for image extraction (prevents /tmp from filling up)
mkdir -p "$PWD/containers/tmp"
export APPTAINER_TMPDIR="$PWD/containers/tmp"
export SINGULARITY_TMPDIR="$PWD/containers/tmp"

# ==============================================================================
# 2. VARIABLES & CLI PARSING
# ==============================================================================
DRY_RUN=false
TOUCH_TARGET=""
UNLOCK=false

# Function to show help
show_help() {
    echo "Usage: ./run_pipeline.sh [options]"
    echo "Options:"
    echo "  -n          Dry-run (print what would happen, don't execute)"
    echo "  -t STRING   Touch mode. Marks specific rules as 'done'."
    echo "  -u          Unlock directory (if Snakemake crashed previously)"
    echo "  -h          Show this help message"
    exit 1
}

while getopts "t:nuh" opt; do
    case ${opt} in
        n) DRY_RUN=true ;;
        t) TOUCH_TARGET=$OPTARG ;;
        u) UNLOCK=true ;;
        h) show_help ;;
        *) show_help ;;
    esac
done

# ==============================================================================
# 3. ENVIRONMENT & SETUP
# ==============================================================================
source $(conda info --base)/etc/profile.d/conda.sh
conda activate snakemake-c

LOG_DIR="results/logs/slurm_logs"
mkdir -p "$LOG_DIR"

CMD=(
    snakemake
    -s Snakefile
    --use-singularity
    --singularity-prefix "/hpc/umc_kaaij/gkarsters/Snakemake/Pipelines/RNA-seq/.snakemake/singularity"
    --singularity-args "--cleanenv --bind /hpc"
    --rerun-incomplete
    --printshellcmds
    --restart-times 3
    --latency-wait 60
    --default-resources "runtime=60" "mem_mb=8000"
    --jobs 500                   # <--- Forces 500 jobs max, even on dry-run!
    --resources mem_mb=200000    # <--- Prevents local memory scaling!
)

# ==============================================================================
# 4. EXECUTION MODES
# ==============================================================================

# Unlock
if [ "$UNLOCK" = true ]; then
    echo "Unlocking directory..."
    "${CMD[@]}" --unlock
    exit 0
fi

# Touch Mode
if [ ! -z "$TOUCH_TARGET" ]; then
    echo "Touch Mode selected for rule(s): '$TOUCH_TARGET'"
    "${CMD[@]}" --touch --cores 1 -R $TOUCH_TARGET
    exit 0
fi

# SLURM execution
CLUSTER_CMD="sbatch \
    --parsable \
    --partition=cpu \
    --time={resources.runtime} \
    --mem={resources.mem_mb}M \
    --cpus-per-task={threads} \
    --output=$LOG_DIR/slurm-%j.out \
    --error=$LOG_DIR/slurm-%j.err \
    --mail-type=FAIL \
    --mail-user=g.j.karsters@umcutrecht.nl \
    --gres=tmpspace:170G"

if [ "$DRY_RUN" = true ]; then
    echo "Performing Dry-Run..."
    "${CMD[@]}" -npr
else
    echo "Submitting to SLURM..."
    # Removed the redundant --jobs 500 here since it is now in CMD
    "${CMD[@]}" --cluster "$CLUSTER_CMD" --cluster-status "$PWD/slurm_status.sh" --cluster-cancel "scancel"
fi