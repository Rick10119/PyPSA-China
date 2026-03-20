#!/bin/bash
#SBATCH --job-name=optimal-point    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=40        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=15G         # memory per cpu-core (4G is default)
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
#SBATCH --qos=high               # high priority queue
#SBATCH --qos=urgent             # urgent priority queue (if available)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=rl8728@princeton.edu


module purge
module load anaconda3/2024.10
conda activate pypsa-china

# Create log directory
mkdir -p logs

echo "=========================================="
echo "Starting optimal-point analysis job"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "=========================================="


echo "=========================================="
echo "Running optimal-point analysis script..."
echo "=========================================="

# Run optimal-point analysis script
python scripts/plot_optimal_point.py

if [ $? -eq 0 ]; then
    echo "✅ Optimal-point analysis script finished successfully"
else
    echo "❌ Optimal-point analysis script FAILED"
    exit 1
fi

echo "=========================================="
echo "Optimal-point analysis script finished!"
echo "End time: $(date)"
echo "Output directory: results/optimal_point_analysis"
echo "=========================================="

# Show generated files
echo "List of generated files:"
ls -la results/optimal_point_analysis/plots/ 2>/dev/null || echo "plots directory does not exist"
ls -la results/optimal_point_analysis/ 2>/dev/null || echo "optimal_point_analysis directory does not exist"

echo "=========================================="
echo "Job finished!"
echo "=========================================="
