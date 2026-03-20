#!/bin/bash
#SBATCH --job-name=scenario-analysis    # create a short name for your job
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
echo "Starting scenario analysis job"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "=========================================="

# Create output directory
mkdir -p results/scenario_analysis

echo "=========================================="
echo "Running scenario cost-change analysis script..."
echo "=========================================="

# Run scenario cost-change analysis script
python scripts/plot_value_scenario_comparison.py \
    --config config.yaml \
    --output results/scenario_analysis \
    --file-type costs

python scripts/plot_value_scenario_comparison_f.py \
    --config config.yaml \
    --output results/scenario_analysis_f \
    --file-type costs

if [ $? -eq 0 ]; then
    echo "✅ Scenario cost-change analysis script finished successfully"
else
    echo "❌ Scenario cost-change analysis script FAILED"
    exit 1
fi

echo "=========================================="
echo "Running scenario capacity-change analysis script..."
echo "=========================================="

# Run scenario capacity-change analysis script
python scripts/plot_value_scenario_comparison.py \
    --config config.yaml \
    --output results/scenario_analysis \
    --file-type capacities

if [ $? -eq 0 ]; then
    echo "✅ Scenario capacity-change analysis script finished successfully"
else
    echo "❌ Scenario capacity-change analysis script FAILED"
    exit 1
fi

echo "=========================================="
echo "All scenario-analysis scripts finished!"
echo "End time: $(date)"
echo "Output directory: results/scenario_analysis"
echo "=========================================="

# Show generated files
echo "List of generated files:"
echo "Scenario comparison plots:"
ls -la results/scenario_analysis/scenario_plots/ 2>/dev/null || echo "scenario_plots directory does not exist"
echo ""
echo "Summary tables:"
ls -la results/scenario_analysis/summary_tables/ 2>/dev/null || echo "summary_tables directory does not exist"
echo ""
echo "Full directory tree:"
ls -la results/scenario_analysis/ 2>/dev/null || echo "scenario_analysis directory does not exist"

echo "=========================================="
echo "Job finished!"
echo "=========================================="
