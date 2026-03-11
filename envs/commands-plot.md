

cd /scratch/gpfs/JENKINS/rl8728/PyPSA-China
module load anaconda3/2024.6
conda activate pypsa-china

git restore .
git pull
snakemake --unlock

sbatch job_fig_3_plot_value_scenario_comparison.slurm
sbatch job_plot_capacity.slurm

sbatch job_scenario_analysis.slurm
