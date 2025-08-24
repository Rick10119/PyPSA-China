l214365l@PU202502
module load anaconda3/2024.6
git clone https://github.com/Rick10119/PyPSA-China.git
conda env update -f envs/environment_plot.yaml
cd Documents/PyPSA-China
cd /scratch/gpfs/rl8728/PyPSA-China
conda activate pypsa-linux
mv ~/Documents/PyPSA-China /scratch/gpfs/rl8728/PyPSA-China-al
cp -R ~/Documents/PyPSA-China /scratch/gpfs/rl8728/PyPSA-China-1

conda deactivate
conda activate pypsa-plot

cd /scratch/gpfs/rl8728/PyPSA-China-0
module load anaconda3/2024.6
conda activate pypsa-plot

git restore .
git pull
snakemake --unlock
sbatch job_scenario_analysis.slurm
sbatch job_plot_capacity.slurm

sbatch jobs/job_HMM_2050_100p.slurm

cd /scratch/gpfs/rl8728/PyPSA-China
module load anaconda3/2024.6
conda activate pypsa-plot

git restore .
git pull
snakemake --unlock
chmod +x submit_multiple_jobs.sh 
./submit_multiple_jobs.sh


取消所有任务：
scancel -u rl8728
