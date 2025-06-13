l214365l@PU202502
module load anaconda3/2024.6
git clone https://github.com/Rick10119/PyPSA-China.git
conda env update -f envs/environment_mac.yaml
cd Documents/PyPSA-China
cd /scratch/gpfs/rl8728/PyPSA-China
conda activate pypsa-linux
mv ~/Documents/PyPSA-China /scratch/gpfs/rl8728/PyPSA-China-al
cp -R ~/Documents/PyPSA-China /scratch/gpfs/rl8728/PyPSA-China-1

snakemake --unlock
git pull
sbatch job.slurm

cd /scratch/gpfs/rl8728/PyPSA-China
module load anaconda3/2024.6
conda activate pypsa-linux
