
module load anaconda3/2024.6
conda env update -f envs/environment.yaml
mv ~/Documents/PyPSA-China /scratch/gpfs/rl8728/PyPSA-China-al
cp -R ~/Documents/PyPSA-China /scratch/gpfs/rl8728/PyPSA-China-1



## scenario test
cd /scratch/gpfs/JENKINS/rl8728/PyPSA-China
module load anaconda3/2024.6
conda activate pypsa-china

git restore .
git pull
snakemake --unlock
chmod +x submit_multiple_jobs.sh 
./submit_multiple_jobs.sh

chmod +x submit_core_scenario.sh 
./submit_core_scenario.sh




Sep 2 run the jobs:

cd /scratch/gpfs/rl8728/PyPSA-China-1
module load anaconda3/2024.6
conda activate pypsa-china

git restore .
git pull
find ./results | xargs touch
snakemake --unlock
chmod +x submit_multiple_jobs.sh 
./submit_multiple_jobs.sh


取消所有任务：
scancel -u rl8728

更新所有结果文件：
find ./results | xargs touch


Sep 2, plot:

find ./results | xargs touch
snakemake --unlock
snakemake --configfile configs/config_MMMF_2050_10p.yaml -np --rerun-incomplete --ignore-incomplete --rerun-triggers mtime
snakemake --configfile configs/config_MMMF_2050_10p.yaml --cores 6 --rerun-incomplete --ignore-incomplete --rerun-triggers mtime
