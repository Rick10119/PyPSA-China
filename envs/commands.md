conda install -c conda-forge pandas=1.4.0
conda install -c gurobi gurobi

# 安装命令
# conda env update -f envs/environment.yaml

# # 核心依赖
# conda install -c conda-forge pypsa=0.21.3 atlite=0.2.14 dask

snakemake solve_network_myopic --config opts=ll topology=current+Neighbor pathway=beta25175 planning_horizons=2020 heating_demand=positive

conda remove

conda env export --no-builds > environment.yaml

conda create --name cloned_env --clone original_env

snakemake --cores 6

./run_pipeline.sh

snakemake --dag | dot -Tpdf > dag.pdf

snakemake -s Snakefile-network --cores 6

snakemake --unlock
snakemake --cores 6

module load anaconda3/2024.6
git clone https://github.com/Rick10119/PyPSA-China.git
conda env update -f envs/environment_mac.yaml
cd Documents/github/PyPSA-China
conda activate pypsa-linux
