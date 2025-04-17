conda install -c conda-forge pandas=1.4.0
conda install -c gurobi gurobi

# 安装命令
# conda env update -f envs/environment.yaml

# conda install -c conda-forge dask xlrd openpyxl pycountry seaborn
# conda install -c conda-forge memory_profiler yaml pytables lxml powerplantmatching pyyaml
# conda install -c conda-forge numpy pandas geopandas xarray netcdf4 networkx scipy shapely progressbar2 pyomo matplotlib
# conda install -c conda-forge proj 'fiona<=1.18.20' geopy tqdm pytz country_converter tabula-py
# conda install -c conda-forge ipython cartopy descartes rasterio
# conda install -c conda-forge -c bioconda snakemake-minimal
# conda install -c conda-forge graphviz
# pip install vresutils tsam
# conda install -c conda-forge xarray
# # 核心依赖
# conda install -c conda-forge pypsa=0.21.1 atlite=0.2.14 dask

snakemake solve_network_myopic --config opts=ll topology=current+Neighbor pathway=beta25175 planning_horizons=2020 heating_demand=positive

conda remove

conda env export --no-builds > environment.yml

snakemake --cores 6

./run_pipeline.sh

snakemake --dag | dot -Tpdf > dag.pdf

snakemake -s Snakefile-network --cores 6

