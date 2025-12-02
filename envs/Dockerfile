FROM registry.codeocean.com/codeocean/py-r:python3.12.8-R4.4.2-JupyterLab4.3.5-RStudiorstudio-server-2024.12.1-563-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

# Install system dependencies required for GIS libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    libpoppler-cpp-dev \
    libgdal-dev \
    libspatialindex-dev \
    && rm -rf /var/lib/apt/lists/*

# Configure conda channels
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels gurobi

# Install conda packages
RUN conda install -y \
    python>=3.9 \
    pip \
    pypsa=0.29.0 \
    atlite \
    dask \
    gurobi \
    xlrd \
    openpyxl \
    pycountry \
    seaborn \
    snakemake-minimal \
    memory_profiler \
    pyyaml \
    pytables \
    lxml \
    numpy=1.24.3 \
    pandas \
    geopandas \
    xarray=2024.2.0 \
    netcdf4 \
    networkx \
    scipy \
    shapely \
    progressbar2 \
    pyomo \
    matplotlib \
    tqdm \
    proj \
    proj-data \
    fiona \
    geopy \
    cartopy \
    descartes \
    rasterio \
    pytz \
    country_converter \
    tabula-py \
    graphviz \
    ipython \
    poppler \
    pyogrio \
    && conda clean -ya

# Set PROJ environment variable to ensure PROJ can find its database
ENV PROJ_LIB=/opt/conda/share/proj

# Install pip packages
RUN pip install --no-cache-dir \
    vresutils>=0.3.1 \
    tsam>=1.1.0 \
    powerplantmatching>=0.4.8
