#!/bin/bash

# Initialize conda for the current shell
source ~/.bash_profile

# Activate the first environment for network preparation
echo "Activating pypsa-china environment for network preparation..."
source activate pypsa-china

# Run the network preparation step
echo "Running network preparation..."
snakemake --cores 6

# Deactivate the first environment
source deactivate

# Activate the second environment for solving
echo "Activating pypsa-gurobi environment for network solving..."
source activate pypsa-gurobi

# Run the network solving step
echo "Running network solving..."
snakemake --cores 6

# Deactivate the second environment
source deactivate

echo "Pipeline completed!" 