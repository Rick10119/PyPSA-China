# PyPSA-China: An Open Optimization Model of the Chinese Energy System

PyPSA-China is an open-source optimization model for the Chinese energy system built on the [PyPSA](https://pypsa.org/) framework. The model enables capacity expansion planning and operational optimization of China's power system, with special focus on aluminum smelting integration and grid flexibility.

## Features

### Core Capabilities

- **Multi-sector Energy System Modeling**: Integrated modeling of electricity, heat, gas, and other energy carriers
- **Renewable Energy Integration**: Detailed modeling of onshore/offshore wind and solar PV with resource potential
- **Transmission Network**: Provincial-level transmission network with multiple voltage levels
- **Storage Technologies**: Battery storage, pumped hydro storage, and hydrogen storage
- **Aluminum Smelting Integration**: Advanced iterative optimization algorithm for aluminum smelter operation
- **Scenario Analysis**: Multi-dimensional scenario framework with 27 scenario combinations
- **Visualization Tools**: Comprehensive plotting and analysis tools for results comparison

### Key Innovations

1. **Iterative Aluminum Optimization**: Novel algorithm that iteratively optimizes aluminum smelter operation based on nodal electricity prices
2. **Three-Dimensional Scenario Framework**: 
   - Smelter operational flexibility (low/mid/high)
   - Primary aluminum demand (low/mid/high)
   - Grid-interaction market opportunity (low/mid/high)
3. **Flexible Capacity Ratios**: Configurable aluminum smelter capacity ratios (100%, 90%, 80%, 70%, 60%)
4. **High-Performance Computing**: SLURM job management for large-scale scenario runs

## Installation

### Prerequisites

- Python 3.8+
- Gurobi solver (or other compatible solver)
- Sufficient memory (20-100 GB depending on network size)

### Environment Setup

1. Clone the repository:
```bash
git clone https://github.com/your-repo/PyPSA-China.git
cd PyPSA-China
```

2. Create and activate conda environment:
```bash
conda env create -f envs/environment_linux.yaml
conda activate pypsa-china
```

3. Install additional dependencies if needed:
```bash
pip install -r requirements.txt
```

## Quick Start

### Basic Usage

1. **Configure the model** by editing `config.yaml`:
   - Set planning horizon (default: 2050)
   - Configure scenario parameters
   - Adjust solver settings

2. **Prepare base network**:
```bash
snakemake -j 1 prepare_base_networks
```

3. **Run optimization**:
```bash
snakemake -j 1 solve_networks
```

4. **Generate summary**:
```bash
snakemake -j 1 make_summary
```

### Running with Aluminum Integration

To enable aluminum smelting integration:

1. Set in `config.yaml`:
```yaml
add_aluminum: True
iterative_optimization: True
aluminum_max_iterations: 10
aluminum_convergence_tolerance: 0.01
```

2. Run the optimization workflow as above.

## Documentation

Comprehensive documentation is available in the `docs/` folder:

### Main Documentation Files

- **[Aluminum Iterative Optimization Guide](docs/README_aluminum_iterative.md)**: Detailed explanation of the iterative optimization algorithm for aluminum smelters, including convergence criteria, network recreation methods, and virtual generator marginal cost settings.

- **[Scenario Dimensions Guide](docs/scenario_dimensions_guide.md)**: Complete guide to configuring and using the three-dimensional scenario framework (smelter flexibility, primary demand, grid interaction).

- **[Scenario Visualization Guide](docs/scenario_visualization_guide.md)**: Instructions for visualizing and comparing scenario results using `plot_scenario_comparison.py`.

- **[SLURM Jobs Guide](docs/slurm_jobs_guide.md)**: Guide for generating and managing SLURM job files for running multiple scenarios on HPC clusters.

## Configuration

### Key Configuration Parameters

#### Aluminum Settings
```yaml
add_aluminum: True                    # Enable aluminum integration
iterative_optimization: True          # Use iterative optimization
aluminum_max_iterations: 10           # Maximum iterations
aluminum_convergence_tolerance: 0.01  # Convergence threshold (1%)
aluminum_capacity_ratio: 1.0          # Capacity ratio (1.0 = 100%)
```

#### Scenario Dimensions
```yaml
aluminum:
  scenario_dimensions:
    smelter_flexibility: "mid"    # low, mid, high
    primary_demand: "mid"         # low, mid, high
    grid_interaction: "mid"       # low, mid, high
```

#### Solver Settings
```yaml
solving:
  solver:
    name: gurobi
  solver_options:
    default:
      Threads: 192
      Method: 2  # barrier method
```

## Project Structure

```
PyPSA-China/
├── config.yaml              # Main configuration file
├── configs/                 # Scenario-specific configs
├── data/                    # Input data
│   ├── costs/              # Technology costs
│   ├── grids/              # Grid topology
│   ├── load/               # Load profiles
│   ├── resources/          # Renewable resource data
│   └── ...
├── scripts/                # Python scripts
│   ├── prepare_base_network*.py
│   ├── solve_network*.py
│   ├── plot_*.py
│   └── ...
├── docs/                   # Documentation
├── results/                # Output results
├── Snakefile              # Snakemake workflow
└── jobs/                  # SLURM job files
```

## Running Scenarios

### Using Snakemake

Run a single scenario:
```bash
snakemake -j 1 solve_networks
```

Run with specific config:
```bash
snakemake -j 1 solve_networks --configfile configs/config_HHH_2050_100p.yaml
```

### Using SLURM (HPC)

1. Generate SLURM jobs:
```bash
python scripts/generate_slurm_jobs_advanced.py
```

2. Submit jobs:
```bash
./submit_multiple_jobs.sh
```

3. Monitor jobs:
```bash
squeue -u $USER
```

## Scenario Analysis

### Generating Scenario Configurations

The model supports generating multiple scenario combinations:
- 3 smelter flexibility levels × 3 demand levels × 3 market opportunity levels = 27 combinations
- Each combination can be run with different capacity ratios

### Visualizing Results

Use the scenario comparison tool:
```bash
python scripts/plot_scenario_comparison.py --file-type costs
python scripts/plot_scenario_comparison.py --file-type capacities
```

This generates:
- 9-panel comparison plots (one per demand-market combination)
- Summary tables with key metrics
- Detailed data files for further analysis

## Key Algorithms

### Iterative Aluminum Optimization

The iterative optimization algorithm:

1. **Initialization**: Start with empty aluminum usage pattern
2. **Iteration Loop**:
   - Step 1: Solve network with continuous aluminum model → get nodal prices
   - Step 2: Solve aluminum optimal operation problem based on nodal prices
   - Step 3: Check objective function change for convergence
3. **Convergence**: Stop when relative objective change < threshold
4. **Output**: Final network results and timing statistics

Key improvements:
- Convergence based on objective function change (not aluminum usage change)
- Network recreation at each iteration for clean state
- Virtual generators use nodal marginal prices
- Fixed aluminum usage via `p_set` (not constraints)

## Output Files

Results are organized by version and scenario:
```
results/
└── version-<version>-<scenario>/
    ├── networks/
    ├── summary/
    │   └── postnetworks/
    │       └── costs.csv
    │       └── capacities.csv
    └── ...
```

## Troubleshooting

### Common Issues

1. **Solver errors**: Check Gurobi license and solver options
2. **Memory issues**: Reduce network size or increase `mem_per_thread`
3. **Convergence issues**: Adjust `aluminum_convergence_tolerance` or `aluminum_max_iterations`
4. **Missing data**: Verify all required data files are present in `data/` directory

### Getting Help

- Check the documentation in `docs/`
- Review example configurations in `configs/`
- Check log files in `logs/`

## License

This project is licensed under multiple licenses:
- Code: MIT License (see `LICENSES/MIT.txt`)
- Data: CC0-1.0 (see `LICENSES/CC0-1.0.txt`)
- Documentation: CC-BY-4.0 (see `LICENSES/CC-BY-4.0.txt`)

## Citation

If you use PyPSA-China in your research, please cite:

```bibtex
@software{pypsa_china,
  title = {PyPSA-China: An Open Optimization Model of the Chinese Energy System},
  author = {PyPSA-China Authors},
  year = {2022},
  url = {https://github.com/your-repo/PyPSA-China}
}
```

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## Acknowledgments

This project builds on the [PyPSA](https://pypsa.org/) framework and is inspired by [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur).

## Contact

For questions and support, please open an issue on GitHub or contact the maintainers.
