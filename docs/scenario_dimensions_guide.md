# Scenario Dimensions Configuration Guide

This document describes how to configure and use the three scenario dimensions in the PyPSA-China project: smelter operational flexibility, primary aluminum demand, and grid-interaction market opportunity.

## Overview

The project supports three independent scenario dimensions, each with three levels (low, mid, high), yielding up to 27 distinct scenario combinations.

### The Three Dimensions

1. **Smelter Operational Flexibility**
   - Controls operational parameters of the aluminum smelter.
   - Includes overcapacity rate, minimum power, restart allowance, and restart cost.

2. **Primary Aluminum Demand**
   - Controls the structure of aluminum demand.
   - Includes domestic demand ratio, export rate, recycling rate, and product lifetime.

3. **Grid-Interaction Market Opportunity**
   - Controls cost parameters related to grid interaction.
   - Includes VRE, battery, and H₂ storage cost reductions, as well as other flexible demand.

## Configuration

### Structure in `config.yaml`

```yaml
aluminum:
  # Scenario dimensions
  scenario_dimensions:
    # 1. Smelter operational flexibility
    smelter_flexibility:
      low:
        overcapacity_rate: 0.9  # p.u.
        p_min: 0.9  # p.u.
        allow_restart: false
        restart_costs: 110000  # $/MW
      mid:
        overcapacity_rate: 0.9  # p.u.
        p_min: 0.9  # p.u.
        allow_restart: true
        restart_costs: 16000  # $/MW
      high:
        overcapacity_rate: 0.7  # p.u.
        p_min: 0.7  # p.u.
        allow_restart: true
        restart_costs: 3200  # $/MW

    # 2. Primary aluminum demand
    primary_demand:
      low:
        domestic_demand_ratio: 0.8  # 80 %
        export_rate: 0.2  # 20 %
        recycling_rate: 0.2  # 20 %
        product_lifetime: 20  # years
      mid:
        domestic_demand_ratio: 0.7  # 70 %
        export_rate: 0.3  # 30 %
        recycling_rate: 0.16  # 16 %
        product_lifetime: 16  # years
      high:
        domestic_demand_ratio: 0.6  # 60 %
        export_rate: 0.4  # 40 %
        recycling_rate: 0.12  # 12 %
        product_lifetime: 12  # years

    # 3. Grid-interaction market opportunity
    grid_interaction:
      low:
        vre_cost_reduction: 0.0
        battery_cost_reduction: 0.0
        h2_storage_cost_reduction: 0.0
        other_flexible_demand: 0.0
      mid:
        vre_cost_reduction: 0.1  # 10 %
        battery_cost_reduction: 0.1  # 10 %
        h2_storage_cost_reduction: 0.1  # 10 %
        other_flexible_demand: 0.05  # 5 %
      high:
        vre_cost_reduction: 0.2  # 20 %
        battery_cost_reduction: 0.2  # 20 %
        h2_storage_cost_reduction: 0.2  # 20 %
        other_flexible_demand: 0.1  # 10 %

  # Currently selected scenario combination (default: mid-mid-mid)
  current_scenario:
    smelter_flexibility: "mid"
    primary_demand: "mid"
    grid_interaction: "mid"
```

## Usage

### 1. Basic usage

```python
from scripts.scenario_utils import load_config, get_scenario_params

# Load configuration
config = load_config()

# Get default scenario parameters
default_params = get_scenario_params(config)

# Get parameters for a specific scenario
high_flex_params = get_scenario_params(
    config,
    smelter_flexibility="high",
    primary_demand="mid",
    grid_interaction="low"
)
```

### 2. Retrieve parameters for a single dimension

```python
from scripts.scenario_utils import (
    get_smelter_params,
    get_demand_params,
    get_grid_interaction_params
)

# Smelter parameters
smelter_params = get_smelter_params(config, "high")

# Demand parameters
demand_params = get_demand_params(config, "low")

# Grid-interaction parameters
grid_params = get_grid_interaction_params(config, "mid")
```

### 3. Generate all scenario combinations

```python
from scripts.scenario_utils import generate_scenario_combinations

# Generate all 27 combinations
combinations = generate_scenario_combinations()

for combo in combinations:
    print(f"Scenario: {combo['name']}")
    print(f"  Smelter flexibility: {combo['smelter_flexibility']}")
    print(f"  Primary demand: {combo['primary_demand']}")
    print(f"  Grid interaction: {combo['grid_interaction']}")
```

### 4. Using scenario parameters in a PyPSA network

```python
# Get scenario parameters
params = get_scenario_params(config, "high", "mid", "low")

# Apply smelter parameters
smelter = params['smelter_flexibility']
network.generators.loc[aluminum_generators, 'p_min_pu'] = smelter['p_min']
network.generators.loc[aluminum_generators, 'start_up_cost'] = smelter['restart_costs']

# Apply demand parameters
demand = params['primary_demand']
aluminum_demand = base_demand * demand['domestic_demand_ratio']

# Apply grid-interaction parameters
grid = params['grid_interaction']
vre_cost_adjusted = base_vre_cost * (1 - grid['vre_cost_reduction'])
```

## Representative Scenario Combinations

1. **Conservative (low-low-low)**
   - Low smelter flexibility
   - Low primary demand
   - Low grid-interaction opportunity

2. **Baseline (mid-mid-mid)**
   - Moderate smelter flexibility
   - Moderate primary demand
   - Moderate grid-interaction opportunity

3. **Aggressive (high-high-high)**
   - High smelter flexibility
   - High primary demand
   - High grid-interaction opportunity

4. **High flexibility – Low demand – High grid interaction (high-low-high)**
   - Suitable for studying aluminum smelters as grid flexibility resources.

5. **Low flexibility – High demand – Low grid interaction (low-high-low)**
   - Suitable for studying traditional smelter operation patterns.

## Parameter Descriptions

### Smelter Operational Flexibility

- **overcapacity_rate**: fraction above nameplate capacity the smelter can operate at.
- **p_min**: minimum power ratio when the smelter is running.
- **allow_restart**: whether the smelter is allowed to restart after a shutdown.
- **restart_costs**: restart cost in $/MW.

### Primary Aluminum Demand

- **domestic_demand_ratio**: domestic demand as a share of total demand.
- **export_rate**: exports as a share of total demand.
- **recycling_rate**: recycled aluminum as a share of total demand.
- **product_lifetime**: average product lifetime in years.

### Grid-Interaction Parameters

- **vre_cost_reduction**: fractional reduction in VRE (variable renewable energy) cost.
- **battery_cost_reduction**: fractional reduction in battery cost.
- **h2_storage_cost_reduction**: fractional reduction in H₂ storage cost.
- **other_flexible_demand**: other flexible demand as a share of total demand.

## Running Examples

Run the example script to see all features:

```bash
python examples/scenario_example.py
```

Run the utility function tests:

```bash
python scripts/scenario_utils.py
```

## Notes

1. **Defaults**: the current default scenario is `mid-mid-mid`.
2. **Validation**: all parameters fall within reasonable ranges, but verification before use is recommended.
3. **Backward compatibility**: the new scenario-dimension configuration does not affect existing functionality.
4. **Extensibility**: new levels or dimensions can be added easily.

## Troubleshooting

If you encounter issues, check the following:

1. Config file format is correct.
2. Scenario names are spelled correctly (`low`, `mid`, `high`).
3. Parameter values are within reasonable ranges.
4. Required packages are installed (e.g., `yaml`).
