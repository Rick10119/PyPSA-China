# Aluminum Smelter Iterative Optimization – Refactoring Notes

> **Potline-based modeling in one sentence:** Each province's aluminum smelter is represented by a single 250 kt/yr potline in a MILP sub-problem; the resulting commitment schedule is scaled back to the provincial total and fixed in the master problem via `p_set`, forming a closed-loop iteration.

## Refactoring Overview

`solve_network_aluminum_iterative.py` has been refactored following the framework of `capacity_expansion_planning_aluminum_iterative.py`. The key improvements are listed below.

## Key Improvements

### 1. Iteration Logic

- **Convergence criterion**: switched from aluminum power-consumption change to relative change in the objective function.
- **Convergence threshold**: defaults to 1 % (0.01), which is more robust than the previous energy-change metric.
- **Timing statistics**: detailed per-iteration timing is now recorded, including total, mean, fastest, and slowest iteration times.

### 2. Network Re-creation ⭐

- **Reload from file**: following `capacity_expansion_planning_aluminum_iterative.py`, the network is reloaded from disk at every iteration instead of deep-copied.
- **Clean state**: avoids accumulated solver artifacts in the network object.
- **More stable solves**: each iteration starts from a pristine network.

### 3. Fixing Aluminum Power Consumption ⭐

- **`p_set` approach**: aluminum power consumption is now fixed via `p_set` rather than explicit equality constraints, following `capacity_expansion_planning_aluminum_iterative.py`.
- **Simultaneous load fixing**: both the smelter link output and the corresponding aluminum load are fixed together.
- **Simpler code**: the complex constraint-injection logic has been removed.

### 4. Virtual Generator Marginal Cost ⭐

- **Nodal marginal prices**: the marginal cost of virtual generators is set to the nodal electricity price instead of zero.
- **Bus mapping**: virtual generators are correctly associated with the corresponding electricity bus.
- **Higher accuracy**: aluminum optimization is driven by realistic nodal prices.

### 5. Aluminum Optimization Network Simplification

- **Component filtering**: simplified the logic for identifying and retaining aluminum-related components.
- **Virtual generators**: improved the way virtual generators are added to electricity buses.
- **Carrier management**: ensures the virtual carrier exists to prevent runtime errors.

### 6. Optimization Method Simplification

- **Removed `safe_optimize`**: the complex wrapper has been replaced by PyPSA's standard `optimize()` call.
- **Simplified error handling**: reduced code complexity for error and status checks.

### 7. Configuration Support

- **Iteration parameters**: maximum iterations and convergence threshold can be set in the config file.
- **Feature toggle**: the `add_aluminum` flag enables or disables the aluminum iterative optimization.

## Configuration Parameters

The following parameters can be set in the config file:

```yaml
# Aluminum feature toggle
add_aluminum: true

# Iteration parameters
aluminum_max_iterations: 10          # Maximum number of iterations
aluminum_convergence_tolerance: 0.01 # Convergence threshold (1 %)

# Unit-commitment toggle
aluminum_commitment: false           # Enable unit-commitment constraints
```

## Main Functions

### `solve_network_iterative()`
- **Purpose**: main entry point for the aluminum iterative optimization.
- **Convergence**: stops when the relative change in the objective falls below the threshold.
- **Network handling**: reloads the network at each iteration to keep state clean.
- **Output**: timing statistics and the final solved network.

### `solve_aluminum_optimization()`
- **Purpose**: solves the aluminum optimal-dispatch sub-problem for a single province.
- **Input**: original network, config, solver settings.
- **Output**: aluminum power-consumption profile.

### `extra_functionality()`
- **Purpose**: injects supplementary constraints (CHP, transmission, retrofit, etc.).
- **Note**: aluminum power-consumption constraints are no longer added here because `p_set` is used instead.

## Usage

1. **Enable aluminum**: set `add_aluminum: true` in the config file.
2. **Set iteration parameters**: adjust `aluminum_max_iterations` and `aluminum_convergence_tolerance` as needed.
3. **Run**: use the standard PyPSA-China workflow.

## Algorithm Flow

1. **Initialization**: aluminum power-consumption profile and objective value are set to empty/null.
2. **Iterative solve**:
   - Step 1: Reload the network and solve the relaxed (continuous) aluminum model to obtain nodal prices.
   - Step 2: Based on nodal prices, solve the aluminum optimal-dispatch MILP sub-problem.
   - Step 3: Check whether the objective has converged.
3. **Convergence check**: stop when the relative change in the objective < threshold.
4. **Output**: return the final network and timing statistics.

## Potline-based Modeling (Representative Potline Method)

The current implementation converts each province's aluminum smelter into a "representative potline + scaling factor" form inside the sub-problem, reducing the MILP size while preserving the effects of start-up/shut-down costs and minimum-output constraints. Specifically, `solve_aluminum_optimization()` reads the province's annual production (in 10 kt/yr) and uses **250 kt/yr (25.0 in 10 kt/yr units)** as the standard potline size: if the provincial capacity does not exceed this value, a single potline of equal size is used directly; otherwise a single 250 kt/yr representative potline is created and `scale_factor = provincial_production / 25` is recorded. The potline capacity is converted to MW via `line_cap_mw = line_cap_10kt * 10000 * 13.3 / 8760`, and operational parameters (`p_min_pu`, `stand_by_cost`, `start_up_cost`) are recomputed by `get_aluminum_smelter_operational_params()` for the representative potline capacity, with `committable=True` enforced to explicitly model unit commitment.

At the solve level, the sub-problem strips all components outside the target province, keeping only the provincial aluminum link, load, and store. A virtual generator is added at the electricity bus with its marginal cost set to the nodal price extracted from the master problem (`nodal_prices`), so that the representative potline optimizes commitment and dispatch against exogenous electricity prices. After the solve, the time series is read from `n.links_t.p0` and multiplied back by `scale_factor` to recover the provincial total; if the solver returns infeasible or hits a time limit without a feasible solution, a flat load curve based on annual production is used as a fallback. Finally, in the master iteration the returned provincial aluminum consumption is written directly into `links_t.p_set` and `loads_t.p_set` (no equality constraints), completing the "potline sub-problem optimization → master problem fixing" closed-loop iteration.

## Network Re-creation

### Old approach (deep copy)
```python
# Deep-copy the network object
n_current = copy.deepcopy(n)
n_current.config = config
n_current.opts = opts
```

### New approach (reload from file) ⭐
```python
# Reload from file – keeps state clean
if original_network_path:
    if "overrides" in kwargs:
        overrides = kwargs["overrides"]
        n_current = pypsa.Network(original_network_path, override_component_attrs=overrides)
    else:
        n_current = pypsa.Network(original_network_path)
else:
    n_current = copy.deepcopy(original_network)

# Re-apply network preparation
n_current = prepare_network(
    n_current,
    kwargs.get("solve_opts", {}),
    kwargs.get("using_single_node", False),
    kwargs.get("single_node_province", "Shandong")
)
```

## Fixing Aluminum Power Consumption

### Old approach (explicit constraints)
```python
# Fix via equality constraints
def add_aluminum_usage_constraints(n, fixed_aluminum_usage):
    p = n.model["Link-p"]
    for smelter in aluminum_smelters:
        lhs = p.loc[:, smelter]
        rhs = fixed_aluminum_usage[smelter].values
        n.model.add_constraints(lhs == rhs, name=f"aluminum-fixed-usage-{smelter}")
```

### New approach (`p_set`) ⭐
```python
# Fix via p_set
for smelter in aluminum_usage.columns:
    if smelter in n_current.links.index:
        fixed_aluminum_power = aluminum_usage[smelter].values
        if not hasattr(n_current.links_t, 'p_set'):
            n_current.links_t.p_set = pd.DataFrame(index=n_current.snapshots, columns=n_current.links.index)
        n_current.links_t.p_set[smelter] = fixed_aluminum_power
        for load in aluminum_loads:
            if hasattr(n_current.loads_t, 'p_set'):
                n_current.loads_t.p_set[load] = fixed_aluminum_power
```

## Virtual Generator Marginal Cost

### Old approach (zero marginal cost)
```python
n_al.add("Generator",
         f"virtual_gen_{bus}",
         bus=bus,
         carrier="virtual",
         p_nom=1e6,
         marginal_cost=0.0)
```

### New approach (nodal price) ⭐
```python
if nodal_prices is not None and bus in nodal_prices.index:
    marginal_cost = nodal_prices[bus]
    logger.info(f"Adding virtual generator at bus {bus} with nodal marginal price")
else:
    marginal_cost = 0.0
    logger.info(f"Adding virtual generator at bus {bus} with zero marginal cost (no nodal price data)")

n_al.add("Generator",
         f"virtual_gen_{bus}",
         bus=bus,
         carrier="virtual",
         p_nom=1e6,
         marginal_cost=marginal_cost)
```

### Nodal Price Extraction
```python
current_nodal_prices = None
if hasattr(n_current, 'buses_t') and hasattr(n_current.buses_t, 'marginal_price'):
    electricity_buses = n_current.buses[n_current.buses.carrier != "aluminum"].index
    if len(electricity_buses) > 0:
        current_nodal_prices = n_current.buses_t.marginal_price[electricity_buses]
        logger.info(f"Extracted nodal prices: {list(electricity_buses)}")
```

## Advantages

1. **More stable convergence**: objective-based convergence is more robust than energy-change-based convergence.
2. **Cleaner network state**: reloading the network at each iteration avoids accumulated solver artifacts.
3. **Higher optimization accuracy**: virtual generators use nodal marginal prices, improving aluminum dispatch quality.
4. **Simpler implementation**: `p_set` fixing replaces explicit constraint injection.
5. **Better performance monitoring**: detailed timing statistics help diagnose algorithm performance.
6. **Flexible configuration**: iteration parameters are adjustable via the config file.
7. **Consistent design**: aligned with `capacity_expansion_planning_aluminum_iterative.py`.
