@ -1,879 +0,0 @@
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot cost analysis for MMMU/MMMF-2050 scenarios.

X-axis: different aluminum capacity ratios (5p–100p)
Y-axes:
- Left: cost savings (Billion CNY)
- Right: emissions reduction (Million tonnes CO2)

The figure shows:
- power-system cost savings
- aluminum operational cost changes
- net cost savings
- emissions reduction

Demand is fixed at M, flexibility at M, market opportunity at M,
and employment transfer at U/F (core scenario: MMMU, comparison: MMMF).
Positive direction: cost reduction and emissions reduction.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
import yaml
import argparse
import copy
import glob
import re

# Font settings
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

# Logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Aluminum cost weighting method (employment scenarios U and F use different weights).
# We specify weights by employment scenario (U=MMMU, F=MMMF) and cost type
# (capital/marginal/standby/other). For example:
# - capital: 0        → ignore capital cost completely
# - marginal: 1       → fully include variable cost
# - standby: 0.5      → include half of standby cost
ALUMINUM_COST_METHODS = {
    "U": {  # MMMU
        "capital": 1,
        "marginal": 1.4227,
        "standby": 0,
        "other": 1.0,
    },
    "F": {  # MMMF (can be tuned separately from MMMU if needed)
        "capital": 0.63,
        "marginal": 1.4227,
        "standby": 0.94847,
        "other": 1.0,
    },
}

# Aluminum cost breakdown ratios (for logging maintenance / labor / restart),
# specified separately for employment scenarios U (MMMU) and F (MMMF):
# - maintenance = capital change × maintenance_ratio
# - labor      = capital change × labor_capital_ratio
#                + standby change × labor_standby_ratio
# - restart_cost = startup + shutdown (if both exist, take the smaller |value| × 2)
ALUMINUM_COST_BREAKDOWN_RATIOS = {
    "U": {  # MMMU
        "maintenance_ratio": 0.625,
        "labor_capital_ratio": 0.375,
        "labor_standby_ratio": 0.0,
    },
    "F": {  # MMMF
        "maintenance_ratio": 1.0,
        "labor_capital_ratio": 0.0,
        "labor_standby_ratio": 0.1,
    },
}

def load_config(config_path):
    """
    Load a YAML config file.
    
    Parameters:
    -----------
    config_path : str or Path
        Path to config file
        
    Returns:
    --------
    dict
        Parsed config
    """
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        logger.error(f"Error loading config file {config_path}: {str(e)}")
        return None

def find_available_years(results_dir, base_version):
    """
    Find available years with postnetwork cost data.
    
    Parameters:
    -----------
    results_dir : str
        Results directory
    base_version : str
        Base version string
        
    Returns:
    --------
    list
        Sorted list of available years
    """
    available_years = []
    results_path = Path(results_dir)
    run_tags = ["FCG", "Neighbor"]  # prefer new convention, fallback to legacy
    
    # Scan for possible years
    for year in [2030, 2040, 2050]:
        # Find version directories containing this year
        year_pattern = f"version-{base_version}-*{year}*"
        version_dirs = list(results_path.glob(year_pattern))
        
        for version_dir in version_dirs:
            # Check if data exists for this year
            summary_dir = version_dir / 'summary' / 'postnetworks' / 'positive'
            if summary_dir.exists():
                for tag in run_tags:
                    year_pattern = f"postnetwork-ll-current+{tag}-linear2050-{year}"
                    year_dir = summary_dir / year_pattern
                    if year_dir.exists() and (year_dir / 'costs.csv').exists():
                        available_years.append(year)
                        break
                if year in available_years:
                    break
    
    # If nothing found, default to 2050
    if not available_years:
        available_years = [2050]
        logger.warning("No year data found, defaulting to 2050")
    
    return sorted(list(set(available_years)))

def load_costs_data(version_name, year, results_dir='results'):
    """
    Load cost data for a given version.
    
    Parameters:
    -----------
    version_name : str
        Version name, e.g. '0814.4H.2-MMMU-2050-100p'
    year : int
        Year
    results_dir : str
        Results directory
        
    Returns:
    --------
    pd.DataFrame or None
        Cost data (or None if missing)
    """
    try:
        # Build file paths (prefer FCG, then Neighbor)
        candidates = [
            Path(
                f"{results_dir}/version-{version_name}/summary/postnetworks/positive/"
                f"postnetwork-ll-current+FCG-linear2050-{year}/costs.csv"
            ),
            Path(
                f"{results_dir}/version-{version_name}/summary/postnetworks/positive/"
                f"postnetwork-ll-current+Neighbor-linear2050-{year}/costs.csv"
            ),
        ]
        file_path = next((p for p in candidates if p.exists()), candidates[0])
        
        if not file_path.exists():
            logger.warning(f"File does not exist: {file_path}")
            return None
        
        # Read CSV file
        df = pd.read_csv(file_path, header=None)
        
        # Handle multi-level index structure
        if len(df.columns) >= 4:
            # First two columns as index, third column as technology name
            df.set_index([0, 1, 2], inplace=True)
            # Rename last column to a numeric column name
            df.columns = [df.columns[0]]
            # Ensure numeric type
            df[df.columns[0]] = pd.to_numeric(df[df.columns[0]], errors='coerce')
        else:
            # Fallback: simple multi-index with two columns
            df = pd.read_csv(file_path, index_col=[0, 1])
            numeric_col = df.columns[0]
            df[numeric_col] = pd.to_numeric(df[numeric_col], errors='coerce')
        
        return df
        
    except Exception as e:
        logger.error(f"Error loading data: {str(e)}")
        return None

def calculate_cost_categories(costs_data):
    """
    Aggregate costs into high-level categories.
    
    Parameters:
    -----------
    costs_data : pd.DataFrame
        Raw cost data
        
    Returns:
    --------
    dict
        Mapping from category name to aggregated cost
    """
    if costs_data is None or costs_data.empty:
        return {}
    
    # Mapping from (cost_type, carrier) to aggregate category
    cost_category_mapping = {
        # variable cost-non-renewable
        ('marginal', 'coal'): 'variable cost-non-renewable',
        ('marginal', 'coal power plant'): 'variable cost-non-renewable',
        ('marginal', 'coal cc'): 'variable cost-non-renewable',
        ('marginal', 'gas'): 'variable cost-non-renewable',
        ('marginal', 'nuclear'): 'variable cost-non-renewable',
        ('marginal', 'CHP coal'): 'variable cost-non-renewable',
        ('marginal', 'CHP gas'): 'variable cost-non-renewable',
        ('marginal', 'OCGT gas'): 'variable cost-non-renewable',
        ('marginal', 'coal boiler'): 'variable cost-non-renewable',
        ('marginal', 'gas boiler'): 'variable cost-non-renewable',
        
        # capital-non-renewable
        ('capital', 'coal'): 'capital-non-renewable',
        ('capital', 'coal power plant'): 'capital-non-renewable',
        ('capital', 'coal cc'): 'capital-non-renewable',
        ('capital', 'gas'): 'capital-non-renewable',
        ('capital', 'nuclear'): 'capital-non-renewable',
        ('capital', 'CHP coal'): 'capital-non-renewable',
        ('capital', 'CHP gas'): 'capital-non-renewable',
        ('capital', 'OCGT gas'): 'capital-non-renewable',
        ('capital', 'coal boiler'): 'capital-non-renewable',
        ('capital', 'gas boiler'): 'capital-non-renewable',
        
        # capital-demand side
        ('capital', 'heat pump'): 'heating-electrification',
        ('capital', 'resistive heater'): 'heating-electrification',
        
        # capital-renewable
        ('capital', 'hydro_inflow'): 'capital-renewable',
        ('capital', 'hydroelectricity'): 'capital-renewable',
        ('capital', 'offwind'): 'capital-renewable',
        ('capital', 'onwind'): 'capital-renewable',
        ('capital', 'solar'): 'capital-renewable',
        ('capital', 'solar thermal'): 'capital-renewable',
        ('capital', 'biomass'): 'capital-renewable',
        ('capital', 'biogas'): 'capital-renewable',
        
        # transmission lines
        ('capital', 'AC'): 'transmission lines',
        ('capital', 'stations'): 'transmission lines',
        
        # batteries
        ('capital', 'battery'): 'batteries',
        ('capital', 'battery discharger'): 'batteries',
        ('marginal', 'battery'): 'batteries',
        ('marginal', 'battery discharger'): 'batteries',
        
        # long-duration storages
        ('capital', 'PHS'): 'long-duration storages',
        ('capital', 'water tanks'): 'long-duration storages',
        ('capital', 'H2'): 'long-duration storages',
        ('capital', 'H2 CHP'): 'long-duration storages',
        ('marginal', 'PHS'): 'long-duration storages',
        ('marginal', 'water tanks'): 'long-duration storages',
        ('marginal', 'H2'): 'long-duration storages',
        ('marginal', 'H2 CHP'): 'long-duration storages',
        
        # other categories
        ('capital', 'CO2 capture'): 'capital-non-renewable',
        ('marginal', 'CO2 capture'): 'variable cost-non-renewable',
        ('capital', 'Sabatier'): 'capital-non-renewable',
        ('marginal', 'Sabatier'): 'variable cost-non-renewable',
        ('capital', 'CO2'): 'capital-non-renewable',
        ('marginal', 'CO2'): 'variable cost-non-renewable',
        ('capital', 'DAC'): 'capital-non-renewable',
        ('marginal', 'DAC'): 'variable cost-non-renewable',
    }
    
    # Aggregate by category
    category_costs = {}
    
    for idx in costs_data.index:
        if len(idx) >= 3:
            component_type, cost_type, carrier = idx[0], idx[1], idx[2]
            
            # 使用分类映射
            category_key = (cost_type, carrier)
            category_name = cost_category_mapping.get(category_key, f"{cost_type} - {carrier}")
            
            if category_name not in category_costs:
                category_costs[category_name] = 0
            
            value = costs_data.loc[idx].iloc[0]
            if not pd.isna(value):
                category_costs[category_name] += value
    
    return category_costs

def calculate_total_emissions_from_costs(costs_data):
    """
    Approximate total emissions from cost data (via coal/gas marginal costs).
    
    Parameters:
    -----------
    costs_data : pd.DataFrame
        Raw cost data
        
    Returns:
    --------
    float
        Total emissions (tonnes CO2)
    """
    if costs_data is None or costs_data.empty:
        return 0
    
    total_emissions = 0
    
    for idx in costs_data.index:
        if len(idx) >= 3:
            component_type, cost_type, carrier = idx[0], idx[1], idx[2]
            if isinstance(carrier, str) and 'coal' in carrier.lower() and cost_type == 'marginal':
                value = costs_data.loc[idx].iloc[0]
                if pd.notna(value):
                    total_emissions += value
            elif isinstance(carrier, str) and 'gas' in carrier.lower() and cost_type == 'marginal':
                value = costs_data.loc[idx].iloc[0]
                if pd.notna(value):
                    total_emissions += value
    
    return total_emissions

def save_plot_data_to_csv(plot_data, output_dir, year, market, scenario_code):
    """
    Save all plot data to CSV files (detailed, summary, placeholder breakdown).
    
    Parameters:
    -----------
    plot_data : dict
        Dictionary with all computed results
    output_dir : str or Path
        Output directory path
    year : int
        Year
    market : str
        Market level
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 创建详细数据表格
    detailed_data = []
    for i, ratio in enumerate(plot_data['capacity_ratios']):
        row = {
            'Capacity_Ratio': ratio,
            'Capacity_Value_Mt': plot_data['capacity_values'][i],
            'Power_Cost_Changes_Billion_CNY': plot_data['power_cost_changes'][i] / 1e9,
            'Aluminum_Cost_Changes_Billion_CNY': plot_data['aluminum_cost_changes'][i] / 1e9,
            'Aluminum_Maintenance_Billion_CNY': plot_data['aluminum_maintenance_cny'][i] / 1e9,
            'Aluminum_Labor_Billion_CNY': plot_data['aluminum_labor_cny'][i] / 1e9,
            'Aluminum_Restart_Billion_CNY': plot_data['aluminum_restart_cny'][i] / 1e9,
            'Net_Cost_Savings_Billion_CNY': plot_data['net_cost_savings'][i] / 1e9,
            'Emissions_Changes_Million_Tonnes_CO2': plot_data['emissions_changes'][i],
            'Is_Max_Savings': i == plot_data['max_saving_index']
        }
        detailed_data.append(row)
    
    # 保存详细数据
    detailed_df = pd.DataFrame(detailed_data)
    prefix = scenario_code.lower()
    detailed_file = output_path / f"{prefix}_{year}_{market}_detailed_data.csv"
    detailed_df.to_csv(detailed_file, index=False, encoding='utf-8')
    logger.info(f"Detailed data saved to: {detailed_file}")
    
    # 创建汇总数据表格
    summary_data = {
        'Year': [year],
        'Market': [market],
        'Base_Version': [plot_data['base_version']],
        'Max_Net_Savings_Billion_CNY': [plot_data['max_saving_value'] / 1e9],
        'Max_Savings_Capacity_Mt': [plot_data['max_saving_capacity']],
        'Max_Savings_Capacity_Ratio': [plot_data['capacity_ratios'][plot_data['max_saving_index']]],
        'Total_Power_Cost_Savings_Billion_CNY': [sum(plot_data['power_cost_changes']) / 1e9],
        'Total_Aluminum_Cost_Changes_Billion_CNY': [sum(plot_data['aluminum_cost_changes']) / 1e9],
        'Total_Aluminum_Maintenance_Billion_CNY': [sum(plot_data['aluminum_maintenance_cny']) / 1e9],
        'Total_Aluminum_Labor_Billion_CNY': [sum(plot_data['aluminum_labor_cny']) / 1e9],
        'Total_Aluminum_Restart_Billion_CNY': [sum(plot_data['aluminum_restart_cny']) / 1e9],
        'Total_Net_Cost_Savings_Billion_CNY': [sum(plot_data['net_cost_savings']) / 1e9],
        'Total_Emissions_Reduction_Million_Tonnes_CO2': [sum(plot_data['emissions_changes'])]
    }
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = output_path / f"{prefix}_{year}_{market}_summary.csv"
    summary_df.to_csv(summary_file, index=False, encoding='utf-8')
    logger.info(f"Summary data saved to: {summary_file}")
    
    # 创建成本分类数据表格（如果可用）
    try:
        # 这里可以添加更详细的成本分类数据保存
        cost_breakdown_file = output_path / f"{prefix}_{year}_{market}_cost_breakdown.csv"
        # 暂时创建一个空的成本分解文件，后续可以扩展
        pd.DataFrame({'Note': ['Cost breakdown data will be added in future versions']}).to_csv(
            cost_breakdown_file, index=False, encoding='utf-8')
        logger.info(f"Cost breakdown placeholder saved to: {cost_breakdown_file}")
    except Exception as e:
        logger.warning(f"Could not save cost breakdown data: {str(e)}")

def plot_single_year_market(
    year,
    market,
    base_version,
    capacity_ratios,
    results_dir,
    ax,
    save_data=True,
    output_dir=None,
    flexibility="M",
    demand="M",
    employment="U",
):
    """
    Plot costs and emissions for a single (year, market) combination.
    
    Parameters:
    -----------
    year : int
        Year
    market : str
        Market level (L, M, H)
    flexibility : str
        Flexibility level (L, M, H, N). Default: M (core scenario)
    demand : str
        Demand level (L, M, H). Default: M (core scenario)
    employment : str
        Employment transfer level (U/F). Default: U (core scenario)
    base_version : str
        Base version string
    capacity_ratios : list
        Capacity ratio labels, e.g. ['5p', '10p', ...]
    results_dir : str
        Results directory
    ax : matplotlib.axes.Axes
        Matplotlib axis
    save_data : bool
        Whether to save data to CSV
    output_dir : str or Path
        Output directory path
        
    Returns:
    --------
    dict
        Dictionary with all computed results
    """
    # EUR → CNY conversion rate
    EUR_TO_CNY = 7.868
    
    # Build version names with the agreed pattern
    version_names = []
    config_versions = {}

    scenario_code = f"{flexibility}{demand}{market}{employment}"

    for ratio in capacity_ratios:
        # 版本号格式: base_version-{scenario_code}-{year}-{ratio}
        version = f"{base_version}-{scenario_code}-{year}-{ratio}"
        version_names.append(version)
        config_versions[ratio] = version
    
    # Baseline versions
    aluminum_baseline_version = f"{base_version}-{scenario_code}-{year}-5p"

    # Power-system baseline: always use U-employment non_flexible version.
    # Employment transfer mainly affects socio-economic impacts, so we keep
    # the power-system baseline fixed at U.
    power_baseline_scenario = f"{flexibility}{demand}{market}U"
    power_baseline_version = f"{base_version}-{power_baseline_scenario}-{year}-non_flexible"
    
    # Data containers
    costs_data = {}
    baseline_data = {}
    
    # Load baseline data
    aluminum_baseline = load_costs_data(aluminum_baseline_version, year, results_dir)
    if aluminum_baseline is not None:
        baseline_data['aluminum'] = aluminum_baseline
    
    power_baseline = load_costs_data(power_baseline_version, year, results_dir)
    if power_baseline is not None:
        baseline_data['power'] = power_baseline
    
    # Load data for each capacity ratio
    for ratio in capacity_ratios:
        version_name = config_versions[ratio]
        costs = load_costs_data(version_name, year, results_dir)
        if costs is not None:
            costs_data[ratio] = costs
    
    if not costs_data or not baseline_data:
        ax.text(0.5, 0.5, f'No data for {year}-{market}', ha='center', va='center',
               transform=ax.transAxes, fontsize=12)
        return
    
    # Cost changes (cost reduction as positive)
    power_cost_changes = []
    aluminum_cost_changes = []
    aluminum_maintenance_cny = []   # breakdown: maintenance (capital×ratio), CNY
    aluminum_labor_cny = []         # breakdown: labor (capital×ratio + standby×ratio), CNY
    aluminum_restart_cny = []       # breakdown: restart (startup+shutdown), CNY
    emissions_changes = []
    capacity_values = []
    
    for ratio in capacity_ratios:
        if ratio in costs_data and 'power' in baseline_data:
            # Power-system cost change (reduction as positive)
            current_costs = calculate_cost_categories(costs_data[ratio])
            baseline_costs = calculate_cost_categories(baseline_data['power'])
            
            # Aggregate all categories except aluminum-related ones
            total_change = 0
            for category, value in current_costs.items():
                if 'aluminum' not in category.lower():
                    baseline_value = baseline_costs.get(category, 0)
                    total_change += (value - baseline_value)
            
            # Cost reduction as positive → take negative value; convert to CNY
            power_cost_changes.append(-total_change * EUR_TO_CNY)
        else:
            power_cost_changes.append(0)
        
        if ratio in costs_data and 'aluminum' in baseline_data:
            # Aluminum cost change (reduction as positive), weighted by type
            current_costs = calculate_cost_categories(costs_data[ratio])
            baseline_costs = calculate_cost_categories(baseline_data['aluminum'])

            method = ALUMINUM_COST_METHODS.get(employment, ALUMINUM_COST_METHODS.get("U", {}))
            ratios = ALUMINUM_COST_BREAKDOWN_RATIOS.get(employment, ALUMINUM_COST_BREAKDOWN_RATIOS.get("U", {}))
            aluminum_change = 0
            aluminum_startup_change = 0
            aluminum_shutdown_change = 0
            has_startup = False
            has_shutdown = False
            # Weighted deltas (apply method weights first, then use for maintenance/labor/restart split)
            capital_delta_weighted = 0.0
            standby_delta_weighted = 0.0

            for category, value in current_costs.items():
                name = category.lower()
                if 'aluminum' not in name:
                    continue

                baseline_value = baseline_costs.get(category, 0)
                delta_raw = value - baseline_value

                # Cost type by name prefix
                if name.startswith('capital'):
                    weight = method.get("capital", 1.0)
                elif name.startswith('marginal'):
                    weight = method.get("marginal", 1.0)
                elif name.startswith('standby'):
                    weight = method.get("standby", 1.0)
                else:
                    weight = method.get("other", 1.0)

                if weight == 0:
                    continue

                delta = weight * delta_raw

                # Use weighted delta for cost breakdown, consistent with aluminum_change
                if name.startswith('capital'):
                    capital_delta_weighted += delta
                elif name.startswith('standby'):
                    standby_delta_weighted += delta

                # startup/shutdown: if both appear, only keep the smaller in absolute value
                if "startup" in name:
                    aluminum_startup_change += delta
                    has_startup = True
                elif "shutdown" in name:
                    aluminum_shutdown_change += delta
                    has_shutdown = True
                else:
                    aluminum_change += delta

            # Merge startup/shutdown → restart_cost (if both exist, 2×min(|startup|, |shutdown|))
            if has_startup and has_shutdown:
                chosen = (
                    aluminum_startup_change
                    if abs(aluminum_startup_change) <= abs(aluminum_shutdown_change)
                    else aluminum_shutdown_change
                )
                restart_eur = 2 * chosen
                aluminum_change += restart_eur
            elif has_startup:
                restart_eur = aluminum_startup_change
                aluminum_change += aluminum_startup_change
            elif has_shutdown:
                restart_eur = aluminum_shutdown_change
                aluminum_change += aluminum_shutdown_change
            else:
                restart_eur = 0.0

            # Breakdown in EUR:
            # maintenance = capital×ratio
            # labor      = capital×ratio + standby×ratio
            maintenance_eur = capital_delta_weighted * ratios.get("maintenance_ratio", 0.0)
            labor_eur = (
                capital_delta_weighted * ratios.get("labor_capital_ratio", 0.0)
                + standby_delta_weighted * ratios.get("labor_standby_ratio", 0.0)
            )
            # Convert to CNY; positive means cost increase
            aluminum_maintenance_cny.append(maintenance_eur * EUR_TO_CNY)
            aluminum_labor_cny.append(labor_eur * EUR_TO_CNY)
            aluminum_restart_cny.append(restart_eur * EUR_TO_CNY)

            # Cost reduction as positive → take negative; convert to CNY
            aluminum_cost_changes.append(-aluminum_change * EUR_TO_CNY)
        else:
            aluminum_cost_changes.append(0)
            aluminum_maintenance_cny.append(0.0)
            aluminum_labor_cny.append(0.0)
            aluminum_restart_cny.append(0.0)
        
        # Emissions change (reduction as positive)
        if ratio in costs_data and 'power' in baseline_data:
            current_emissions = calculate_total_emissions_from_costs(costs_data[ratio])
            baseline_emissions_total = calculate_total_emissions_from_costs(baseline_data['power'])
            
            # Emissions reduction as positive → take negative; convert to Mt CO2
            emissions_change = -(current_emissions - baseline_emissions_total)
            emissions_changes.append(emissions_change / 1e6)
        else:
            emissions_changes.append(0)
        
        # Read capacity ratio from config
        config_file = f"configs/config_{scenario_code}_{year}_{ratio}.yaml"
        config = load_config(config_file)
        if config is None:
            raise FileNotFoundError(
                f"Required config YAML not found or unreadable: {config_file}. "
                f"Please generate it first (scenario_code={scenario_code}, year={year}, ratio={ratio})."
            )

        capacity_ratio = config.get('aluminum_capacity_ratio', 1.0)
        if 'aluminum' in config and 'capacity_ratio' in config['aluminum']:
            capacity_ratio = config['aluminum']['capacity_ratio']
        
        # Actual capacity (4500 * capacity_ratio)
        actual_capacity = 4500 * capacity_ratio
        capacity_values.append(actual_capacity)
    
    # Secondary y-axis (emissions)
    ax2 = ax.twinx()
    
    # X positions and bar width
    x = capacity_values
    bar_width = 150  # 减少柱子宽度
    
    # Slightly offset positions for power vs aluminum bars
    x_power = [pos - bar_width/6 for pos in x]   # power-system costs: slight left shift
    x_aluminum = [pos + bar_width/6 for pos in x]  # aluminum costs: slight right shift
    
    # Power-system cost savings (bottom)
    bars1 = ax.bar(x_power, power_cost_changes, bar_width*0.8, color='#1f77b4', alpha=0.8, 
                   label='Power System Cost Savings')
    
    # Aluminum operational cost change, stacked on top of power-system costs
    bars2 = ax.bar(x_aluminum, aluminum_cost_changes, bar_width*0.8, bottom=power_cost_changes, 
                   color='#ff7f0e', alpha=0.8, label='Aluminum Operation Cost Increase')
    
    # Net cost savings
    net_cost_savings = [power_cost_changes[i] + aluminum_cost_changes[i] for i in range(len(capacity_values))]
    
    # Net savings curve (black line, original x positions)
    ax.plot(x, net_cost_savings, 'k-', linewidth=3, label='Net Cost Savings', marker='o', markersize=8, zorder=20)
    
    # Find maximum net savings
    max_saving_index = np.argmax(net_cost_savings)
    max_saving_value = net_cost_savings[max_saving_index]
    max_saving_capacity = capacity_values[max_saving_index]
    
    # Mark maximum with a star
    ax.plot(max_saving_capacity, max_saving_value, 'r*', markersize=15, 
            label=f'Highest Net Savings: {max_saving_value/1e9:.0f}B CNY', zorder=30)
    
    # Annotate the maximum net savings
    ax.annotate(f'{max_saving_value/1e9:.0f}B',
                xy=(max_saving_capacity, max_saving_value),
                xytext=(0, 20),
                textcoords="offset points",
                ha='center', va='bottom', 
                fontsize=12, weight='bold', color='red',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Emissions reduction (right y-axis, original x positions)
    line1 = ax2.plot(x, emissions_changes, linewidth=2, marker='o', 
                     markersize=6, label='Emissions Reduction', color='red')
    
    # Axis labels
    ax.set_xlabel('Aluminum Smelting Capacity (Mt)', fontsize=12)
    ax.set_ylabel('Cost Savings (Billion CNY)', fontsize=12, color='blue')
    ax2.set_ylabel('Emissions Reduction (Million Tonnes CO2)', fontsize=12, color='red')
    # ax.set_title(f'Year {year}, Market {market}', fontsize=14, fontweight='bold')
    
    # Zero lines
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
    ax2.axhline(y=0, color='red', linestyle='--', alpha=0.5, linewidth=1)
    
    # Grid on primary y-axis
    ax.grid(True, alpha=0.3, axis='y')
    
    # X ticks and labels
    ax.set_xticks(x)
    ax.set_xticklabels([f'{cap/100:.0f}' for cap in capacity_values], fontsize=12)
    
    # Left y-axis ticks in Billion CNY
    y_ticks = ax.get_yticks()
    y_tick_labels = [f'{tick/1e9:.0f}' for tick in y_ticks]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels, fontsize=12, color='blue')
    
    # Right y-axis ticks in Million tonnes CO2
    y2_ticks = ax2.get_yticks()
    y2_tick_labels = [f'{tick:.0f}' for tick in y2_ticks]
    ax2.set_yticks(y2_ticks)
    ax2.set_yticklabels(y2_tick_labels, fontsize=12, color='red')
    
    # Do not add legend here; a single shared legend is added at figure level
    
    # 准备返回的数据
    plot_data = {
        'capacity_values': capacity_values,
        'power_cost_changes': power_cost_changes,
        'aluminum_cost_changes': aluminum_cost_changes,
        'aluminum_maintenance_cny': aluminum_maintenance_cny,
        'aluminum_labor_cny': aluminum_labor_cny,
        'aluminum_restart_cny': aluminum_restart_cny,
        'net_cost_savings': net_cost_savings,
        'emissions_changes': emissions_changes,
        'max_saving_index': max_saving_index,
        'max_saving_value': max_saving_value,
        'max_saving_capacity': max_saving_capacity,
        'year': year,
        'market': market,
        'base_version': base_version,
        'capacity_ratios': capacity_ratios
    }
    
    # 保存数据到CSV文件
    if save_data and output_dir is not None:
        save_plot_data_to_csv(plot_data, output_dir, year, market, scenario_code)
    
    return plot_data

def plot_mmm_2050_analysis(employment="U"):
    """
    Plot cost analysis for MMMU/MMMF-2050 scenarios for a given employment case.
    """
    # Base version from main config
    main_config = load_config('config.yaml')
    if main_config is None:
        logger.error("Unable to load main config file config.yaml")
        return
    
    base_version = main_config.get('version', '0815.1H.1')
    logger.info(f"Base version read from main config file: {base_version}")
    
    # Capacity ratios
    capacity_ratios = ['5p', '10p', '15p', '20p', '30p', '40p', '50p', '60p', '70p', '80p', '90p', '100p']
    
    # Fixed MMM*-2050 setting (M demand, M flexibility, M market; employment set by argument)
    year = 2050
    market = 'M'
    flexibility = "M"
    demand = "M"
    
    # Single subplot
    fig, ax = plt.subplots(1, 1, figsize=(10, 8.5))
    
    scenario_code = f"{flexibility}{demand}{market}{employment}"
    logger.info(f"Plotting {scenario_code}-2050 scenario chart...")
    
    # Output directory: distinguish MMMU vs MMMF by scenario_tag
    scenario_tag = f"mmm{employment.lower()}"  # e.g. U→mmmu, F→mmmf
    output_dir = Path(f"results/{scenario_tag}_2050_analysis")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate plot and save data
    plot_data = plot_single_year_market(
        year,
        market,
        base_version,
        capacity_ratios,
        'results',
        ax,
        save_data=True,
        output_dir=output_dir,
        flexibility=flexibility,
        demand=demand,
        employment=employment,
    )
    
    # Legend entries
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor='#1f77b4', alpha=0.8, label='Power System Cost Savings'),
        plt.Rectangle((0,0),1,1, facecolor='#ff7f0e', alpha=0.8, label='Aluminum Operation Cost Increase'),
        plt.Line2D([0], [0], color='black', linewidth=3, marker='o', markersize=8, label='Net Cost Savings'),
        plt.Line2D([0], [0], marker='*', color='red', markersize=15, linestyle='', label='Highest Net Savings'),
        plt.Line2D([0], [0], color='red', linewidth=2, marker='o', markersize=6, label='Emissions Reduction')
    ]
    
    # Add legend
    ax.legend(handles=legend_elements, loc='lower left', fontsize=15)
    
    # Optional: overall title
    # fig.suptitle('MMM-2050 Scenario Analysis\n(Demand: M, Flexibility: M)\nCost Savings (Positive) & Emissions Reduction (Positive)', 
    #              fontsize=16, fontweight='bold', y=0.95)
    
    plt.tight_layout()
    
    # Save figure (output directory already created)
    plot_file = output_dir / f"{scenario_tag}_2050_analysis.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"{scenario_code}-2050 scenario analysis chart saved to: {plot_file}")
    
    # Log saved data filenames
    logger.info(f"Data files saved to: {output_dir}")
    logger.info(f"- Detailed data: {scenario_tag}_{year}_{market}_detailed_data.csv")
    logger.info(f"- Summary data: {scenario_tag}_{year}_{market}_summary.csv")
    logger.info(f"- Cost breakdown: {scenario_tag}_{year}_{market}_cost_breakdown.csv")
    
    # plt.show()

def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description='Plot cost analysis chart for MMMU/MMMF-2050 scenario')
    parser.add_argument('--results-dir', default='results', help='Results directory path (default: results)')
    parser.add_argument(
        '--output',
        default=None,
        help='(Deprecated; kept for backward compatibility. Output dir is derived from scenario.)',
    )
    parser.add_argument(
        '--employment',
        default='U',
        choices=['U', 'F'],
        help='Employment transfer level: U (core MMMU) or F (MMMF). Default: U.',
    )
    
    args = parser.parse_args()
    
    logger.info("Starting MMM*-2050 scenario analysis")
    logger.info(f"Results directory: {args.results_dir}")
    logger.info(f"Employment level: {args.employment}")
    
    # Plot MMM*-2050 scenario analysis chart
    plot_mmm_2050_analysis(employment=args.employment)
    
    logger.info("Analysis completed!")

if __name__ == "__main__":
    main()