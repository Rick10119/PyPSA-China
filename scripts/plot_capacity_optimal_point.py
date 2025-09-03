#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制最优点对应的容量和净价值散点图
横轴：铝产能容量 (MW)
纵轴：净价值 (十亿人民币)
不同年份用不同颜色表示
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

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def load_config(config_path):
    """
    加载配置文件
    
    Parameters:
    -----------
    config_path : str or Path
        配置文件路径
        
    Returns:
    --------
    dict
        配置内容
    """
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        logger.error(f"加载配置文件 {config_path} 时出错: {str(e)}")
        return None

def find_available_years(results_dir, base_version):
    """
    查找可用的年份数据
    
    Parameters:
    -----------
    results_dir : str
        结果目录
    base_version : str
        基础版本号
        
    Returns:
    --------
    list
        可用的年份列表
    """
    available_years = []
    results_path = Path(results_dir)
    
    # 查找所有可能的年份目录
    for year in [2030, 2040, 2050]:
        # 查找包含该年份的版本目录
        year_pattern = f"version-{base_version}-*{year}*"
        version_dirs = list(results_path.glob(year_pattern))
        
        for version_dir in version_dirs:
            # 检查是否有该年份的数据
            summary_dir = version_dir / 'summary' / 'postnetworks' / 'positive'
            if summary_dir.exists():
                year_pattern = f"postnetwork-ll-current+Neighbor-linear2050-{year}"
                year_dir = summary_dir / year_pattern
                if year_dir.exists() and (year_dir / 'costs.csv').exists():
                    available_years.append(year)
                    break
    
    # 如果没有找到任何年份，默认使用2050
    if not available_years:
        available_years = [2050]
        logger.warning("未找到任何年份的数据，默认使用2050年")
    
    return sorted(list(set(available_years)))

def load_costs_data(version_name, year, results_dir='results'):
    """
    加载指定版本的成本数据
    
    Parameters:
    -----------
    version_name : str
        版本名称，如 '0814.4H.2-MMM-2050-100p'
    year : int
        年份
    results_dir : str
        结果目录
        
    Returns:
    --------
    pd.DataFrame or None
        成本数据
    """
    try:
        # 构建文件路径
        file_path = Path(f"{results_dir}/version-{version_name}/summary/postnetworks/positive/postnetwork-ll-current+Neighbor-linear2050-{year}/costs.csv")
        
        if not file_path.exists():
            logger.warning(f"文件不存在: {file_path}")
            return None
        
        # 读取CSV文件
        df = pd.read_csv(file_path, header=None)
        
        # 处理多级索引结构
        if len(df.columns) >= 4:
            # 设置多级索引：前两列作为索引，第三列作为技术名称
            df.set_index([0, 1, 2], inplace=True)
            # 重命名最后一列为数值列名
            df.columns = [df.columns[0]]
            # 将数值列转换为数值类型
            df[df.columns[0]] = pd.to_numeric(df[df.columns[0]], errors='coerce')
        else:
            # 如果列数不足，使用默认的多级索引
            df = pd.read_csv(file_path, index_col=[0, 1])
            numeric_col = df.columns[0]
            df[numeric_col] = pd.to_numeric(df[numeric_col], errors='coerce')
        
        return df
        
    except Exception as e:
        logger.error(f"加载数据时出错: {str(e)}")
        return None

def calculate_cost_categories(costs_data):
    """
    计算成本分类
    
    Parameters:
    -----------
    costs_data : pd.DataFrame
        成本数据
        
    Returns:
    --------
    dict
        按分类组织的成本数据
    """
    if costs_data is None or costs_data.empty:
        return {}
    
    # 定义成本类型和资源组合的分类映射
    cost_category_mapping = {
        # variable cost-non-renewable - 非可再生能源可变成本
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
        
        # capital-non-renewable - 非可再生能源资本成本
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
        
        # capital-demand side - 需求侧资本成本
        ('capital', 'heat pump'): 'heating-electrification',
        ('capital', 'resistive heater'): 'heating-electrification',
        
        # capital-renewable - 可再生能源资本成本
        ('capital', 'hydro_inflow'): 'capital-renewable',
        ('capital', 'hydroelectricity'): 'capital-renewable',
        ('capital', 'offwind'): 'capital-renewable',
        ('capital', 'onwind'): 'capital-renewable',
        ('capital', 'solar'): 'capital-renewable',
        ('capital', 'solar thermal'): 'capital-renewable',
        ('capital', 'biomass'): 'capital-renewable',
        ('capital', 'biogas'): 'capital-renewable',
        
        # transmission lines - 输电线路
        ('capital', 'AC'): 'transmission lines',
        ('capital', 'stations'): 'transmission lines',
        
        # batteries - 电池储能
        ('capital', 'battery'): 'batteries',
        ('capital', 'battery discharger'): 'batteries',
        ('marginal', 'battery'): 'batteries',
        ('marginal', 'battery discharger'): 'batteries',
        
        # long-duration storages - 长时储能
        ('capital', 'PHS'): 'long-duration storages',
        ('capital', 'water tanks'): 'long-duration storages',
        ('capital', 'H2'): 'long-duration storages',
        ('capital', 'H2 CHP'): 'long-duration storages',
        ('marginal', 'PHS'): 'long-duration storages',
        ('marginal', 'water tanks'): 'long-duration storages',
        ('marginal', 'H2'): 'long-duration storages',
        ('marginal', 'H2 CHP'): 'long-duration storages',
        
        # 其他分类
        ('capital', 'CO2 capture'): 'capital-non-renewable',
        ('marginal', 'CO2 capture'): 'variable cost-non-renewable',
        ('capital', 'Sabatier'): 'capital-non-renewable',
        ('marginal', 'Sabatier'): 'variable cost-non-renewable',
        ('capital', 'CO2'): 'capital-non-renewable',
        ('marginal', 'CO2'): 'variable cost-non-renewable',
        ('capital', 'DAC'): 'capital-non-renewable',
        ('marginal', 'DAC'): 'variable cost-non-renewable',
    }
    
    # 按成本分类组织数据
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
    从成本数据中计算总碳排放（通过coal和gas的marginal成本估算）
    
    Parameters:
    -----------
    costs_data : pd.DataFrame
        成本数据
        
    Returns:
    --------
    float
        总碳排放量（吨CO2）
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



def calculate_actual_capacity_ratio(year: int, cap_ratio: float, demand_level: str) -> float:
    """
    Calculate actual capacity ratio
    
    Args:
        year: Year
        cap_ratio: Excess capacity retention ratio (e.g., 0.1 means 10%)
        demand_level: Demand level ('mid')
        
    Returns:
        Actual capacity ratio
    """
    # Total capacity
    total_capacity = 4500
    
    # Set demand by year (using actual data)
    demand_by_year = {
        "2030": 2902.417177819193,
        "2040": 1508.1703393209764,
        "2050": 1166.6836345743664,
    }
    
    demand = demand_by_year.get(str(year), 0)
    
    # Calculate actual capacity ratio: demand/capacity × (1-cap) + cap
    actual_ratio = (demand / total_capacity) * (1 - cap_ratio) + cap_ratio
    
    return actual_ratio

def find_optimal_points(base_version, capacity_ratios, results_dir='results'):
    """
    Find optimal points for each year-market-flexibility combination
    
    Parameters:
    -----------
    base_version : str
        Base version number
    capacity_ratios : list
        List of capacity ratios
    results_dir : str
        Results directory
        
    Returns:
    --------
    list
        List containing optimal point information, each element is (year, market, flexibility, capacity, net_value)
    """
    # Euro to CNY conversion rate
    EUR_TO_CNY = 7.8
    
    # Find available years
    available_years = find_available_years(results_dir, base_version)
    logger.info(f"Found available years: {available_years}")
    
    # Define market opportunities and flexibility levels
    markets = ['L', 'M', 'H']
    flexibilities = ['L', 'M', 'H', 'N']  # low, mid, high, non_constrained
    
    optimal_points = []
    
    for year in available_years:
        for market in markets:
            for flexibility in flexibilities:
                logger.info(f"Analyzing optimal point for {year}-{market}-{flexibility}...")
                
                # Build version names - use different base_version for different scenarios
                # Extract the base part and construct scenario-specific base_version
                base_parts = base_version.split('-')
                if len(base_parts) >= 2:
                    # Use the scenario part (e.g., MMM, HML, etc.) from base_version
                    scenario_part = base_parts[1]  # e.g., "MMM", "HML", etc.
                    scenario_base_version = f"{base_parts[0]}-{scenario_part}"
                else:
                    # Fallback to original base_version
                    scenario_base_version = base_version
                
                version_names = []
                config_versions = {}
                
                for ratio in capacity_ratios:
                    # Version format: scenario_base_version-{flexibility}{demand}{market}-{year}-{ratio}
                    # Demand is fixed as 'M' (mid)
                    version = f"{scenario_base_version}-{flexibility}M{market}-{year}-{ratio}"
                    version_names.append(version)
                    config_versions[ratio] = version
                
                # Baseline versions
                aluminum_baseline_version = f"{scenario_base_version}-{flexibility}M{market}-{year}-5p"
                power_baseline_version = f"{scenario_base_version}-{flexibility}M{market}-{year}-non_flexible"
                
                # Collect data
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
                    logger.warning(f"No data found for {year}-{market}-{flexibility}")
                    continue
                
                # Calculate net values (same method as plot_capacity_multi_year_market_comparison)
                power_cost_changes = []
                aluminum_cost_changes = []
                capacity_values = []
                
                for ratio in capacity_ratios:
                    if ratio in costs_data and 'power' in baseline_data:
                        # Calculate power system cost changes (cost reduction is positive)
                        current_costs = calculate_cost_categories(costs_data[ratio])
                        baseline_costs = calculate_cost_categories(baseline_data['power'])
                        
                        # Calculate total cost change (excluding aluminum related)
                        power_cost_change = 0
                        for category, value in current_costs.items():
                            if 'aluminum' not in category.lower():
                                baseline_value = baseline_costs.get(category, 0)
                                power_cost_change += (value - baseline_value)
                        
                        # Cost reduction is positive direction, so take negative value
                        power_cost_changes.append(-power_cost_change * EUR_TO_CNY)
                    else:
                        power_cost_changes.append(0)
                    
                    if ratio in costs_data and 'aluminum' in baseline_data:
                        # Calculate aluminum cost changes (cost reduction is positive)
                        current_costs = calculate_cost_categories(costs_data[ratio])
                        baseline_costs = calculate_cost_categories(baseline_data['aluminum'])
                        
                        # Calculate aluminum related cost changes
                        aluminum_cost_change = 0
                        for category, value in current_costs.items():
                            if 'aluminum' in category.lower():
                                baseline_value = baseline_costs.get(category, 0)
                                aluminum_cost_change += (value - baseline_value)
                        
                        # Cost reduction is positive direction, so take negative value
                        aluminum_cost_changes.append(-aluminum_cost_change * EUR_TO_CNY)
                    else:
                        aluminum_cost_changes.append(0)
                    
                    # Read capacity ratio values and calculate actual capacity
                    scenario_suffix = f"{flexibility}M{market}"
                    config_file = f"configs/config_{scenario_suffix}_{year}_{ratio}.yaml"
                    config = load_config(config_file)
                    capacity_ratio = config.get('aluminum_capacity_ratio', 1.0)
                    if 'aluminum' in config and 'capacity_ratio' in config['aluminum']:
                        capacity_ratio = config['aluminum']['capacity_ratio']
                    
                    # Calculate actual capacity using the same method as generate_capacity_test_configs
                    # Convert cap_ratio from percentage to decimal (e.g., "10p" -> 0.1)
                    cap_ratio_decimal = float(ratio.replace('p', '')) / 100.0
                    actual_capacity_ratio = calculate_actual_capacity_ratio(year, cap_ratio_decimal, 'mid')
                    # Calculate actual capacity in 10,000 tons/year (4500 * actual_capacity_ratio)
                    actual_capacity = 4500 * actual_capacity_ratio
                    capacity_values.append(actual_capacity)
                
                # Calculate net cost savings (same as plot_capacity_multi_year_market_comparison)
                net_cost_savings = [power_cost_changes[i] + aluminum_cost_changes[i] for i in range(len(capacity_values))]
                
                # Find the point with maximum net savings
                if net_cost_savings:
                    max_saving_index = np.argmax(net_cost_savings)
                    optimal_capacity = capacity_values[max_saving_index]
                    optimal_net_value = net_cost_savings[max_saving_index]
                    
                    optimal_points.append({
                        'year': year,
                        'market': market,
                        'flexibility': flexibility,
                        'capacity': optimal_capacity,
                        'net_value': optimal_net_value / 1e9  # Convert to billion CNY
                    })
                    
                    logger.info(f"{year}-{market}-{flexibility} optimal point: capacity={optimal_capacity:.1f} 10,000 tons/year, net_value={optimal_net_value/1e9:.2f}B CNY")
    
    return optimal_points

def plot_optimal_points_scatter():
    """
    Plot scatter chart of optimal points showing capacity and net value
    """
    # Load base version from main config file
    main_config = load_config('config.yaml')
    if main_config is None:
        logger.error("Cannot load main config file config.yaml")
        return
    
    base_version = main_config.get('version', '0815.1H.1')
    logger.info(f"Loaded base version from main config: {base_version}")
    
    # Define capacity ratios
    capacity_ratios = ['5p', '10p', '20p', '30p', '40p', '50p', '60p', '70p', '80p', '90p', '100p']
    
    # Find all optimal points
    optimal_points = find_optimal_points(base_version, capacity_ratios, 'results')
    
    if not optimal_points:
        logger.error("No optimal point data found")
        return
    
    # Create scatter plot
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Group data by year
    years = sorted(list(set([point['year'] for point in optimal_points])))
    colors = plt.cm.viridis(np.linspace(0, 1, len(years)))
    year_colors = dict(zip(years, colors))
    
    # Create different markers for market opportunities
    markets = ['L', 'M', 'H']
    market_markers = {'L': 'o', 'M': 's', 'H': '^'}
    
    # Create different colors for flexibility levels
    flexibilities = ['L', 'M', 'H', 'N']
    flexibility_colors = {'L': 'blue', 'M': 'green', 'H': 'orange', 'N': 'red'}
    
    # Plot scatter points
    for point in optimal_points:
        year = point['year']
        market = point['market']
        flexibility = point['flexibility']
        capacity = point['capacity']
        net_value = point['net_value']
        
        # Use year color for main color, flexibility for edge color
        ax.scatter(capacity, net_value, 
                  c=year_colors[year], 
                  marker=market_markers[market],
                  s=200, alpha=0.7, 
                  edgecolors=flexibility_colors[flexibility], 
                  linewidth=2)
    
    # Add labels
    ax.set_xlabel('Aluminum Capacity (10,000 tons/year)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Net Value (Billion CNY)', fontsize=14, fontweight='bold')
    ax.set_title('Optimal Points: Capacity vs Net Value Analysis', fontsize=16, fontweight='bold')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add vertical line for annual demand (1166.6836345743664 万吨)
    annual_demand = 1166.6836345743664  # 万吨
    ax.axvline(x=annual_demand, color='red', linestyle='--', linewidth=2, alpha=0.8, 
               label=f'Annual Demand: {annual_demand:.1f} 10,000 tons/year')
    
    # Create legend
    legend_elements = []
    
    # Year legend
    for year in years:
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor=year_colors[year], 
                                        markersize=12, label=f'Year {year}'))
    
    # Market legend
    for market in markets:
        market_desc = {'L': 'Low Market', 'M': 'Mid Market', 'H': 'High Market'}
        legend_elements.append(plt.Line2D([0], [0], marker=market_markers[market], 
                                        color='w', markerfacecolor='gray', 
                                        markersize=12, label=market_desc[market]))
    
    # Flexibility legend
    for flexibility in flexibilities:
        flex_desc = {'L': 'Low Flexibility', 'M': 'Mid Flexibility', 'H': 'High Flexibility', 'N': 'Non-constrained'}
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor='w', markeredgecolor=flexibility_colors[flexibility],
                                        markersize=12, markeredgewidth=2, label=flex_desc[flexibility]))
    
    # Add annual demand line to legend
    legend_elements.append(plt.Line2D([0], [0], color='red', linestyle='--', linewidth=2, 
                                    label=f'Annual Demand: {annual_demand:.1f} 10,000 tons/year'))
    
    # Add legend
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10, ncol=2)
    
    # Add labels for each point
    for point in optimal_points:
        year = point['year']
        market = point['market']
        flexibility = point['flexibility']
        capacity = point['capacity']
        net_value = point['net_value']
        
        ax.annotate(f'{year}-{market}-{flexibility}', 
                   xy=(capacity, net_value),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=9, alpha=0.8)
    
    plt.tight_layout()
    
    # Save plot
    output_dir = Path('results/optimal_points_analysis')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    plot_file = output_dir / "optimal_points_scatter.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Optimal points scatter plot saved to: {plot_file}")
    
    # Show plot
    plt.show()

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Plot scatter chart of optimal points showing capacity and net value')
    parser.add_argument('--results-dir', default='results', help='Results directory path (default: results)')
    parser.add_argument('--output', default='results/optimal_points_analysis', help='Output directory')
    
    args = parser.parse_args()
    
    logger.info(f"Starting analysis of optimal points for capacity and net value")
    logger.info(f"Results directory: {args.results_dir}")
    logger.info(f"Output directory: {args.output}")
    
    # Plot optimal points scatter chart
    plot_optimal_points_scatter()
    
    logger.info("Analysis completed!")

if __name__ == "__main__":
    main()
