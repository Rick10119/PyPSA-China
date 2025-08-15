#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制2050年份、M市场机会下的成本分析图表
横轴：不同产能容量（wton/year）
纵轴：电力系统节约的成本（上方）和电解铝的运行成本（下方）
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
import yaml

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

def load_costs_data(version_name, year=2050):
    """
    加载指定版本的成本数据
    
    Parameters:
    -----------
    version_name : str
        版本名称，如 '0814.4H.1-MMM-100p'
    year : int
        年份
        
    Returns:
    --------
    pd.DataFrame or None
        成本数据
    """
    try:
        # 构建文件路径
        file_path = Path(f"results/version-{version_name}/summary/postnetworks/positive/postnetwork-ll-current+Neighbor-linear2050-{year}/costs.csv")
        
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
        分类后的成本数据
    """
    if costs_data is None or costs_data.empty:
        return {}
    
    # 成本分类映射
    cost_categories = {
        'Power System Cost': 0.0,
        'Aluminum Operation Cost': 0.0
    }
    
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    for idx in costs_data.index:
        if len(idx) >= 3:
            component_type, cost_type, carrier = idx[0], idx[1], idx[2]
            value = costs_data.loc[idx].iloc[0]
            
            if pd.isna(value):
                value = 0.0
            
            # 转换为人民币
            value_cny = value * EUR_TO_CNY
            
            # 分类
            if isinstance(carrier, str) and 'aluminum' in carrier.lower():
                cost_categories['Aluminum Operation Cost'] += value_cny
            else:
                cost_categories['Power System Cost'] += value_cny
    
    return cost_categories

def get_scenario_info_from_config(config_file):
    """
    从配置文件中提取场景信息
    
    Parameters:
    -----------
    config_file : str
        配置文件路径
        
    Returns:
    --------
    dict
        场景信息
    """
    config = load_config(config_file)
    if config is None:
        return {}
    
    scenario_info = {}
    
    # 尝试从版本号中提取场景信息
    version = config.get('version', '')
    if 'MMM' in version:
        scenario_info['demand'] = 'M'  # 中等需求
        scenario_info['market'] = 'M'  # 中等市场机会
        scenario_info['flexibility'] = 'M'  # 中等灵活性
    elif 'MML' in version:
        scenario_info['demand'] = 'M'
        scenario_info['market'] = 'M'
        scenario_info['flexibility'] = 'L'  # 低灵活性
    elif 'MMH' in version:
        scenario_info['demand'] = 'M'
        scenario_info['market'] = 'M'
        scenario_info['flexibility'] = 'H'  # 高灵活性
    elif 'LMM' in version:
        scenario_info['demand'] = 'L'  # 低需求
        scenario_info['market'] = 'M'
        scenario_info['flexibility'] = 'M'
    elif 'HMM' in version:
        scenario_info['demand'] = 'H'  # 高需求
        scenario_info['market'] = 'M'
        scenario_info['flexibility'] = 'M'
    
    # 从配置文件中读取其他场景信息
    if 'add_aluminum' in config:
        scenario_info['add_aluminum'] = config['add_aluminum']
    
    if 'aluminum_capacity_ratio' in config:
        scenario_info['aluminum_capacity_ratio'] = config['aluminum_capacity_ratio']
    
    return scenario_info

def plot_2050_m_market_costs():
    """
    Plot cost analysis chart for 2050 year, M market opportunity scenario
    """
    # Read base version number from main config file
    main_config = load_config('config.yaml')
    if main_config is None:
        logger.error("Cannot load main config file config.yaml")
        return
    
    base_version = main_config.get('version', '0814.4H.1')
    logger.info(f"Read base version number from main config: {base_version}")
    
    # Define capacity ratios
    capacity_ratios = ['5p', '10p', '20p', '30p', '40p', '50p', '60p', '70p', '80p', '90p', '100p']
    
    # Build version names (M market opportunity, medium flexibility)
    version_names = []
    config_versions = {}
    
    for ratio in capacity_ratios:
        # Version number format: base_version-MMM-ratio
        version = f"{base_version}-MMM-{ratio}"
        version_names.append(version)
        config_versions[ratio] = version
        logger.info(f"Built version number: {version}")
    
    # Baseline versions: aluminum cost baseline is 5p, power system cost baseline is non-flexible
    aluminum_baseline_version = f"{base_version}-MMM-5p"
    power_baseline_version = f"{base_version}-MMM-non-flexible"
    
    logger.info(f"Aluminum cost baseline version: {aluminum_baseline_version}")
    logger.info(f"Power system cost baseline version: {power_baseline_version}")
    
    # Scenario information (fixed as M market opportunity, medium flexibility, 2050)
    scenario_info = {
        'demand': 'M',      # Medium demand
        'market': 'M',      # Medium market opportunity
        'flexibility': 'M', # Medium flexibility
        'year': 2050
    }
    logger.info(f"Scenario information: {scenario_info}")
    
    # Collect data
    costs_data = {}
    baseline_data = {}
    
    # Load baseline version data
    logger.info(f"Loading aluminum cost baseline version: {aluminum_baseline_version}")
    aluminum_baseline = load_costs_data(aluminum_baseline_version)
    if aluminum_baseline is not None:
        baseline_data['aluminum'] = aluminum_baseline
        logger.info("Successfully loaded aluminum cost baseline data")
    else:
        logger.warning(f"Cannot load aluminum cost baseline version {aluminum_baseline_version}")
    
    logger.info(f"Loading power system cost baseline version: {power_baseline_version}")
    power_baseline = load_costs_data(power_baseline_version)
    if power_baseline is not None:
        baseline_data['power'] = power_baseline
        logger.info("Successfully loaded power system cost baseline data")
    else:
        logger.warning(f"Cannot load power system cost baseline version {power_baseline_version}")
    
    # Load data for each capacity ratio
    for ratio, version_name in zip(capacity_ratios, version_names):
        logger.info(f"Loading version: {version_name}")
        data = load_costs_data(version_name)
        if data is not None:
            costs_data[ratio] = data
        else:
            logger.warning(f"Cannot load data for version {version_name}, will use empty data")
            # Create empty DataFrame to ensure data structure consistency
            costs_data[ratio] = pd.DataFrame()
    
    # Load non-flexible baseline data for emissions comparison
    non_flexible_version = f"{base_version}-MMM-non-flexible"
    logger.info(f"Loading non-flexible version for emissions baseline: {non_flexible_version}")
    non_flexible_data = load_costs_data(non_flexible_version)
    if non_flexible_data is not None:
        costs_data['non-flexible'] = non_flexible_data
        logger.info("Successfully loaded non-flexible baseline data for emissions")
    else:
        logger.warning(f"Cannot load non-flexible version {non_flexible_version}, emissions comparison will not be available")
        costs_data['non-flexible'] = pd.DataFrame()
    
    if not costs_data:
        logger.error("No valid data found")
        return
    
    # Calculate categorized costs and changes relative to baseline for each version
    categorized_costs = {}
    cost_changes = {}
    
    for ratio in capacity_ratios:
        if ratio in costs_data:
            current_costs = calculate_cost_categories(costs_data[ratio])
            categorized_costs[ratio] = current_costs
            
            # Calculate changes relative to baseline
            cost_changes[ratio] = {}
            
            # Aluminum cost changes (relative to 5p)
            # 检查基准数据和当前数据是否都存在且有效
            if ('aluminum' in baseline_data and not baseline_data['aluminum'].empty and 
                'Aluminum Operation Cost' in current_costs and current_costs['Aluminum Operation Cost'] is not None):
                aluminum_baseline_cost = calculate_cost_categories(baseline_data['aluminum']).get('Aluminum Operation Cost', 0)
                current_aluminum_cost = current_costs.get('Aluminum Operation Cost', 0)
                # 确保两个值都不是NaN或None
                if pd.notna(aluminum_baseline_cost) and pd.notna(current_aluminum_cost):
                    cost_changes[ratio]['Aluminum Operation Cost'] = aluminum_baseline_cost - current_aluminum_cost
                    logger.info(f"{ratio}: 电解铝成本变化计算成功 - 基准: {aluminum_baseline_cost:.2f}, 当前: {current_aluminum_cost:.2f}, 变化: {cost_changes[ratio]['Aluminum Operation Cost']:.2f}")
                else:
                    cost_changes[ratio]['Aluminum Operation Cost'] = 0
                    logger.warning(f"{ratio}: 电解铝成本数据缺失，变化值设为0 - 基准: {aluminum_baseline_cost}, 当前: {current_aluminum_cost}")
            else:
                cost_changes[ratio]['Aluminum Operation Cost'] = 0
                if 'aluminum' not in baseline_data or baseline_data['aluminum'].empty:
                    logger.warning(f"{ratio}: 电解铝基准数据缺失，变化值设为0")
                if 'Aluminum Operation Cost' not in current_costs or current_costs['Aluminum Operation Cost'] is None:
                    logger.warning(f"{ratio}: 电解铝当前成本数据缺失，变化值设为0")
            
            # Power system cost changes (relative to non-flexible)
            # 检查基准数据和当前数据是否都存在且有效
            if ('power' in baseline_data and not baseline_data['power'].empty and 
                'Power System Cost' in current_costs and current_costs['Power System Cost'] is not None):
                power_baseline_cost = calculate_cost_categories(baseline_data['power']).get('Power System Cost', 0)
                current_power_cost = current_costs.get('Power System Cost', 0)
                # 确保两个值都不是NaN或None
                if pd.notna(power_baseline_cost) and pd.notna(current_power_cost):
                    cost_changes[ratio]['Power System Cost'] = power_baseline_cost - current_power_cost
                    logger.info(f"{ratio}: 电力系统成本变化计算成功 - 基准: {power_baseline_cost:.2f}, 当前: {current_power_cost:.2f}, 变化: {cost_changes[ratio]['Power System Cost']:.2f}")
                else:
                    cost_changes[ratio]['Power System Cost'] = 0
                    logger.warning(f"{ratio}: 电力系统成本数据缺失，变化值设为0 - 基准: {power_baseline_cost}, 当前: {current_power_cost}")
            else:
                cost_changes[ratio]['Power System Cost'] = 0
                if 'power' not in baseline_data or baseline_data['power'].empty:
                    logger.warning(f"{ratio}: 电力系统基准数据缺失，变化值设为0")
                if 'Power System Cost' not in current_costs or current_costs['Power System Cost'] is None:
                    logger.warning(f"{ratio}: 电力系统当前成本数据缺失，变化值设为0")
        else:
            # If data is missing, create empty data
            categorized_costs[ratio] = {'Power System Cost': 0, 'Aluminum Operation Cost': 0}
            cost_changes[ratio] = {'Power System Cost': 0, 'Aluminum Operation Cost': 0}
    
    if not categorized_costs:
        logger.error("Cannot calculate cost categories")
        return
    
    # Prepare plotting data (using cost changes, positive values indicate cost reduction)
    ratios = list(categorized_costs.keys())
    power_cost_changes = [cost_changes[ratio].get('Power System Cost', 0) / 1e9 for ratio in ratios]  # Convert to billion CNY
    aluminum_cost_changes = [cost_changes[ratio].get('Aluminum Operation Cost', 0) / 1e9 for ratio in ratios]  # Convert to billion CNY
    
    # Read capacity ratio values from config files
    capacity_values = []
    for ratio in ratios:
        config_file = f"configs/config_{ratio}.yaml"
        config = load_config(config_file)
        if config is not None:
            # Read aluminum_capacity_ratio from config file
            capacity_ratio = config.get('aluminum_capacity_ratio', 1.0)
            # If not in config file, try to read from aluminum section
            if 'aluminum' in config and 'capacity_ratio' in config['aluminum']:
                capacity_ratio = config['aluminum']['capacity_ratio']
            
            # Calculate actual capacity (4500 * capacity ratio)
            actual_capacity = 4500 * capacity_ratio
            capacity_values.append(actual_capacity)
            logger.info(f"Read capacity ratio from {config_file}: {capacity_ratio}, actual capacity: {actual_capacity:.0f}MW")
        else:
            # If config file doesn't exist, use default values
            default_ratio = float(ratio.replace('p', '')) / 100.0
            default_capacity = 4500 * default_ratio
            capacity_values.append(default_capacity)
            logger.warning(f"Cannot load config file {config_file}, using default capacity ratio: {default_ratio}, capacity: {default_capacity:.0f}MW")
    
    # Create single chart with dual y-axes
    fig, ax = plt.subplots(figsize=(14, 8))
    ax2 = ax.twinx()  # Create secondary y-axis
    
    # Don't set ylim here - let each axis find its own range
    # We'll align the 0 points later using zero lines
    
    # Use actual capacity values for x-axis
    x = capacity_values
    bar_width = 225  # Adjust bar width based on capacity scale
    
    # Plot power system cost savings (stacked on top)
    bars1 = ax.bar(x, power_cost_changes, bar_width, color='#1f77b4', alpha=0.8, label='Power System Cost Savings')
    
    # Plot aluminum operation cost savings (stacked below)
    bars2 = ax.bar(x, aluminum_cost_changes, bar_width, bottom=power_cost_changes, color='#ff7f0e', alpha=0.8, label='Aluminum Operation Cost Savings')
    
    # Calculate net cost savings (total savings)
    net_cost_savings = [power_cost_changes[i] + aluminum_cost_changes[i] for i in range(len(ratios))]
    
    # Plot net cost curve (black line)
    ax.plot(x, net_cost_savings, 'k-', linewidth=3, label='Net Cost Savings', marker='o', markersize=8, zorder=20)
    
    # Calculate emissions changes relative to non-flexible baseline
    emissions_changes = []
    for ratio in ratios:
        if ratio in costs_data:
            # Get emissions data for current ratio
            current_emissions = costs_data[ratio]
            # Get emissions data for non-flexible baseline
            baseline_emissions = costs_data.get('non-flexible', pd.DataFrame())
            
            # 检查基准数据和当前数据是否都存在且有效
            if not current_emissions.empty and not baseline_emissions.empty:
                # Calculate total emissions (coal + gas)
                current_total = 0
                baseline_total = 0
                
                # Look for emissions data in the costs file
                for idx in current_emissions.index:
                    if len(idx) >= 3:
                        component_type, cost_type, carrier = idx[0], idx[1], idx[2]
                        if isinstance(carrier, str) and 'coal' in carrier.lower() and cost_type == 'marginal':
                            current_value = current_emissions.loc[idx].iloc[0]
                            if pd.notna(current_value):
                                current_total += current_value
                        elif isinstance(carrier, str) and 'gas' in carrier.lower() and cost_type == 'marginal':
                            current_value = current_emissions.loc[idx].iloc[0]
                            if pd.notna(current_value):
                                current_total += current_value
                
                for idx in baseline_emissions.index:
                    if len(idx) >= 3:
                        component_type, cost_type, carrier = idx[0], idx[1], idx[2]
                        if isinstance(carrier, str) and 'coal' in carrier.lower() and cost_type == 'marginal':
                            baseline_value = baseline_emissions.loc[idx].iloc[0]
                            if pd.notna(baseline_value):
                                baseline_total += baseline_value
                        elif isinstance(carrier, str) and 'gas' in carrier.lower() and cost_type == 'marginal':
                            baseline_value = baseline_emissions.loc[idx].iloc[0]
                            if pd.notna(baseline_value):
                                baseline_total += baseline_value
                
                # 只有当两个值都有效时才计算变化量
                if pd.notna(current_total) and pd.notna(baseline_total):
                    # Calculate emissions change (baseline - current, positive means reduction)
                    emissions_change = baseline_total - current_total
                    emissions_changes.append(emissions_change / 1e6)  # Convert to million tonnes CO2
                    logger.info(f"{ratio}: 排放变化计算成功 - 基准: {baseline_total:.2f}, 当前: {current_total:.2f}, 变化: {emissions_change/1e6:.2f}M tonnes CO2")
                else:
                    emissions_changes.append(0)
                    logger.warning(f"{ratio}: 排放数据缺失，变化值设为0 - 基准: {baseline_total}, 当前: {current_total}")
            else:
                emissions_changes.append(0)
                if current_emissions.empty:
                    logger.warning(f"{ratio}: 当前排放数据缺失，变化值设为0")
                if baseline_emissions.empty:
                    logger.warning(f"{ratio}: 排放基准数据缺失，变化值设为0")
        else:
            emissions_changes.append(0)
            logger.warning(f"{ratio}: 排放数据文件缺失，变化值设为0")
    
    # Plot emissions changes on right y-axis (green line)
    ax2.plot(x, emissions_changes, 'g-', linewidth=3, label='Emissions Change', marker='s', markersize=8, zorder=25)
    
    # Find the point with highest total savings
    max_saving_index = np.argmax(net_cost_savings)
    max_saving_value = net_cost_savings[max_saving_index]
    max_saving_ratio = ratios[max_saving_index]
    max_saving_capacity = capacity_values[max_saving_index]
    
    # Mark the highest point with star
    ax.plot(max_saving_capacity, max_saving_value, 'r*', markersize=15, label=f'Highest Savings Point: {max_saving_ratio}', zorder=30)
    
    # Add value label for the highest point
    ax.annotate(f'{max_saving_value:.1f}B',
                xy=(max_saving_capacity, max_saving_value),
                xytext=(0, 20),
                textcoords="offset points",
                ha='center', va='bottom', 
                fontsize=12, weight='bold', color='red',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Set chart properties
    ax.set_xlabel('Capacity (wton/year)', fontsize=12)
    ax.set_ylabel('Cost Savings (Billion CNY)', fontsize=12, color='black')
    ax2.set_ylabel('Emissions Change (Million Tonnes CO2)', fontsize=12, color='green')
    
    # Generate title dynamically with scenario information
    title = f'2050 Year M Market Opportunity - Cost Savings and Emissions Analysis for Different Capacities\n'
    title += f'(Baseline: Aluminum cost relative to 5p, Power system cost relative to non-flexible, Emissions relative to non-flexible)'
    title += f'\nScenario: Demand={scenario_info["demand"]}, Market={scenario_info["market"]}, Flexibility={scenario_info["flexibility"]}, Year={scenario_info["year"]}'
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Combine legends from both axes
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
    
    # Add zero lines for both axes to ensure alignment
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.3, linewidth=1)
    ax2.axhline(y=0, color='green', linestyle='-', alpha=0.3, linewidth=1)
    
    # Only align the 0 points, keep different scales for better readability
    # The zero lines will visually show the alignment
    
    # Set x-axis ticks and labels
    ax.set_xticks(capacity_values)
    # Capacity values are already in wton/year
    x_labels = [f"{capacity:.0f}" for capacity in capacity_values]
    ax.set_xticklabels(x_labels, rotation=45, ha='right')
    
    # Add value labels
    for i, (capacity, power_change, aluminum_change) in enumerate(zip(capacity_values, power_cost_changes, aluminum_cost_changes)):
        # Power system cost labels (at the top of power system bars)
        if abs(power_change) > 0.01:
            ax.text(capacity, power_change + 0.1, f'{power_change:.1f}B', 
                   ha='center', va='bottom', fontweight='bold', fontsize=9, color='#1f77b4')
        
        # Aluminum cost labels (in the center of aluminum bars)
        if abs(aluminum_change) > 0.01:
            ax.text(capacity, power_change + aluminum_change/2, f'{aluminum_change:.1f}B', 
                   ha='center', va='center', fontweight='bold', fontsize=9, color='white')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save chart
    output_dir = Path("results/plots")
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_file = output_dir / "2050_m_market_costs_analysis.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Chart saved to: {plot_file}")
    
    # Display chart
    plt.show()
    
    # Save data
    data_rows = []
    for ratio, capacity in zip(ratios, capacity_values):
        power_cost_change = cost_changes[ratio].get('Power System Cost', 0) / 1e9
        aluminum_cost_change = cost_changes[ratio].get('Aluminum Operation Cost', 0) / 1e9
        power_cost_absolute = categorized_costs[ratio].get('Power System Cost', 0) / 1e9
        aluminum_cost_absolute = categorized_costs[ratio].get('Aluminum Operation Cost', 0) / 1e9
        
        # Get emissions change for this ratio
        emissions_change = 0
        if ratio in costs_data and 'non-flexible' in costs_data:
            current_emissions = costs_data[ratio]
            baseline_emissions = costs_data['non-flexible']
            
            if not current_emissions.empty and not baseline_emissions.empty:
                current_total = 0
                baseline_total = 0
                
                # Calculate emissions from coal and gas
                for idx in current_emissions.index:
                    if len(idx) >= 3:
                        component_type, cost_type, carrier = idx[0], idx[1], idx[2]
                        if isinstance(carrier, str) and 'coal' in carrier.lower() and cost_type == 'marginal':
                            current_value = current_emissions.loc[idx].iloc[0]
                            if pd.notna(current_value):
                                current_total += current_value
                        elif isinstance(carrier, str) and 'gas' in carrier.lower() and cost_type == 'marginal':
                            current_value = current_emissions.loc[idx].iloc[0]
                            if pd.notna(current_value):
                                current_total += current_value
                
                for idx in baseline_emissions.index:
                    if len(idx) >= 3:
                        component_type, cost_type, carrier = idx[0], idx[1], idx[2]
                        if isinstance(carrier, str) and 'coal' in carrier.lower() and cost_type == 'marginal':
                            baseline_value = baseline_emissions.loc[idx].iloc[0]
                            if pd.notna(baseline_value):
                                baseline_total += baseline_value
                        elif isinstance(carrier, str) and 'gas' in carrier.lower() and cost_type == 'marginal':
                            baseline_value = baseline_emissions.loc[idx].iloc[0]
                            if pd.notna(baseline_value):
                                baseline_total += baseline_value
                
                # 只有当两个值都有效时才计算变化量
                if pd.notna(current_total) and pd.notna(baseline_total):
                    emissions_change = (baseline_total - current_total) / 1e6  # Convert to million tonnes CO2
                    logger.info(f"{ratio}: 数据导出 - 排放变化计算成功: {emissions_change:.2f}M tonnes CO2")
                else:
                    emissions_change = 0
                    logger.warning(f"{ratio}: 数据导出 - 排放数据缺失，变化值设为0")
        
        data_rows.append({
            'Capacity Ratio': ratio,
            'Capacity (wton/year)': capacity,
            'Power System Cost Change (Billion CNY)': power_cost_change,
            'Aluminum Operation Cost Change (Billion CNY)': aluminum_cost_change,
            'Power System Cost Absolute (Billion CNY)': power_cost_absolute,
            'Aluminum Operation Cost Absolute (Billion CNY)': aluminum_cost_absolute,
            'Total Cost Change (Billion CNY)': power_cost_change + aluminum_cost_change,
            'Emissions Change (Million Tonnes CO2)': emissions_change
        })
    
    data_df = pd.DataFrame(data_rows)
    data_file = output_dir / "2050_m_market_costs_analysis.csv"
    data_df.to_csv(data_file, index=False)
    logger.info(f"Data saved to: {data_file}")
    
    # Print summary information
    print("\n=== 2050 Year M Market Opportunity Cost Savings and Emissions Analysis Summary ===")
    print(f"{'Capacity Ratio':<15} {'Capacity(wton/year)':<15} {'Power System Savings(B CNY)':<25} {'Aluminum Savings(B CNY)':<25} {'Total Savings(B CNY)':<20} {'Emissions Change(M Tonnes CO2)':<25}")
    print("-" * 130)
    for row in data_rows:
        print(f"{row['Capacity Ratio']:<15} {row['Capacity (wton/year)']:<15.0f} {row['Power System Cost Change (Billion CNY)']:<25.2f} {row['Aluminum Operation Cost Change (Billion CNY)']:<25.2f} {row['Total Cost Change (Billion CNY)']:<20.2f} {row['Emissions Change (Million Tonnes CO2)']:<25.2f}")
    
    # Print baseline information
    print(f"\n=== Baseline Version Information ===")
    print(f"Aluminum cost baseline: {aluminum_baseline_version} (5p)")
    print(f"Power system cost baseline: {power_baseline_version} (non-flexible)")
    print(f"Emissions baseline: {power_baseline_version} (non-flexible)")
    print("Positive values indicate cost savings and emissions reduction, negative values indicate cost increases and emissions increase")

if __name__ == "__main__":
    plot_2050_m_market_costs()
