#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制MMM-2050情景下的成本分析图表
横轴：不同电解铝容量比例（5p-100p）
纵轴：成本节约（十亿人民币）和碳排放减少（百万吨CO2）
显示电力系统成本节约、电解铝运行成本变化、净成本节约和碳排放减少
需求设为M，灵活性设为M，市场机会设为M
成本减少为正方向，碳排放减少为正方向
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
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'sans-serif']
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
        logger.error(f"Error loading config file {config_path}: {str(e)}")
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
        logger.warning("No year data found, defaulting to 2050")
    
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
            logger.warning(f"File does not exist: {file_path}")
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
        logger.error(f"Error loading data: {str(e)}")
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

def save_plot_data_to_csv(plot_data, output_dir, year, market):
    """
    将画图数据保存为CSV文件
    
    Parameters:
    -----------
    plot_data : dict
        包含所有计算结果的字典
    output_dir : str or Path
        输出目录路径
    year : int
        年份
    market : str
        市场机会级别
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 创建详细数据表格
    detailed_data = []
    for i, ratio in enumerate(plot_data['capacity_ratios']):
        detailed_data.append({
            'Capacity_Ratio': ratio,
            'Capacity_Value_Mt': plot_data['capacity_values'][i],
            'Power_Cost_Changes_Billion_CNY': plot_data['power_cost_changes'][i] / 1e9,
            'Aluminum_Cost_Changes_Billion_CNY': plot_data['aluminum_cost_changes'][i] / 1e9,
            'Net_Cost_Savings_Billion_CNY': plot_data['net_cost_savings'][i] / 1e9,
            'Emissions_Changes_Million_Tonnes_CO2': plot_data['emissions_changes'][i],
            'Is_Max_Savings': i == plot_data['max_saving_index']
        })
    
    # 保存详细数据
    detailed_df = pd.DataFrame(detailed_data)
    detailed_file = output_path / f"mmm_{year}_{market}_detailed_data.csv"
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
        'Total_Net_Cost_Savings_Billion_CNY': [sum(plot_data['net_cost_savings']) / 1e9],
        'Total_Emissions_Reduction_Million_Tonnes_CO2': [sum(plot_data['emissions_changes'])]
    }
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = output_path / f"mmm_{year}_{market}_summary.csv"
    summary_df.to_csv(summary_file, index=False, encoding='utf-8')
    logger.info(f"Summary data saved to: {summary_file}")
    
    # 创建成本分类数据表格（如果可用）
    try:
        # 这里可以添加更详细的成本分类数据保存
        cost_breakdown_file = output_path / f"mmm_{year}_{market}_cost_breakdown.csv"
        # 暂时创建一个空的成本分解文件，后续可以扩展
        pd.DataFrame({'Note': ['Cost breakdown data will be added in future versions']}).to_csv(
            cost_breakdown_file, index=False, encoding='utf-8')
        logger.info(f"Cost breakdown placeholder saved to: {cost_breakdown_file}")
    except Exception as e:
        logger.warning(f"Could not save cost breakdown data: {str(e)}")

def plot_single_year_market(year, market, base_version, capacity_ratios, results_dir, ax, save_data=True, output_dir=None):
    """
    绘制单个年份-市场组合的图表
    
    Parameters:
    -----------
    year : int
        年份
    market : str
        市场机会级别 (L, M, H)
    base_version : str
        基础版本号
    capacity_ratios : list
        容量比例列表
    results_dir : str
        结果目录
    ax : matplotlib.axes.Axes
        子图对象
    save_data : bool
        是否保存数据到CSV文件
    output_dir : str or Path
        输出目录路径
        
    Returns:
    --------
    dict
        包含所有计算结果的字典
    """
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    # 构建版本名称 - 使用正确的格式
    version_names = []
    config_versions = {}
    
    for ratio in capacity_ratios:
        # 版本号格式: base_version-MM{market}-{year}-{ratio}
        version = f"{base_version}-MM{market}-{year}-{ratio}"
        version_names.append(version)
        config_versions[ratio] = version
    
    # 基准版本
    aluminum_baseline_version = f"{base_version}-MM{market}-{year}-5p"
    power_baseline_version = f"{base_version}-MM{market}-{year}-non_flexible"
    
    # 收集数据
    costs_data = {}
    baseline_data = {}
    
    # 加载基准版本数据
    aluminum_baseline = load_costs_data(aluminum_baseline_version, year, results_dir)
    if aluminum_baseline is not None:
        baseline_data['aluminum'] = aluminum_baseline
    
    power_baseline = load_costs_data(power_baseline_version, year, results_dir)
    if power_baseline is not None:
        baseline_data['power'] = power_baseline
    
    # 加载各容量比例的数据
    for ratio in capacity_ratios:
        version_name = config_versions[ratio]
        costs = load_costs_data(version_name, year, results_dir)
        if costs is not None:
            costs_data[ratio] = costs
    
    if not costs_data or not baseline_data:
        ax.text(0.5, 0.5, f'No data for {year}-{market}', ha='center', va='center', 
               transform=ax.transAxes, fontsize=12)
        return
    
    # 计算成本变化（成本减少为正方向）
    power_cost_changes = []
    aluminum_cost_changes = []
    emissions_changes = []
    capacity_values = []
    
    for ratio in capacity_ratios:
        if ratio in costs_data and 'power' in baseline_data:
            # 计算电力系统成本变化（成本减少为正方向）
            current_costs = calculate_cost_categories(costs_data[ratio])
            baseline_costs = calculate_cost_categories(baseline_data['power'])
            
            # 计算总成本变化（排除aluminum相关）
            total_change = 0
            for category, value in current_costs.items():
                if 'aluminum' not in category.lower():
                    baseline_value = baseline_costs.get(category, 0)
                    total_change += (value - baseline_value)
            
            # 成本减少为正方向，所以取负值
            power_cost_changes.append(-total_change * EUR_TO_CNY)  # 转换为人民币，成本减少为正
        else:
            power_cost_changes.append(0)
        
        if ratio in costs_data and 'aluminum' in baseline_data:
            # 计算电解铝成本变化（成本减少为正方向）
            current_costs = calculate_cost_categories(costs_data[ratio])
            baseline_costs = calculate_cost_categories(baseline_data['aluminum'])
            
            # 计算aluminum相关成本变化
            aluminum_change = 0
            for category, value in current_costs.items():
                if 'aluminum' in category.lower():
                    baseline_value = baseline_costs.get(category, 0)
                    aluminum_change += (value - baseline_value)
            
            # 成本减少为正方向，所以取负值
            aluminum_cost_changes.append(-aluminum_change * EUR_TO_CNY)  # 转换为人民币，成本减少为正
        else:
            aluminum_cost_changes.append(0)
        
        # 计算碳排放变化（碳排放减少为正方向）
        if ratio in costs_data and 'power' in baseline_data:
            current_emissions = calculate_total_emissions_from_costs(costs_data[ratio])
            baseline_emissions_total = calculate_total_emissions_from_costs(baseline_data['power'])
            
            # 碳排放减少为正方向，所以取负值
            emissions_change = -(current_emissions - baseline_emissions_total)  # 碳排放减少为正
            # 转换为百万吨CO2
            emissions_changes.append(emissions_change / 1e6)
        else:
            emissions_changes.append(0)
        
        # 读取容量比例值
        # 从base_version中提取市场机会部分，格式：0815.1H.1-MML-2030-5p -> MML
        market_part = base_version.split('-')[1] if len(base_version.split('-')) > 1 else 'MML'
        config_file = f"configs/config_{market_part}_{year}_{ratio}.yaml"
        config = load_config(config_file)
        if config is not None:
            capacity_ratio = config.get('aluminum_capacity_ratio', 1.0)
            if 'aluminum' in config and 'capacity_ratio' in config['aluminum']:
                capacity_ratio = config['aluminum']['capacity_ratio']
            
            # 计算实际容量 (4500 * capacity ratio)
            actual_capacity = 4500 * capacity_ratio
            capacity_values.append(actual_capacity)
        else:
            # 如果配置文件不存在，使用默认值
            default_ratio = float(ratio.replace('p', '')) / 100.0
            default_capacity = 4500 * default_ratio
            capacity_values.append(default_capacity)
    
    # 创建双y轴图表
    ax2 = ax.twinx()
    
    # 创建图表
    x = capacity_values
    bar_width = 150  # 减少柱子宽度
    
    # 为电力系统成本和电解铝成本创建稍微错开的位置
    x_power = [pos - bar_width/6 for pos in x]  # 电力系统成本稍微向左
    x_aluminum = [pos + bar_width/6 for pos in x]  # 电解铝成本稍微向右
    
    # 绘制电力系统成本节约（底部）
    bars1 = ax.bar(x_power, power_cost_changes, bar_width*0.8, color='#1f77b4', alpha=0.8, 
                   label='Power System Cost Savings')
    
    # 绘制电解铝运行成本节约（从电力系统成本顶部往上堆叠，稍微错开）
    bars2 = ax.bar(x_aluminum, aluminum_cost_changes, bar_width*0.8, bottom=power_cost_changes, 
                   color='#ff7f0e', alpha=0.8, label='Aluminum Operation Cost Increase')
    
    # 计算净值（总成本节约）
    net_cost_savings = [power_cost_changes[i] + aluminum_cost_changes[i] for i in range(len(capacity_values))]
    
    # 绘制净值曲线（黑色线，使用原始x轴位置）
    ax.plot(x, net_cost_savings, 'k-', linewidth=3, label='Net Cost Savings', marker='o', markersize=8, zorder=20)
    
    # 找出净值最大的位置
    max_saving_index = np.argmax(net_cost_savings)
    max_saving_value = net_cost_savings[max_saving_index]
    max_saving_capacity = capacity_values[max_saving_index]
    
    # 用星号标出净值最大处
    ax.plot(max_saving_capacity, max_saving_value, 'r*', markersize=15, 
            label=f'Highest Net Savings: {max_saving_value/1e9:.1f}B CNY', zorder=30)
    
    # 为净值最大点添加数值标签
    ax.annotate(f'{max_saving_value/1e9:.1f}B',
                xy=(max_saving_capacity, max_saving_value),
                xytext=(0, 20),
                textcoords="offset points",
                ha='center', va='bottom', 
                fontsize=12, weight='bold', color='red',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # 绘制碳排放变化（右y轴，使用原始x轴位置）
    line1 = ax2.plot(x, emissions_changes, linewidth=2, marker='o', 
                     markersize=6, label='Emissions Reduction', color='red')
    
    # 设置标签
    ax.set_xlabel('Aluminum Smelting Capacity (Mt)', fontsize=12)
    ax.set_ylabel('Cost Savings (Billion CNY)', fontsize=12, color='blue')
    ax2.set_ylabel('Emissions Reduction (Million Tonnes CO2)', fontsize=12, color='red')
    # ax.set_title(f'Year {year}, Market {market}', fontsize=14, fontweight='bold')
    
    # 添加零线
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
    ax2.axhline(y=0, color='red', linestyle='--', alpha=0.5, linewidth=1)
    
    # 添加网格
    ax.grid(True, alpha=0.3, axis='y')
    
    # 设置x轴刻度和标签
    ax.set_xticks(x)
    ax.set_xticklabels([f'{cap/100:.0f}' for cap in capacity_values], fontsize=12)
    
    # 设置y轴标签为十亿人民币单位
    y_ticks = ax.get_yticks()
    y_tick_labels = [f'{tick/1e9:.0f}' for tick in y_ticks]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels, fontsize=12, color='blue')
    
    # 设置右y轴标签为百万吨CO2单位
    y2_ticks = ax2.get_yticks()
    y2_tick_labels = [f'{tick:.0f}' for tick in y2_ticks]
    ax2.set_yticks(y2_ticks)
    ax2.set_yticklabels(y2_tick_labels, fontsize=12, color='red')
    
    # 不在这里添加图例，只在总图上添加一个统一的图例
    
    # 准备返回的数据
    plot_data = {
        'capacity_values': capacity_values,
        'power_cost_changes': power_cost_changes,
        'aluminum_cost_changes': aluminum_cost_changes,
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
        save_plot_data_to_csv(plot_data, output_dir, year, market)
    
    return plot_data

def plot_mmm_2050_analysis():
    """
    绘制MMM-2050情景下的成本分析图表
    """
    # 从主配置文件读取基础版本号
    main_config = load_config('config.yaml')
    if main_config is None:
        logger.error("Unable to load main config file config.yaml")
        return
    
    base_version = main_config.get('version', '0815.1H.1')
    logger.info(f"Base version read from main config file: {base_version}")
    
    # 定义容量比例
    capacity_ratios = ['5p', '10p', '20p', '30p', '40p', '50p', '60p', '70p', '80p', '90p', '100p']
    
    # 固定为MMM-2050情景
    year = 2050
    market = 'M'
    
    # 创建单个图表
    fig, ax = plt.subplots(1, 1, figsize=(10, 8.5))
    
    logger.info(f"Plotting MMM-2050 scenario chart...")
    
    # 设置输出目录
    output_dir = Path('results/mmm_2050_analysis')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 调用绘图函数并保存数据
    plot_data = plot_single_year_market(year, market, base_version, capacity_ratios, 'results', ax, 
                                       save_data=True, output_dir=output_dir)
    
    # 创建图例
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor='#1f77b4', alpha=0.8, label='Power System Cost Savings'),
        plt.Rectangle((0,0),1,1, facecolor='#ff7f0e', alpha=0.8, label='Aluminum Operation Cost Increase'),
        plt.Line2D([0], [0], color='black', linewidth=3, marker='o', markersize=8, label='Net Cost Savings'),
        plt.Line2D([0], [0], marker='*', color='red', markersize=15, linestyle='', label='Highest Net Savings'),
        plt.Line2D([0], [0], color='red', linewidth=2, marker='o', markersize=6, label='Emissions Reduction')
    ]
    
    # 添加图例
    ax.legend(handles=legend_elements, loc='lower left', fontsize=15)
    
    # # 添加总标题
    # fig.suptitle('MMM-2050 Scenario Analysis\n(Demand: M, Flexibility: M)\nCost Savings (Positive) & Emissions Reduction (Positive)', 
    #              fontsize=16, fontweight='bold', y=0.95)
    
    plt.tight_layout()
    
    # 保存图表（输出目录已在前面创建）
    plot_file = output_dir / "mmm_2050_analysis.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"MMM-2050 scenario analysis chart saved to: {plot_file}")
    
    # 打印数据保存信息
    logger.info(f"Data files saved to: {output_dir}")
    logger.info(f"- Detailed data: mmm_{year}_{market}_detailed_data.csv")
    logger.info(f"- Summary data: mmm_{year}_{market}_summary.csv")
    logger.info(f"- Cost breakdown: mmm_{year}_{market}_cost_breakdown.csv")
    
    # plt.show()

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='Plot cost analysis chart for MMM-2050 scenario')
    parser.add_argument('--results-dir', default='results', help='Results directory path (default: results)')
    parser.add_argument('--output', default='results/mmm_2050_analysis', help='Output directory')
    
    args = parser.parse_args()
    
    logger.info(f"Starting MMM-2050 scenario analysis")
    logger.info(f"Results directory: {args.results_dir}")
    logger.info(f"Output directory: {args.output}")
    
    # Plot MMM-2050 scenario analysis chart
    plot_mmm_2050_analysis()
    
    logger.info("Analysis completed!")

if __name__ == "__main__":
    main()
