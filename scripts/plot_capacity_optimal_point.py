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



def find_optimal_points(base_version, capacity_ratios, results_dir='results'):
    """
    找到每个年份-市场组合的最优点
    
    Parameters:
    -----------
    base_version : str
        基础版本号
    capacity_ratios : list
        容量比例列表
    results_dir : str
        结果目录
        
    Returns:
    --------
    list
        包含最优点信息的列表，每个元素为(year, market, capacity, net_value)
    """
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    # 查找可用的年份
    available_years = find_available_years(results_dir, base_version)
    logger.info(f"找到可用的年份: {available_years}")
    
    # 定义市场机会
    markets = ['L', 'M', 'H']
    
    optimal_points = []
    
    for year in available_years:
        for market in markets:
            logger.info(f"正在分析 {year}年-{market}市场 的最优点...")
            
            # 构建版本名称
            version_names = []
            config_versions = {}
            
            for ratio in capacity_ratios:
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
                logger.warning(f"没有找到 {year}年-{market}市场 的数据")
                continue
            
            # 计算净价值
            net_values = []
            capacity_values = []
            
            for ratio in capacity_ratios:
                if ratio in costs_data and 'power' in baseline_data:
                    # 计算电力系统成本变化（成本减少为正方向）
                    current_costs = calculate_cost_categories(costs_data[ratio])
                    baseline_costs = calculate_cost_categories(baseline_data['power'])
                    
                    # 计算总成本变化（排除aluminum相关）
                    power_cost_change = 0
                    for category, value in current_costs.items():
                        if 'aluminum' not in category.lower():
                            baseline_value = baseline_costs.get(category, 0)
                            power_cost_change += (value - baseline_value)
                    
                    # 计算电解铝成本变化
                    aluminum_cost_change = 0
                    if 'aluminum' in baseline_data:
                        aluminum_baseline_costs = calculate_cost_categories(baseline_data['aluminum'])
                        for category, value in current_costs.items():
                            if 'aluminum' in category.lower():
                                baseline_value = aluminum_baseline_costs.get(category, 0)
                                aluminum_cost_change += (value - baseline_value)
                    
                    # 计算净价值（成本减少为正方向，所以取负值）
                    net_value = -(power_cost_change + aluminum_cost_change) * EUR_TO_CNY
                    net_values.append(net_value)
                    
                    # 读取容量比例值
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
                else:
                    net_values.append(0)
                    capacity_values.append(0)
            
            # 找到净价值最大的点
            if net_values:
                max_index = np.argmax(net_values)
                optimal_capacity = capacity_values[max_index]
                optimal_net_value = net_values[max_index]
                
                optimal_points.append({
                    'year': year,
                    'market': market,
                    'capacity': optimal_capacity,
                    'net_value': optimal_net_value / 1e9  # 转换为十亿人民币
                })
                
                logger.info(f"{year}年-{market}市场 最优点: 容量={optimal_capacity:.0f}MW, 净价值={optimal_net_value/1e9:.2f}B CNY")
    
    return optimal_points

def plot_optimal_points_scatter():
    """
    绘制最优点对应的容量和净价值散点图
    """
    # 从主配置文件读取基础版本号
    main_config = load_config('config.yaml')
    if main_config is None:
        logger.error("无法加载主配置文件 config.yaml")
        return
    
    base_version = main_config.get('version', '0814.4H.2')
    logger.info(f"从主配置文件读取到基础版本号: {base_version}")
    
    # 定义容量比例
    capacity_ratios = ['5p', '10p', '20p', '30p', '40p', '50p', '60p', '70p', '80p', '90p', '100p']
    
    # 找到所有最优点
    optimal_points = find_optimal_points(base_version, capacity_ratios, 'results')
    
    if not optimal_points:
        logger.error("没有找到任何最优点数据")
        return
    
    # 创建散点图
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 按年份分组数据
    years = sorted(list(set([point['year'] for point in optimal_points])))
    colors = plt.cm.viridis(np.linspace(0, 1, len(years)))
    year_colors = dict(zip(years, colors))
    
    # 为每个市场机会创建不同的标记
    markets = ['L', 'M', 'H']
    markers = ['o', 's', '^']
    market_markers = dict(zip(markets, markers))
    
    # 绘制散点
    for point in optimal_points:
        year = point['year']
        market = point['market']
        capacity = point['capacity']
        net_value = point['net_value']
        
        ax.scatter(capacity, net_value, 
                  c=[year_colors[year]], 
                  marker=market_markers[market],
                  s=150, alpha=0.8, edgecolors='black', linewidth=1)
    
    # 添加标签
    ax.set_xlabel('铝产能容量 (MW)', fontsize=14, fontweight='bold')
    ax.set_ylabel('净价值 (十亿人民币)', fontsize=14, fontweight='bold')
    ax.set_title('最优点对应的容量和净价值分析', fontsize=16, fontweight='bold')
    
    # 添加网格
    ax.grid(True, alpha=0.3)
    
    # 创建图例
    legend_elements = []
    
    # 年份图例
    for year in years:
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor=year_colors[year], 
                                        markersize=10, label=f'{year}年'))
    
    # 市场图例
    for market in markets:
        market_desc = {'L': '低市场机会', 'M': '中市场机会', 'H': '高市场机会'}
        legend_elements.append(plt.Line2D([0], [0], marker=market_markers[market], 
                                        color='w', markerfacecolor='gray', 
                                        markersize=10, label=market_desc[market]))
    
    # 添加图例
    ax.legend(handles=legend_elements, loc='upper left', fontsize=12)
    
    # 为每个点添加标签
    for point in optimal_points:
        year = point['year']
        market = point['market']
        capacity = point['capacity']
        net_value = point['net_value']
        
        ax.annotate(f'{year}-{market}', 
                   xy=(capacity, net_value),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=10, alpha=0.8)
    
    plt.tight_layout()
    
    # 保存图表
    output_dir = Path('results/optimal_points_analysis')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    plot_file = output_dir / "optimal_points_scatter.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"最优点散点图已保存到: {plot_file}")
    
    # 显示图表
    plt.show()

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='绘制最优点对应的容量和净价值散点图')
    parser.add_argument('--results-dir', default='results', help='结果目录路径 (默认: results)')
    parser.add_argument('--output', default='results/optimal_points_analysis', help='输出目录')
    
    args = parser.parse_args()
    
    logger.info(f"开始分析最优点对应的容量和净价值")
    logger.info(f"结果目录: {args.results_dir}")
    logger.info(f"输出目录: {args.output}")
    
    # 绘制最优点散点图
    plot_optimal_points_scatter()
    
    logger.info("分析完成！")

if __name__ == "__main__":
    main()
