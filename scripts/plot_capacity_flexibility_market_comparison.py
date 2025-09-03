#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制多灵活性多市场机会的容量分析图表
竖轴：不同灵活性级别（L、M、H、N）
横轴：不同市场机会（L、M、H）
每个子图显示不同容量比例下的结果
Demand设为M，年份设为2050
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

def get_base_version_from_config(config_path):
    """
    从配置文件中获取基准版本号
    
    Parameters:
    -----------
    config_path : str or Path
        配置文件路径
        
    Returns:
    --------
    str
        基准版本号
    """
    config = load_config(config_path)
    if config is None:
        logger.warning(f"无法加载配置文件 {config_path}，使用默认版本号")
        return None
    
    version = config.get('version', '')
    return version

def find_scenario_results(results_dir, base_version):
    """
    查找所有场景的结果目录，基于基准版本号
    
    Parameters:
    -----------
    results_dir : str
        结果目录路径
    base_version : str
        基准版本号
        
    Returns:
    --------
    dict
        按场景组织的版本信息
    """
    results_path = Path(results_dir)
    if not results_path.exists():
        logger.error(f"结果目录不存在: {results_dir}")
        return {}
    
    # 定义所有可能的场景组合
    flexibility_levels = ['L', 'M', 'H', 'N']  # Low, Mid, High, Non-constrained
    demand_levels = ['M']  # 固定为Mid
    market_levels = ['L', 'M', 'H']  # Low, Mid, High
    year = '2050'
    
    scenarios = {}
    
    # 为每个flexibility-market组合创建100p场景
    for flexibility in flexibility_levels:
        for market in market_levels:
            scenario_code = f"{flexibility}M{market}"  # 固定demand为M
            scenarios[scenario_code] = {}
            
            # 创建100p版本信息
            version_name_100p = f"{base_version}-{scenario_code}-{year}-100p"
            version_dir_100p = results_path / f"version-{version_name_100p}"
            
            scenarios[scenario_code]['100p'] = {
                'version_name': version_name_100p,
                'version_dir': version_dir_100p,
                'year': year,
                'flexibility': flexibility,
                'demand': 'M',
                'market': market,
                'config_type': '100p'
            }
    
    # 为每个场景分配对应的non-flexible基准配置
    for flex in flexibility_levels:
        for mar in market_levels:
            full_scenario_code = f"{flex}M{mar}"
            if full_scenario_code in scenarios:
                # 对于non-flexible情景，使用中等水平的flex和demand，但保持相同的market
                baseline_flexibility = 'M'  # 使用中等灵活性
                baseline_demand = 'M'  # 使用中等需求
                baseline_market = mar  # 保持相同的market
                
                baseline_scenario_code = f"{baseline_flexibility}{baseline_demand}{baseline_market}"
                
                # 创建non-flexible版本信息
                version_name_non_flex = f"{base_version}-{baseline_scenario_code}-{year}-non_flexible"
                version_dir_non_flex = results_path / f"version-{version_name_non_flex}"
                
                scenarios[full_scenario_code]['non_flexible'] = {
                    'version_name': version_name_non_flex,
                    'version_dir': version_dir_non_flex,
                    'year': year,
                    'flexibility': baseline_flexibility,  # 基准使用中等flex
                    'demand': baseline_demand,  # 基准使用中等demand
                    'market': baseline_market,  # 基准使用相同的market
                    'config_type': 'non_flexible'
                }
    
    return scenarios

def load_scenario_data(scenario_info, file_type='costs'):
    """
    加载单个场景的数据
    
    Parameters:
    -----------
    scenario_info : dict
        场景信息
    file_type : str
        文件类型 (costs 或 capacities)
        
    Returns:
    --------
    dict
        包含100p和non_flexible数据的字典
    """
    data = {}
    
    for config_type, info in scenario_info.items():
        version_dir = info['version_dir']
        year = info['year']
        
        # 构建summary目录路径
        summary_dir = version_dir / 'summary' / 'postnetworks' / 'positive'
        
        if summary_dir.exists():
            # 查找对应年份的目录
            year_pattern = f"postnetwork-ll-current+Neighbor-linear2050-{year}"
            year_dir = summary_dir / year_pattern
            
            if year_dir.exists():
                file_path = year_dir / f"{file_type}.csv"
                if file_path.exists():
                    try:
                        df = load_single_csv_file(file_path)
                        if df is not None:
                            data[config_type] = df
                        else:
                            logger.warning(f"无法加载 {config_type} 的 {file_type} 数据，将使用空数据")
                            # 创建空的DataFrame，确保数据结构一致
                            data[config_type] = pd.DataFrame()
                    except Exception as e:
                        logger.warning(f"加载 {config_type} 的 {file_type} 数据时出错: {str(e)}，将使用空数据")
                        # 创建空的DataFrame，确保数据结构一致
                        data[config_type] = pd.DataFrame()
                else:
                    # 创建空的DataFrame，确保数据结构一致
                    data[config_type] = pd.DataFrame()
            else:
                # 创建空的DataFrame，确保数据结构一致
                data[config_type] = pd.DataFrame()
        else:
            # 创建空的DataFrame，确保数据结构一致
            data[config_type] = pd.DataFrame()
    
    return data

def load_single_csv_file(file_path):
    """
    加载单个CSV文件
    
    Parameters:
    -----------
    file_path : Path
        CSV文件路径
        
    Returns:
    --------
    pd.DataFrame or None
        加载的数据
    """
    try:
        # 首先尝试读取文件的前几行来判断结构
        with open(file_path, 'r') as f:
            first_lines = [f.readline().strip() for _ in range(5)]
        
        # 检查是否有多级索引（通过检查前几行是否有空列）
        has_multiindex = any(',' in line and line.split(',')[1] == '' for line in first_lines)
        
        if has_multiindex:
            # 多级索引文件 - 读取时保留所有列
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
        else:
            # 普通单级索引文件
            df = pd.read_csv(file_path, index_col=0)
            # 尝试将所有列转换为数值类型
            for col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        return df
        
    except Exception as e:
        logger.warning(f"加载 {file_path} 时出错: {str(e)}")
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

def calculate_cost_difference(costs_100p, costs_non_flex):
    """
    计算100p和non_flexible之间的成本差异
    
    Parameters:
    -----------
    costs_100p : pd.DataFrame
        100p配置的成本数据
    costs_non_flex : pd.DataFrame
        non_flexible配置的成本数据
        
    Returns:
    --------
    dict
        按成本分类组织的差异数据
    """
    if costs_100p is None or costs_non_flex is None:
        return {}
    
    # 过滤掉aluminum相关的数据
    def filter_aluminum_data(df):
        if df is None or df.empty:
            return df
        
        # 创建过滤后的DataFrame
        filtered_df = df.copy()
        
        # 删除包含aluminum的行
        aluminum_mask = []
        for idx in filtered_df.index:
            if len(idx) >= 3:
                # 检查任何索引级别是否包含aluminum
                has_aluminum = any('aluminum' in str(level).lower() for level in idx)
                aluminum_mask.append(not has_aluminum)
            else:
                aluminum_mask.append(True)
        
        return filtered_df[aluminum_mask]
    
    # 过滤两个数据集
    costs_100p = filter_aluminum_data(costs_100p)
    costs_non_flex = filter_aluminum_data(costs_non_flex)
    
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
    category_changes = {}
    total_change = 0
    
    # 计算每个成本分类的变化量
    for idx in costs_100p.index:
        if len(idx) >= 3:
            component_type, cost_type, carrier = idx[0], idx[1], idx[2]
            
            # 使用分类映射
            category_key = (cost_type, carrier)
            category_name = cost_category_mapping.get(category_key, f"{cost_type} - {carrier}")
            
            if category_name not in category_changes:
                category_changes[category_name] = 0
            
            # 获取两个版本的对应值
            v1_value = costs_100p.loc[idx].iloc[0]
            v2_value = 0
            
            # 在non_flexible版本中查找对应值
            for idx2 in costs_non_flex.index:
                if len(idx2) >= 3 and idx2[0] == component_type and idx2[1] == cost_type and idx2[2] == carrier:
                    v2_value = costs_non_flex.loc[idx2].iloc[0]
                    break
            
            # 处理NaN值
            if pd.isna(v1_value):
                v1_value = 0
            if pd.isna(v2_value):
                v2_value = 0
            
            # 计算变化量（100p - non_flexible，节约为正，增加为负）
            change = v1_value - v2_value
            category_changes[category_name] += change
            total_change += change
    
    # 过滤掉不需要展示的分类
    exclude_categories = {
        'marginal - renewable',
        'marginal - heat pump',
        'marginal - resistive heater',
        'marginal - onwind',
        'marginal - offwind',
        'marginal - solar',
        'marginal - solar thermal',
        'marginal - biomass',
        'marginal - biogas',
        'marginal - H2',
        'marginal - H2 CHP',
    }
    
    filtered_categories = {}
    filtered_total_change = 0  # 只计算过滤后的分类的总变化量
    
    for category, change in category_changes.items():
        # 过滤掉nan相关的分类和不需要展示的分类
        if (category not in exclude_categories and 
            abs(change) > 1e-6 and 
            'nan' not in str(category).lower() and
            not pd.isna(category)):
            filtered_categories[category] = change
            filtered_total_change += change
    
    # 添加过滤后的总变化量
    filtered_categories['Total Change'] = filtered_total_change
    
    return filtered_categories

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

def plot_single_flexibility_market(flexibility, market, base_version, capacity_ratios, results_dir, ax):
    """
    绘制单个灵活性-市场组合的图表
    
    Parameters:
    -----------
    flexibility : str
        灵活性级别 (L, M, H, N)
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
    """
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    year = 2050
    
    # 构建版本名称 - 使用正确的格式
    version_names = []
    config_versions = {}
    
    for ratio in capacity_ratios:
        # 版本号格式: base_version-{flexibility}M{market}-{year}-{ratio}
        version = f"{base_version}-{flexibility}M{market}-{year}-{ratio}"
        version_names.append(version)
        config_versions[ratio] = version
    
    # 基准版本
    aluminum_baseline_version = f"{base_version}-{flexibility}M{market}-{year}-5p"
    power_baseline_version = f"{base_version}-{flexibility}M{market}-{year}-non_flexible"
    
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
        ax.text(0.5, 0.5, f'No data for {flexibility}-{market}', ha='center', va='center', 
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
        market_part = base_version.split('-')[1] if len(base_version.split('-')) > 1 else f'{flexibility}M{market}'
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
    ax.set_xlabel('Aluminum Capacity (MW)', fontsize=12)
    ax.set_ylabel('Cost Savings (Billion CNY)', fontsize=12, color='blue')
    ax2.set_ylabel('Emissions Reduction (Million Tonnes CO2)', fontsize=12, color='red')
    ax.set_title(f'Flexibility: {flexibility}, Market: {market}', fontsize=14, fontweight='bold')
                        
                        # 添加零线
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
    ax2.axhline(y=0, color='red', linestyle='--', alpha=0.5, linewidth=1)
    # 添加网格
    ax.grid(True, alpha=0.3, axis='y')
                        
    # 设置x轴刻度和标签
    ax.set_xticks(x)
    ax.set_xticklabels([f'{cap:.0f}' for cap in capacity_values], fontsize=10)
    
    # 设置y轴标签为十亿人民币单位
    y_ticks = ax.get_yticks()
    y_tick_labels = [f'{tick/1e9:.1f}B' for tick in y_ticks]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels, fontsize=10, color='blue')
    
    # 设置右y轴标签为百万吨CO2单位
    y2_ticks = ax2.get_yticks()
    y2_tick_labels = [f'{tick:.1f}M' for tick in y2_ticks]
    ax2.set_yticks(y2_ticks)
    ax2.set_yticklabels(y2_tick_labels, fontsize=10, color='red')
    
    # 不在这里添加图例，只在总图上添加一个统一的图例

def generate_scenario_plots(scenarios, output_dir, file_type='costs'):
    """
    生成场景对比图表
    
    Parameters:
    -----------
    scenarios : dict
        场景信息
    output_dir : Path
        输出目录
    file_type : str
        文件类型
    """
    # 创建输出目录
    plots_dir = output_dir / "scenario_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # 定义场景代码到描述的映射
    scenario_descriptions = {
        'L': 'Low',
        'M': 'Mid', 
        'H': 'High',
        'N': 'Non-constrained'
    }
    
    # 从主配置文件读取基础版本号
    main_config = load_config('config.yaml')
    if main_config is None:
        logger.error("无法加载主配置文件 config.yaml")
        return
    
    base_version = main_config.get('version', '0814.4H.2')
    logger.info(f"从主配置文件读取到基础版本号: {base_version}")
    
    # 定义容量比例
    capacity_ratios = ['5p', '10p', '20p', '30p', '40p', '50p', '60p', '70p', '80p', '90p', '100p']
    
    # 定义市场机会和灵活性级别
    markets = ['L', 'M', 'H']
    flexibility_levels = ['L', 'M', 'H', 'N']
    
    # 创建子图
    n_markets = len(markets)
    n_flexibility = len(flexibility_levels)
    
    fig, axes = plt.subplots(n_markets, n_flexibility, figsize=(6*n_flexibility, 5*n_markets))
    
    # 如果只有一个市场，确保axes是二维数组
    if n_markets == 1:
        axes = axes.reshape(1, -1)
    
    # 设置子图之间的间距
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    
    # 为每个市场-灵活性组合绘制图表
    for i, market in enumerate(markets):
        for j, flexibility in enumerate(flexibility_levels):
            ax = axes[i, j]
            
            logger.info(f"正在绘制 {market}市场-{flexibility}灵活性 的图表...")
            plot_single_flexibility_market(flexibility, market, base_version, capacity_ratios, 'results', ax)
    
    # 创建统一的图例
    # 从第一个子图获取图例元素
    first_ax = axes[0, 0] if n_markets == 1 else axes[0, 0]
    first_ax2 = first_ax.twinx()
        
        # 创建图例元素
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor='#1f77b4', alpha=0.8, label='Power System Cost Savings'),
        plt.Rectangle((0,0),1,1, facecolor='#ff7f0e', alpha=0.8, label='Aluminum Operation Cost Savings'),
        plt.Line2D([0], [0], color='black', linewidth=3, marker='o', markersize=8, label='Net Cost Savings'),
        plt.Line2D([0], [0], marker='*', color='red', markersize=15, linestyle='', label='Highest Net Savings'),
        plt.Line2D([0], [0], color='red', linewidth=2, marker='o', markersize=6, label='Emissions Reduction')
    ]
    
    # 在总图右侧添加统一图例
    fig.legend(handles=legend_elements, loc='center right', bbox_to_anchor=(1.15, 0.5),
               title='Legend', fontsize=12, title_fontsize=14)
    
    # 添加总标题
    fig.suptitle('Multi-Flexibility Market Opportunity Analysis\n(Demand: M, Year: 2050)\nCost Savings (Positive) & Emissions Reduction (Positive)', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # 添加行标签（市场机会）
    for i, market in enumerate(markets):
        market_desc = {'L': 'Low', 'M': 'Mid', 'H': 'High'}
        fig.text(0.02, 0.8 - i*(0.8/(n_markets-1)) if n_markets > 1 else 0.8, f'Market: {market_desc[market]}', 
                fontsize=14, fontweight='bold', rotation=90, ha='center', va='center')
    
    # 添加列标签（灵活性级别）
    for j, flexibility in enumerate(flexibility_levels):
        flexibility_desc = {'L': 'Low', 'M': 'Mid', 'H': 'High', 'N': 'Non-constrained'}
        fig.text(0.2 + j*(0.6/(n_flexibility-1)) if n_flexibility > 1 else 0.2, 0.95, 
                f'Flexibility: {flexibility_desc[flexibility]}', fontsize=14, fontweight='bold', ha='center')
    
    plt.tight_layout()
    
    # 保存图表
    plot_file = plots_dir / f"flexibility_market_comparison_{file_type}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"灵活性市场对比图表已保存到: {plot_file}")
    
    plt.close()

def validate_scenario_matching(scenarios, file_type='costs'):
    """
    验证场景匹配的正确性
    
    Parameters:
    -----------
    scenarios : dict
        场景信息
    file_type : str
        文件类型
    """
    total_scenarios = len(scenarios)
    valid_scenarios = 0
    invalid_scenarios = 0
    
    for scenario_code, scenario_info in scenarios.items():
        if len(scenario_code) == 3:
            # 检查目录是否存在
            dir_100p = scenario_info['100p']['version_dir']
            dir_non_flex = scenario_info['non_flexible']['version_dir']
            
            if dir_100p.exists() and dir_non_flex.exists():
                valid_scenarios += 1
            else:
                invalid_scenarios += 1
                logger.warning(f"场景 {scenario_code}: 目录缺失")
    
    logger.info(f"场景验证完成: 有效 {valid_scenarios} 个, 无效 {invalid_scenarios} 个")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='可视化多灵活性多市场机会的容量分析结果')
    parser.add_argument('--results-dir', default='results', help='结果目录路径 (默认: results)')
    parser.add_argument('--output', default='results/scenario_analysis', help='输出目录')
    parser.add_argument('--file-type', choices=['costs', 'capacities'], default='costs', 
                       help='分析的文件类型 (默认: costs)')
    parser.add_argument('--config', default='config.yaml', help='配置文件路径 (默认: config.yaml)')
    
    args = parser.parse_args()
    
    # 设置日志级别
    log_level = logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    
    logger.info(f"开始分析多灵活性多市场机会结果，文件类型: {args.file_type}")
    logger.info(f"结果目录: {args.results_dir}")
    logger.info(f"输出目录: {args.output}")
    logger.info(f"配置文件: {args.config}")
    
    # 从配置文件读取基准版本号
    base_version = get_base_version_from_config(args.config)
    if base_version is None:
        logger.error("无法获取基准版本号，请检查配置文件")
        return
    
    # 查找场景结果
    scenarios = find_scenario_results(args.results_dir, base_version)
    
    if not scenarios:
        logger.error("没有找到任何场景结果")
        return
    
    logger.info(f"基于基准版本 {base_version} 构建了 {len(scenarios)} 个场景")
    
    # 创建输出目录
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 验证场景匹配正确性
    logger.info("验证场景匹配正确性...")
    validate_scenario_matching(scenarios, args.file_type)
    
    # 生成场景对比图表
    logger.info("生成场景对比图表...")
    generate_scenario_plots(scenarios, output_path, args.file_type)
    
    logger.info("分析完成！")
    logger.info(f"结果保存在: {output_path}")

if __name__ == "__main__":
    main()
