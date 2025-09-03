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
    
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    # 第一步：收集所有数据并保存到CSV
    all_plot_data = []
    market_levels = ['L', 'M', 'H']
    flexibility_levels = ['L', 'M', 'H', 'N']
    
    logger.info(f"开始收集绘图数据，总共需要处理 {len(flexibility_levels) * len(market_levels)} 个场景组合")
    for flexibility in flexibility_levels:
        for market in market_levels:
            scenario_code = f"{flexibility}M{market}"
            logger.info(f"处理场景 {scenario_code} (F:{flexibility}, D:M, M:{market})")
            if scenario_code in scenarios:
                scenario_data = load_scenario_data(scenarios[scenario_code], file_type)
                
                # 检查是否有数据
                has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
                has_non_flex_data = 'non_flexible' in scenario_data and not scenario_data['non_flexible'].empty
                
                # 如果数据不可用，直接跳过，不打印任何信息
                if not has_100p_data or not has_non_flex_data:
                    continue
                
                # 计算成本差异
                cost_diff = calculate_cost_difference(scenario_data['100p'], scenario_data['non_flexible'])
                
                if cost_diff:
                    # 转换为人民币并排除Total Change
                    for k, v in cost_diff.items():
                        if k != 'Total Change':
                            if not pd.isna(v):
                                value_cny = v * EUR_TO_CNY
                                all_plot_data.append({
                                    'Market': market,
                                    'Flexibility': flexibility,
                                    'Category': k,
                                    'Value (CNY)': value_cny,
                                    'Value (Billion CNY)': value_cny / 1e9,
                                    'Scenario_Code': scenario_code
                                })
                    logger.info(f"场景 {scenario_code}: 添加了 {len([k for k, v in cost_diff.items() if k != 'Total Change' and not pd.isna(v)])} 条数据")
                else:
                    logger.warning(f"场景 {scenario_code}: cost_diff 为空或False")
    
    logger.info(f"数据收集循环完成，all_plot_data 长度: {len(all_plot_data)}")
    
    # 保存汇总的绘图数据
    logger.info(f"收集到的绘图数据条数: {len(all_plot_data)}")
    if not all_plot_data:
        logger.warning("没有找到有效的绘图数据，无法生成图表")
        logger.warning("请检查场景数据是否正确加载")
        return
    
    all_plot_df = pd.DataFrame(all_plot_data)
    summary_plot_csv = plots_dir / f"all_plot_data_{file_type}.csv"
    all_plot_df.to_csv(summary_plot_csv, index=False)
    
    # 第二步：基于CSV数据生成图表
    # 为每个market-flexibility组合创建子图
    fig, axes = plt.subplots(3, 4, figsize=(20, 16), sharey=True)
    
    # 设置子图之间的间距，移除内部边框，调整行间距
    # 使用gridspec来精确控制行间距
    from matplotlib import gridspec
    gs = gridspec.GridSpec(3, 4, figure=fig, height_ratios=[1, 1, 1], hspace=0.1, wspace=0.1)
    
    # 重新分配axes
    axes = []
    for i in range(3):
        row_axes = []
        for j in range(4):
            ax = fig.add_subplot(gs[i, j])
            row_axes.append(ax)
        axes.append(row_axes)
    axes = np.array(axes)
    
    for i, market in enumerate(market_levels):
        for j, flexibility in enumerate(flexibility_levels):
            ax = axes[i, j]
            
            # 移除内部边框，只保留外部边框
            if i < 2:  # 不是最后一行
                ax.spines['bottom'].set_visible(False)
            if j < 3:  # 不是最后一列
                ax.spines['right'].set_visible(False)
            if i > 0:  # 不是第一行
                ax.spines['top'].set_visible(False)
            if j > 0:  # 不是第一列
                ax.spines['left'].set_visible(False)
            
            # 从CSV数据中筛选当前market-flexibility组合的数据
            current_data = all_plot_df[
                (all_plot_df['Market'] == market) & 
                (all_plot_df['Flexibility'] == flexibility)
            ]
            
            if not current_data.empty:
                # 获取该market-flexibility组合下所有分类的数据
                category_data_dict = {}
                for _, row in current_data.iterrows():
                    category_data_dict[row['Category']] = row['Value (CNY)']
                
                if category_data_dict:
                    # 获取所有分类
                    categories = list(category_data_dict.keys())
                    
                    if categories:
                        # 从配置文件读取成本分类颜色
                        config = load_config('config.yaml')
                        category_colors = config.get('cost_category_colors', {}) if config else {}
                        
                        # 定义资源分类的优先级顺序，用于在正负号相同时进行排序
                        category_priority = {
                            "variable cost-non-renewable": 1,
                            "capital-non-renewable": 2,
                            "heating-electrification": 3,
                            "capital-renewable": 4,
                            "transmission lines": 5,
                            "batteries": 6,
                            "long-duration storages": 7,
                            "carbon capture": 8,
                            "synthetic fuels": 9,
                            "carbon management": 10,
                        }
                        
                        # 按照优先级排序分类
                        def sort_key(category):
                            priority = category_priority.get(category, 999)
                            return priority
                        
                        categories = sorted(categories, key=sort_key)
                        
                        # 为每个分类分配颜色，如果不在预定义中则使用默认颜色
                        colors = []
                        for category in categories:
                            if category in category_colors:
                                colors.append(category_colors[category])
                            else:
                                # 使用默认颜色映射
                                colors.append(plt.cm.tab20(len(colors) % 20))
                        
                        # 准备堆叠数据
                        positive_changes = []
                        negative_changes = []
                        
                        for category in categories:
                            value = category_data_dict.get(category, 0)
                            # 调整方向：成本减少（负值）显示在上方，成本增加（正值）显示在下方
                            if value < 0:  # 成本减少，显示在上方
                                positive_changes.append(abs(value))  # 取绝对值
                                negative_changes.append(0)
                            elif value > 0:  # 成本增加，显示在下方
                                positive_changes.append(0)
                                negative_changes.append(value)
                            else:
                                positive_changes.append(0)
                                negative_changes.append(0)
                        
                        # 用于跟踪已经添加到legend的分类
                        added_to_legend = set()
                        
                        # 创建堆叠柱状图
                        x_pos = [0]  # 只有一个柱子
                        width = 0.6  # 柱子宽度
                        
                        # 绘制正值堆叠（成本减少，在横轴上面）
                        bottom_positive = 0
                        for cat_idx, category in enumerate(categories):
                            category_value = positive_changes[cat_idx]
                            if category_value > 0:
                                # 只在第一次遇到该分类时添加到legend
                                label = category if category not in added_to_legend else ""
                                ax.bar(x_pos, [category_value], width, 
                                       bottom=[bottom_positive],
                                       color=colors[cat_idx], 
                                       alpha=0.8,
                                       label=label)
                                bottom_positive += category_value
                                added_to_legend.add(category)
                        
                        # 绘制负值堆叠（成本增加，在横轴下面）
                        bottom_negative = 0
                        for cat_idx, category in enumerate(categories):
                            category_value = negative_changes[cat_idx]
                            if category_value > 0:
                                # 负值部分也添加到legend，但只在第一次遇到该分类时添加
                                label = category if category not in added_to_legend else ""
                                ax.bar(x_pos, [-category_value], width,  # 使用负值
                                       bottom=[bottom_negative],
                                       color=colors[cat_idx], 
                                       alpha=0.8,
                                       label=label)
                                bottom_negative += category_value
                                added_to_legend.add(category)
                        
                        # 设置标签
                        ax.set_xticks(x_pos)
                        ax.set_xticklabels([''], fontsize=14)
                        
                        # 只在最左边的子图显示y轴标签
                        if j == 0:
                            ax.set_ylabel('Cost Change (Billion CNY)', fontsize=14)
                        else:
                            ax.set_ylabel('')
                        
                        # 添加零线
                        ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1.5)
                        
                        # 设置标题，将情景编号放在外面
                        if i == 0:  # 第一行显示market标签
                            ax.set_title(f'Flexibility: {scenario_descriptions[flexibility]}', 
                                       fontsize=14, fontweight='bold', pad=10)
                        if j == 0:  # 第一列显示market标签
                            ax.text(-0.2, 0.5, f'Market: {scenario_descriptions[market]}', 
                                   fontsize=14, fontweight='bold', rotation=90, 
                                   ha='center', va='center', transform=ax.transAxes)
                        
                        # 添加网格
                        ax.grid(True, alpha=0.3, axis='y')
                        
                        # 设置y轴标签为十亿人民币单位，所有子图都显示
                        y_ticks = ax.get_yticks()
                        y_tick_labels = [f'{tick/1e9:.1f}B' for tick in y_ticks]
                        ax.set_yticks(y_ticks)
                        ax.set_yticklabels(y_tick_labels, fontsize=12)
                        
                        # 设置统一的y轴范围
                        ax.set_ylim(-40e9, 100e9)
                    else:
                        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', 
                               transform=ax.transAxes, fontsize=10)
                else:
                    ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', 
                           transform=ax.transAxes, fontsize=10)
            else:
                ax.text(0.5, 0.5, f'No data for\nMarket:{market}, Flexibility:{flexibility}', ha='center', va='center', 
                       transform=ax.transAxes, fontsize=10)
    
    # 创建图例
    if not all_plot_df.empty:
        # 获取所有唯一的分类
        all_categories = all_plot_df['Category'].unique()
        
        # 从配置文件读取成本分类颜色
        config = load_config('config.yaml')
        category_colors = config.get('cost_category_colors', {}) if config else {}
        
        # 创建图例元素
        legend_elements = []
        for category in all_categories:
            if category in category_colors:
                color = category_colors[category]
            else:
                # 使用默认颜色
                color = plt.cm.tab20(len(legend_elements) % 20)
            
            legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=color, 
                                               label=category, alpha=0.8))
        
        # 在图表外部添加统一图例
        fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.5),
                   title='Resource Categories', fontsize=14, title_fontsize=16)
    
    plt.tight_layout()
    
    # 保存图表
    plot_file = plots_dir / f"scenario_comparison_{file_type}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Scenario comparison plot saved to: {plot_file}")
    
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
