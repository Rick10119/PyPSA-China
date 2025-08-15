#!/usr/bin/env python3
"""
可视化generate_value_test_configs生成的场景结果
生成9个图（3market*3demand），每个图显示不同flexibility场景下的结果
"""

import logging
import os
import sys
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import re
import yaml

logger = logging.getLogger(__name__)

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

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
    # logger.info(f"从配置文件 {config_path} 读取到版本号: {version}")
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
    demand_levels = ['L', 'M', 'H']  # Low, Mid, High
    market_levels = ['L', 'M', 'H']  # Low, Mid, High
    year = '2050'
    
    scenarios = {}
    
    # 为每个flexibility-demand-market组合创建100p场景
    for flexibility in flexibility_levels:
        for demand in demand_levels:
            for market in market_levels:
                scenario_code = f"{flexibility}{demand}{market}"
                scenarios[scenario_code] = {}
                
                # 创建100p版本信息
                version_name_100p = f"{base_version}-{scenario_code}-{year}-100p"
                version_dir_100p = results_path / f"version-{version_name_100p}"
                
                scenarios[scenario_code]['100p'] = {
                    'version_name': version_name_100p,
                    'version_dir': version_dir_100p,
                    'year': year,
                    'flexibility': flexibility,
                    'demand': demand,
                    'market': market,
                    'config_type': '100p'
                }
    
    # 为non-flexible情景，每个market只创建一个配置（使用中等水平的flex和demand）
    # logger.info("为non-flexible情景创建基准配置（每个market一个，使用中等flex和demand）")
    
    # 为每个场景分配对应的non-flexible基准配置
    for flex in flexibility_levels:
        for dem in demand_levels:
            for mar in market_levels:
                full_scenario_code = f"{flex}{dem}{mar}"
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
    
    # logger.info(f"基于基准版本 {base_version} 构建了 {len(scenarios)} 个场景")
    # logger.info(f"100p配置: {len(scenarios)} 个 (4种flex × 3种demand × 3种market)")
    # logger.info(f"non-flexible基准配置: 每个场景使用对应market的中等flex和demand配置")
    # logger.info(f"注意：对于non-flexible情景，相同market的场景共享相同的基准配置")
    
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
            
            # 添加调试信息：显示大额变化
            # if abs(v1_value - v2_value) > 1e9:  # 大于1B的变化
            #     logger.info(f"大额变化: {component_type}-{cost_type}-{carrier}: 100p={v1_value/1e9:.2f}B, non_flex={v2_value/1e9:.2f}B, 差值={(v1_value-v2_value)/1e9:.2f}B")
            
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
    
    # 添加调试信息：显示原始总变化量和过滤后的总变化量
    # logger.info(f"原始总变化量: {total_change/1e9:.2f}B, 过滤后总变化量: {filtered_total_change/1e9:.2f}B")
    
    # 显示被过滤掉的分类信息
    # excluded_changes = {}
    # for category, change in category_changes.items():
    #     if category in exclude_categories and abs(change) > 1e6:  # 只显示大于1M的变化
    #         excluded_changes[category] = change
    
    # if excluded_changes:
    #     logger.info(f"被过滤掉的分类变化量:")
    #     for category, change in sorted(excluded_changes.items(), key=lambda x: abs(x[1]), reverse=True):
    #         logger.info(f"  {category}: {change/1e9:.2f}B")
    
    return filtered_categories

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
    
    # logger.info(f"图表输出目录: {plots_dir}")
    
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
    demand_levels = ['L', 'M', 'H']
    market_levels = ['L', 'M', 'H']
    flexibility_levels = ['L', 'M', 'H', 'N']
    
    # logger.info("正在收集绘图数据...")
    # logger.info("注意：对于non-flexible情景，所有flex-demand组合都使用相同的基准配置")
    # logger.info("基准配置使用中等水平的flex和demand，每个market一个")
    
    for demand in demand_levels:
        for market in market_levels:
            for flex in flexibility_levels:
                scenario_code = f"{flex}{demand}{market}"
                if scenario_code in scenarios:
                    scenario_data = load_scenario_data(scenarios[scenario_code], file_type)
                    
                    # 检查是否有数据
                    has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
                    has_non_flex_data = 'non_flexible' in scenario_data and not scenario_data['non_flexible'].empty
                    
                    # 如果数据不可用，直接跳过，不打印任何信息
                    if not has_100p_data or not has_non_flex_data:
                        continue
                    
                    # 数据可用，打印场景验证信息
                    # logger.info(f"验证场景 {scenario_code}:")
                    # logger.info(f"  版本名称: {scenarios[scenario_code]['100p']['version_name']}")
                    # logger.info(f"  版本目录: {scenarios[scenario_code]['100p']['version_dir']}")
                    # logger.info(f"  non_flexible版本名称: {scenarios[scenario_code]['non_flexible']['version_name']}")
                    # logger.info(f"  non_flexible版本目录: {scenarios[scenario_code]['non_flexible']['version_dir']}")
                    # logger.info(f"  注意：non_flexible基准使用中等flex和demand，适用于所有场景")
                    
                    # 数据可用，继续处理
                    # 添加调试信息：检查数据内容
                    # logger.info(f"场景 {scenario_code} 数据检查:")
                    # logger.info(f"  100p数据行数: {len(scenario_data['100p'])}")
                    # logger.info(f"  non_flexible数据行数: {len(scenario_data['non_flexible'])}")
                    
                    # 计算成本差异
                    cost_diff = calculate_cost_difference(scenario_data['100p'], scenario_data['non_flexible'])
                    
                    # 为所有场景添加关键成本分类的调试信息
                    # if cost_diff:
                    #     key_categories = ['capital - onwind', 'capital - solar', 'capital - battery']
                    #     for category in key_categories:
                    #         if category in cost_diff:
                    #             logger.info(f"  {category}: 差异 = {cost_diff[category]/1e9:.2f}B EUR")
                    if cost_diff:
                        # 转换为人民币并排除Total Change
                        for k, v in cost_diff.items():
                            if k != 'Total Change':
                                if not pd.isna(v):
                                    value_cny = v * EUR_TO_CNY
                                    all_plot_data.append({
                                        'Demand': demand,
                                        'Market': market,
                                        'Flexibility': flex,
                                        'Category': k,
                                        'Value (CNY)': value_cny,
                                        'Value (Billion CNY)': value_cny / 1e9,
                                        'Scenario_Code': scenario_code
                                    })
    
    # 保存汇总的绘图数据
    if all_plot_data:
        all_plot_df = pd.DataFrame(all_plot_data)
        summary_plot_csv = plots_dir / f"all_plot_data_{file_type}.csv"
        all_plot_df.to_csv(summary_plot_csv, index=False)
        # logger.info(f"汇总绘图数据已保存到: {summary_plot_csv}")
        
        # 打印数据统计信息
        # logger.info(f"绘图数据统计:")
        # logger.info(f"  总数据行数: {len(all_plot_df)}")
        # logger.info(f"  唯一分类数: {all_plot_df['Category'].nunique()}")
        # logger.info(f"  唯一场景数: {all_plot_df['Scenario_Code'].nunique()}")
        # logger.info(f"  数值范围: {all_plot_df['Value (Billion CNY)'].min():.2f}B 到 {all_plot_df['Value (Billion CNY)'].max():.2f}B CNY")
    
    # 第二步：基于CSV数据生成图表
    # logger.info("正在基于CSV数据生成图表...")
    
    # 为每个demand-market组合创建子图
    fig, axes = plt.subplots(3, 3, figsize=(20, 16))
    fig.suptitle(f'Scenario Comparison: {file_type.capitalize()} Changes by Flexibility Level', fontsize=16, y=0.98)
    
    for i, demand in enumerate(demand_levels):
        for j, market in enumerate(market_levels):
            ax = axes[i, j]
            
            # 从CSV数据中筛选当前demand-market组合的数据
            current_data = all_plot_df[
                (all_plot_df['Demand'] == demand) & 
                (all_plot_df['Market'] == market)
            ]
            
            if not current_data.empty:
                # 获取该demand-market组合下所有flexibility级别的数据
                flex_data_dict = {}
                for flex in flexibility_levels:
                    flex_data = current_data[current_data['Flexibility'] == flex]
                    if not flex_data.empty:
                        flex_data_dict[flex] = {}
                        for _, row in flex_data.iterrows():
                            flex_data_dict[flex][row['Category']] = row['Value (CNY)']
                
                if flex_data_dict:
                    # 获取所有分类
                    all_categories = set()
                    for flex_data in flex_data_dict.values():
                        all_categories.update(flex_data.keys())
                    
                    if all_categories:
                        categories = list(all_categories)
                        flex_names = list(flex_data_dict.keys())
                        
                        # 创建堆叠柱状图
                        x_pos = np.arange(len(flex_names))
                        width = 0.8
                        
                        # 从配置文件读取成本分类颜色
                        config = load_config('config.yaml')
                        category_colors = config.get('cost_category_colors', {}) if config else {}
                        
                        # 为每个分类分配颜色，如果不在预定义中则使用默认颜色
                        colors = []
                        for category in categories:
                            if category in category_colors:
                                colors.append(category_colors[category])
                            else:
                                # 使用默认颜色映射
                                colors.append(plt.cm.tab20(len(colors) % 20))
                        
                        # 准备堆叠数据 - 为每个flexibility级别创建完整的数组
                        positive_changes = []
                        negative_changes = []
                        
                        for flex in flex_names:
                            flex_data = flex_data_dict.get(flex, {})
                            flex_positive = []
                            flex_negative = []
                            
                            for category in categories:
                                value = flex_data.get(category, 0)
                                # 调整方向：成本减少（负值）显示在上方，成本增加（正值）显示在下方
                                if value < 0:  # 成本减少，显示在上方
                                    flex_positive.append(abs(value))  # 取绝对值
                                    flex_negative.append(0)
                                elif value > 0:  # 成本增加，显示在下方
                                    flex_positive.append(0)
                                    flex_negative.append(value)
                                else:
                                    flex_positive.append(0)
                                    flex_negative.append(0)
                            
                            positive_changes.append(flex_positive)
                            negative_changes.append(flex_negative)
                        
                        # 用于跟踪已经添加到legend的分类
                        added_to_legend = set()
                        
                        # 绘制正值堆叠（成本减少，在横轴上面）
                        bottom_positive = np.zeros(len(x_pos))
                        for cat_idx, category in enumerate(categories):
                            category_values = [changes[cat_idx] for changes in positive_changes]
                            if any(val > 0 for val in category_values):
                                # 只在第一次遇到该分类时添加到legend
                                label = category if category not in added_to_legend else ""
                                ax.bar(x_pos, category_values, width, 
                                       bottom=bottom_positive,
                                       color=colors[cat_idx], 
                                       alpha=0.8,
                                       label=label)
                                bottom_positive += np.array(category_values)
                                added_to_legend.add(category)
                        
                        # 绘制负值堆叠（成本增加，在横轴下面）
                        bottom_negative = np.zeros(len(x_pos))
                        for cat_idx, category in enumerate(categories):
                            category_values = [changes[cat_idx] for changes in negative_changes]
                            if any(val > 0 for val in category_values):
                                # 负值部分也添加到legend，但只在第一次遇到该分类时添加
                                label = category if category not in added_to_legend else ""
                                ax.bar(x_pos, -np.array(category_values), width,  # 使用负值
                                       bottom=bottom_negative,
                                       color=colors[cat_idx], 
                                       alpha=0.8,
                                       label=label)
                                bottom_negative += np.array(category_values)
                                added_to_legend.add(category)
                        
                        # 设置标签
                        ax.set_xticks(x_pos)
                        ax.set_xticklabels([f'{flex}' for flex in flex_names], fontsize=10)
                        ax.set_ylabel('Cost Change (Billion CNY)', fontsize=9)
                        
                        # 添加零线
                        ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1.5)
                        
                        # 设置标题
                        ax.set_title(f'Demand: {scenario_descriptions[demand]}, Market: {scenario_descriptions[market]}', 
                                   fontsize=10, fontweight='bold')
                        
                        # 添加网格
                        ax.grid(True, alpha=0.3, axis='y')
                        
                        # 设置y轴标签为十亿人民币单位
                        y_ticks = ax.get_yticks()
                        y_tick_labels = [f'{tick/1e9:.1f}B' for tick in y_ticks]
                        ax.set_yticks(y_ticks)
                        ax.set_yticklabels(y_tick_labels, fontsize=8)
                        
                        # 设置统一的y轴范围，最大值为40十亿人民币
                        ax.set_ylim(-40e9, 80e9)
                    else:
                        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', 
                               transform=ax.transAxes, fontsize=10)
                else:
                    ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', 
                           transform=ax.transAxes, fontsize=10)
            else:
                ax.text(0.5, 0.5, f'No data for\nDemand:{demand}, Market:{market}', ha='center', va='center', 
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
                   title='Resource Categories', fontsize=10)
        
        # logger.info(f"图例包含 {len(legend_elements)} 个分类")
    
    plt.tight_layout()
    # plt.show()
    
    # 保存图表
    plot_file = plots_dir / f"scenario_comparison_{file_type}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    # logger.info(f"Scenario comparison plot saved to: {plot_file}")
    
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
    # logger.info("=== 验证场景匹配正确性 ===")
    
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
    
    # logger.info(f"场景验证完成: 有效 {valid_scenarios} 个, 无效 {invalid_scenarios} 个")
    # logger.info("")

def generate_summary_table(scenarios, output_dir, file_type='costs'):
    """
    生成场景对比摘要表格
    
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
    tables_dir = output_dir / "summary_tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    # 准备表格数据
    table_data = []
    
    for scenario_code, scenario_info in scenarios.items():
        if len(scenario_code) == 3:  # 确保是有效的3位数场景代码
            flexibility, demand, market = scenario_code[0], scenario_code[1], scenario_code[2]
            
            # 加载场景数据
            scenario_data = load_scenario_data(scenario_info, file_type)
            
            # 检查是否有数据
            has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
            has_non_flex_data = 'non_flexible' in scenario_data and not scenario_data['non_flexible'].empty
            
            # 如果数据不可用，直接跳过，不打印任何信息，也不添加到表格中
            if not has_100p_data or not has_non_flex_data:
                continue
            
            # 数据可用，继续处理
            # 添加调试信息：验证数据匹配
            # logger.info(f"场景 {scenario_code}: 100p数据行数={len(scenario_data['100p'])}, non_flexible数据行数={len(scenario_data['non_flexible'])}")
            
            # 计算成本差异
            cost_diff = calculate_cost_difference(scenario_data['100p'], scenario_data['non_flexible'])
            
            if cost_diff and 'Total Change' in cost_diff:
                # 注意：cost_diff中的数据已经是欧元单位，需要转换为人民币
                total_change = cost_diff['Total Change'] * EUR_TO_CNY
                
                                    # 添加调试信息：显示总变化量
                    # logger.info(f"场景 {scenario_code}: 过滤后总成本变化 = {total_change/1e9:.2f}B CNY (原始欧元值: {cost_diff['Total Change']/1e9:.2f}B EUR)")
                
                # 计算各分类的变化量
                category_changes = {}
                for category, change in cost_diff.items():
                    if category != 'Total Change':
                        category_changes[category] = change * EUR_TO_CNY
                
                # 找出变化最大的前3个分类
                top_changes = sorted(category_changes.items(), key=lambda x: abs(x[1]), reverse=True)[:3]
                
                row = {
                    'Scenario': scenario_code,
                    'Flexibility': flexibility,
                    'Demand': demand,
                    'Market': market,
                    'Top Category 1': f"{top_changes[0][0]}: {top_changes[0][1]/1e9:.2f}B" if len(top_changes) > 0 else "N/A",
                    'Top Category 2': f"{top_changes[1][0]}: {top_changes[1][1]/1e9:.2f}B" if len(top_changes) > 1 else "N/A",
                    'Top Category 3': f"{top_changes[2][0]}: {top_changes[2][1]/1e9:.2f}B" if len(top_changes) > 2 else "N/A"
                }
                table_data.append(row)
            else:
                row = {
                    'Scenario': scenario_code,
                    'Flexibility': flexibility,
                    'Demand': demand,
                    'Market': market,
                    'Top Category 1': 'N/A',
                    'Top Category 2': 'N/A',
                    'Top Category 3': 'N/A'
                }
                table_data.append(row)
    
    if table_data:
        # 创建DataFrame并排序
        df = pd.DataFrame(table_data)
        df = df.sort_values(['Demand', 'Market'])
        
        # 保存表格
        table_file = tables_dir / f"scenario_summary_{file_type}.csv"
        df.to_csv(table_file, index=False)
        # logger.info(f"Summary table saved to: {table_file}")
        
        # 打印表格
        print(f"\n=== {file_type.capitalize()} Summary Table ===")
        print(df.to_string(index=False))
    else:
        logger.warning("No data available for summary table")

def generate_raw_cost_comparison_table(scenarios, output_dir, file_type='costs'):
    """
    生成原始成本对比表格，包含每个情景下100p和non-flexible的原始成本数据
    
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
    tables_dir = output_dir / "raw_cost_tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    
    # 准备表格数据
    table_data = []
    
    for scenario_code, scenario_info in scenarios.items():
        if len(scenario_code) == 3:  # 确保是有效的3位数场景代码
            flexibility, demand, market = scenario_code[0], scenario_code[1], scenario_code[2]
            
            # 加载场景数据
            scenario_data = load_scenario_data(scenario_info, file_type)
            
            # 检查是否有数据
            has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
            has_non_flex_data = 'non_flexible' in scenario_data and not scenario_data['non_flexible'].empty
            
            # 如果数据不可用，直接跳过
            if not has_100p_data or not has_non_flex_data:
                continue
            
            # 数据可用，处理原始成本数据
            costs_100p = scenario_data['100p']
            costs_non_flex = scenario_data['non_flexible']
            
            # 收集所有唯一的成本分类
            all_categories = set()
            
            # 从100p数据中收集分类
            for idx in costs_100p.index:
                if len(idx) >= 3:
                    component_type, cost_type, carrier = idx[0], idx[1], idx[2]
                    category_name = f"{cost_type} - {carrier}"
                    all_categories.add(category_name)
            
            # 从non_flexible数据中收集分类
            for idx in costs_non_flex.index:
                if len(idx) >= 3:
                    component_type, cost_type, carrier = idx[0], idx[1], idx[2]
                    category_name = f"{cost_type} - {carrier}"
                    all_categories.add(category_name)
            
            # 为每个成本分类创建一行数据
            for category in sorted(all_categories):
                cost_type, carrier = category.split(" - ", 1)
                
                # 获取100p数据
                value_100p = 0
                for idx in costs_100p.index:
                    if len(idx) >= 3 and idx[1] == cost_type and idx[2] == carrier:
                        value_100p = costs_100p.loc[idx].iloc[0]
                        if pd.isna(value_100p):
                            value_100p = 0
                        break
                
                # 获取non_flexible数据
                value_non_flex = 0
                for idx in costs_non_flex.index:
                    if len(idx) >= 3 and idx[1] == cost_type and idx[2] == carrier:
                        value_non_flex = costs_non_flex.loc[idx].iloc[0]
                        if pd.isna(value_non_flex):
                            value_non_flex = 0
                        break
                
                # 计算差异
                difference = value_100p - value_non_flex
                
                # 添加到表格数据
                row = {
                    'Scenario': scenario_code,
                    'Flexibility': flexibility,
                    'Demand': demand,
                    'Market': market,
                    'Category': category,
                    'Cost_Type': cost_type,
                    'Carrier': carrier,
                    'Value_100p_EUR': value_100p,
                    'Value_100p_Billion_EUR': value_100p / 1e9,
                    'Value_Non_Flex_EUR': value_non_flex,
                    'Value_Non_Flex_Billion_EUR': value_non_flex / 1e9,
                    'Difference_EUR': difference,
                    'Difference_Billion_EUR': difference / 1e9,
                    'Difference_Percent': (difference / value_non_flex * 100) if value_non_flex != 0 else 0
                }
                table_data.append(row)
    
    if table_data:
        # 创建DataFrame并排序
        df = pd.DataFrame(table_data)
        df = df.sort_values(['Scenario', 'Category'])
        
        # 保存表格
        table_file = tables_dir / f"raw_cost_comparison_{file_type}.csv"
        df.to_csv(table_file, index=False)
        # logger.info(f"Raw cost comparison table saved to: {table_file}")
        
        # 创建汇总表格（按场景汇总）
        summary_data = []
        for scenario in df['Scenario'].unique():
            scenario_df = df[df['Scenario'] == scenario]
            
            # 计算总成本
            total_100p = scenario_df['Value_100p_EUR'].sum()
            total_non_flex = scenario_df['Value_Non_Flex_EUR'].sum()
            total_diff = total_100p - total_non_flex
            
            # 找出变化最大的前5个分类
            top_changes = scenario_df.nlargest(5, 'Difference_Billion_EUR')[['Category', 'Difference_Billion_EUR']]
            top_changes_str = "; ".join([f"{row['Category']}: {row['Difference_Billion_EUR']:.2f}B" for _, row in top_changes.iterrows()])
            
            summary_row = {
                'Scenario': scenario,
                'Total_100p_Billion_EUR': total_100p / 1e9,
                'Total_Non_Flex_Billion_EUR': total_non_flex / 1e9,
                'Total_Difference_Billion_EUR': total_diff / 1e9,
                'Total_Difference_Percent': (total_diff / total_non_flex * 100) if total_non_flex != 0 else 0,
                'Top_5_Changes': top_changes_str
            }
            summary_data.append(summary_row)
        
        # 保存汇总表格
        summary_df = pd.DataFrame(summary_data)
        summary_df = summary_df.sort_values('Scenario')
        summary_file = tables_dir / f"raw_cost_summary_{file_type}.csv"
        summary_df.to_csv(summary_file, index=False)
        # logger.info(f"Raw cost summary table saved to: {summary_file}")
    else:
        # logger.warning("No data available for raw cost comparison table")
        pass

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='可视化generate_value_test_configs生成的场景结果')
    parser.add_argument('--results-dir', default='results', help='结果目录路径 (默认: results)')
    parser.add_argument('--output', default='results/scenario_analysis', help='输出目录')
    parser.add_argument('--file-type', choices=['costs', 'capacities'], default='costs', 
                       help='分析的文件类型 (默认: costs)')

    parser.add_argument('--config', default='config.yaml', help='配置文件路径 (默认: config.yaml)')
    
    args = parser.parse_args()
    
    # 设置日志级别
    log_level = logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    
    # logger.info(f"开始分析场景结果，文件类型: {args.file_type}")
    # logger.info(f"结果目录: {args.results_dir}")
    # logger.info(f"输出目录: {args.output}")
    # logger.info(f"配置文件: {args.config}")
    
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
    
    # logger.info(f"基于基准版本 {base_version} 构建了 {len(scenarios)} 个场景")
    
    # 创建输出目录
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 验证场景匹配正确性
    # logger.info("验证场景匹配正确性...")
    validate_scenario_matching(scenarios, args.file_type)
    
    # 生成场景对比图表
    # logger.info("生成场景对比图表...")
    generate_scenario_plots(scenarios, output_path, args.file_type)
    
    # 生成摘要表格
    # logger.info("生成摘要表格...")
    generate_summary_table(scenarios, output_path, args.file_type)
    
    # 生成原始成本对比表格
    # logger.info("生成原始成本对比表格...")
    generate_raw_cost_comparison_table(scenarios, output_path, args.file_type)
    
    # 打印数据完整性统计
    # print_data_completeness_stats(scenarios, args.file_type)
    
    # logger.info("分析完成！")
    # logger.info(f"结果保存在: {output_path}")

def print_data_completeness_stats(scenarios, file_type):
    """
    打印数据完整性统计
    
    Parameters:
    -----------
    scenarios : dict
        场景信息
    file_type : str
        文件类型
    """
    # print(f"\n=== 数据完整性统计 ({file_type}) ===")
    
    # total_scenarios = len(scenarios)
    # complete_scenarios = 0
    # missing_100p = 0
    # missing_non_flex = 0
    # missing_both = 0
    
    # for scenario_code, scenario_info in scenarios.items():
    #     scenario_data = load_scenario_data(scenario_info, file_type)
        
    #     has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
    #     has_non_flex_data = 'non_flexible' in scenario_info and not scenario_data['non_flexible'].empty
        
    #     if has_100p_data and has_non_flex_data:
    #         complete_scenarios += 1
    #     elif not has_100p_data and not has_non_flex_data:
    #         missing_both += 1
    #     elif not has_100p_data:
    #         missing_100p += 1
    #     else:
    #         missing_non_flex += 1
    
    # print(f"总场景数: {total_scenarios}")
    # print(f"数据完整: {complete_scenarios} ({complete_scenarios/total_scenarios*100:.1f}%)")
    # print(f"100p数据缺失: {missing_100p} ({missing_100p/total_scenarios*100:.1f}%)")
    # print(f"non_flexible数据缺失: {missing_non_flex} ({missing_non_flex/total_scenarios*100:.1f}%)")
    # print(f"两者都缺失: {missing_both} ({missing_both/total_scenarios*100:.1f}%)")
    
    # # 按flexibility-demand-market组合显示详细统计
    # print(f"\n=== 按场景组合的数据完整性 ===")
    # flexibility_levels = ['L', 'M', 'H', 'N']
    # demand_levels = ['L', 'M', 'H']
    # market_levels = ['L', 'M', 'H']
    
    # for flexibility in flexibility_levels:
    #     for demand in demand_levels:
    #         for market in market_levels:
    #             scenario_code = f"{flexibility}{demand}{market}"
    #             if scenario_code in scenarios:
    #                 scenario_data = load_scenario_data(scenarios[scenario_code], file_type)
                    
    #                 has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
    #                 has_non_flex_data = 'non_flexible' in scenario_data and not scenario_data['non_flexible'].empty
                    
    #                 if has_100p_data and has_non_flex_data:
    #                     status = "✓ Complete"
    #                 elif not has_100p_data and not has_non_flex_data:
    #                     status = "✗ Both missing"
    #                 elif not has_100p_data:
    #                     status = "✗ 100p missing"
    #                 else:
    #                     status = "✗ Non_flex missing"
                    
    #                 print(f"  {scenario_code} (F:{flexibility}, D:{demand}, M:{market}): {status}")

if __name__ == "__main__":
    main()
