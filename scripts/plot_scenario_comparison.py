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
    logger.info(f"从配置文件 {config_path} 读取到版本号: {version}")
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
    config_types = ['100p', 'non_flexible']
    year = '2050'
    
    scenarios = {}
    
    # 为每个flexibility-demand-market组合创建场景
    for flexibility in flexibility_levels:
        for demand in demand_levels:
            for market in market_levels:
                scenario_code = f"{flexibility}{demand}{market}"
                scenarios[scenario_code] = {}
                
                # 为每个config_type创建版本信息
                for config_type in config_types:
                    # 构建版本名称
                    version_name = f"{base_version}-{scenario_code}-{year}-{config_type}"
                    version_dir = results_path / f"version-{version_name}"
                    
                    scenarios[scenario_code][config_type] = {
                        'version_name': version_name,
                        'version_dir': version_dir,
                        'year': year,
                        'flexibility': flexibility,
                        'demand': demand,
                        'market': market,
                        'config_type': config_type
                    }
    
    logger.info(f"基于基准版本 {base_version} 构建了 {len(scenarios)} 个场景")
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
                            logger.debug(f"成功加载 {config_type} 的 {file_type} 数据")
                        else:
                            logger.warning(f"无法加载 {config_type} 的 {file_type} 数据，将使用空数据")
                            # 创建空的DataFrame，确保数据结构一致
                            data[config_type] = pd.DataFrame()
                    except Exception as e:
                        logger.warning(f"加载 {config_type} 的 {file_type} 数据时出错: {str(e)}，将使用空数据")
                        # 创建空的DataFrame，确保数据结构一致
                        data[config_type] = pd.DataFrame()
                else:
                    logger.warning(f"文件不存在: {file_path}，将使用空数据")
                    # 创建空的DataFrame，确保数据结构一致
                    data[config_type] = pd.DataFrame()
            else:
                logger.warning(f"年份目录不存在: {year_dir}，将使用空数据")
                # 创建空的DataFrame，确保数据结构一致
                data[config_type] = pd.DataFrame()
        else:
            logger.warning(f"Summary目录不存在: {summary_dir}，将使用空数据")
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
        
        # capital–demand side - 需求侧资本成本
        ('capital', 'heat pump'): 'capital–demand side',
        ('capital', 'resistive heater'): 'capital–demand side',
        
        # capital–renewable - 可再生能源资本成本
        ('capital', 'hydro_inflow'): 'capital–renewable',
        ('capital', 'hydroelectricity'): 'capital–renewable',
        ('capital', 'offwind'): 'capital–renewable',
        ('capital', 'onwind'): 'capital–renewable',
        ('capital', 'solar'): 'capital–renewable',
        ('capital', 'solar thermal'): 'capital–renewable',
        ('capital', 'biomass'): 'capital–renewable',
        ('capital', 'biogas'): 'capital–renewable',
        
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
        ('capital', 'CO2 capture'): 'carbon capture',
        ('marginal', 'CO2 capture'): 'carbon capture',
        ('capital', 'Sabatier'): 'synthetic fuels',
        ('marginal', 'Sabatier'): 'synthetic fuels',
        ('capital', 'CO2'): 'carbon management',
        ('marginal', 'CO2'): 'carbon management',
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
        'synthetic fuels',
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
    for category, change in category_changes.items():
        if category not in exclude_categories and abs(change) > 1e-6:
            filtered_categories[category] = change
    
    # 添加总变化量
    filtered_categories['Total Change'] = total_change
    
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
    
    # 定义场景代码到描述的映射
    scenario_descriptions = {
        'L': 'Low',
        'M': 'Mid', 
        'H': 'High',
        'N': 'Non-constrained'
    }
    
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    # 为每个demand-market组合创建子图
    fig, axes = plt.subplots(3, 3, figsize=(20, 16))
    fig.suptitle(f'Scenario Comparison: {file_type.capitalize()} Changes by Flexibility Level', fontsize=16, y=0.98)
    
    # 定义demand和market级别
    demand_levels = ['L', 'M', 'H']
    market_levels = ['L', 'M', 'H']
    
    # 为每个demand-market组合创建子图
    for i, demand in enumerate(demand_levels):
        for j, market in enumerate(market_levels):
            ax = axes[i, j]
            
            # 收集该demand-market组合下所有flexibility级别的数据
            flexibility_levels = ['L', 'M', 'H', 'N']
            all_flex_data = {}
            
            for flex in flexibility_levels:
                scenario_code = f"{flex}{demand}{market}"
                if scenario_code in scenarios:
                    scenario_data = load_scenario_data(scenarios[scenario_code], file_type)
                    
                    # 检查是否有数据
                    has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
                    has_non_flex_data = 'non_flexible' in scenario_data and not scenario_data['non_flexible'].empty
                    
                    if has_100p_data and has_non_flex_data:
                        # 计算成本差异
                        cost_diff = calculate_cost_difference(scenario_data['100p'], scenario_data['non_flexible'])
                        if cost_diff:
                            # 转换为人民币并排除aluminum相关数据
                            cost_diff_cny = {}
                            for k, v in cost_diff.items():
                                if 'aluminum' not in k.lower():  # 排除aluminum相关数据
                                    cost_diff_cny[k] = v * EUR_TO_CNY
                            
                            if cost_diff_cny:  # 如果有非aluminum数据
                                all_flex_data[flex] = cost_diff_cny
            
            if all_flex_data:
                # 准备堆叠图数据
                all_categories = set()
                for flex_data in all_flex_data.values():
                    all_categories.update(flex_data.keys())
                
                if all_categories:
                    # 按flexibility级别组织数据
                    flex_names = list(all_flex_data.keys())
                    categories = list(all_categories)
                    
                    # 创建堆叠条形图数据
                    positive_data = []
                    negative_data = []
                    
                    for flex in flex_names:
                        if flex in all_flex_data:
                            flex_positive = []
                            flex_negative = []
                            
                            for category in categories:
                                value = all_flex_data[flex].get(category, 0)
                                if value > 0:
                                    flex_positive.append(value)
                                    flex_negative.append(0)
                                else:
                                    flex_positive.append(0)
                                    flex_negative.append(abs(value))
                            
                            positive_data.append(flex_positive)
                            negative_data.append(flex_negative)
                    
                    # 创建堆叠条形图
                    x_pos = np.arange(len(categories))
                    width = 0.6
                    
                    # 绘制正值的堆叠图
                    bottom_positive = np.zeros(len(categories))
                    for i_flex, flex in enumerate(flex_names):
                        if flex in all_flex_data:
                            ax.bar(x_pos, positive_data[i_flex], width, bottom=bottom_positive, 
                                   label=f'{flex} (Flexibility)', alpha=0.7)
                            bottom_positive += positive_data[i_flex]
                    
                    # 绘制负值的堆叠图
                    bottom_negative = np.zeros(len(categories))
                    for i_flex, flex in enumerate(flex_names):
                        if flex in all_flex_data:
                            ax.bar(x_pos, [-v for v in negative_data[i_flex]], width, bottom=bottom_negative, 
                                   alpha=0.7)
                            bottom_negative += negative_data[i_flex]
                    
                    # 设置标签
                    ax.set_xticks(x_pos)
                    ax.set_xticklabels(categories, fontsize=8, rotation=45, ha='right')
                    ax.set_ylabel('Cost Change (CNY)', fontsize=9)
                    
                    # 添加零线
                    ax.axhline(y=0, color='black', linestyle='-', alpha=0.5)
                    
                    # 设置标题
                    ax.set_title(f'Demand: {scenario_descriptions[demand]}, Market: {scenario_descriptions[market]}', 
                               fontsize=10, fontweight='bold')
                    
                    # 添加图例
                    ax.legend(fontsize=8, loc='upper right')
                    
                    # 添加网格
                    ax.grid(True, alpha=0.3, axis='y')
                else:
                    ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', 
                           transform=ax.transAxes, fontsize=10)
            else:
                ax.text(0.5, 0.5, f'No data for\nDemand:{demand}, Market:{market}', ha='center', va='center', 
                       transform=ax.transAxes, fontsize=10)
    
    plt.tight_layout()
    
    # 保存图表
    plot_file = plots_dir / f"scenario_comparison_{file_type}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Scenario comparison plot saved to: {plot_file}")
    
    plt.close()

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
        if len(scenario_code) == 2:  # 确保是有效的场景代码
            demand, market = scenario_code[0], scenario_code[1]
            
            # 加载场景数据
            scenario_data = load_scenario_data(scenario_info, file_type)
            
            # 检查是否有数据
            has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
            has_non_flex_data = 'non_flexible' in scenario_data and not scenario_data['non_flexible'].empty
            
            if has_100p_data and has_non_flex_data:
                # 计算成本差异
                cost_diff = calculate_cost_difference(scenario_data['100p'], scenario_data['non_flexible'])
                
                if cost_diff and 'Total Change' in cost_diff:
                    total_change = cost_diff['Total Change'] * EUR_TO_CNY
                    
                    # 计算各分类的变化量
                    category_changes = {}
                    for category, change in cost_diff.items():
                        if category != 'Total Change':
                            category_changes[category] = change * EUR_TO_CNY
                    
                    # 找出变化最大的前3个分类
                    top_changes = sorted(category_changes.items(), key=lambda x: abs(x[1]), reverse=True)[:3]
                    
                    row = {
                        'Scenario': scenario_code,
                        'Demand': demand,
                        'Market': market,
                        'Total Change (Billion CNY)': total_change / 1e9,
                        'Top Category 1': f"{top_changes[0][0]}: {top_changes[0][1]/1e9:.2f}B" if len(top_changes) > 0 else "N/A",
                        'Top Category 2': f"{top_changes[1][0]}: {top_changes[1][1]/1e9:.2f}B" if len(top_changes) > 1 else "N/A",
                        'Top Category 3': f"{top_changes[2][0]}: {top_changes[2][1]/1e9:.2f}B" if len(top_changes) > 2 else "N/A"
                    }
                    table_data.append(row)
                else:
                    row = {
                        'Scenario': scenario_code,
                        'Demand': demand,
                        'Market': market,
                        'Total Change (Billion CNY)': 'N/A',
                        'Top Category 1': 'N/A',
                        'Top Category 2': 'N/A',
                        'Top Category 3': 'N/A'
                    }
                    table_data.append(row)
            else:
                # 确定缺失的数据类型
                if not has_100p_data and not has_non_flex_data:
                    missing_info = "Both missing"
                elif not has_100p_data:
                    missing_info = "100p missing"
                else:
                    missing_info = "Non_flexible missing"
                
                row = {
                    'Scenario': scenario_code,
                    'Demand': demand,
                    'Market': market,
                    'Total Change (Billion CNY)': f'Data missing ({missing_info})',
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
        logger.info(f"Summary table saved to: {table_file}")
        
        # 打印表格
        print(f"\n=== {file_type.capitalize()} Summary Table ===")
        print(df.to_string(index=False))
    else:
        logger.warning("No data available for summary table")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='可视化generate_value_test_configs生成的场景结果')
    parser.add_argument('--results-dir', default='results', help='结果目录路径 (默认: results)')
    parser.add_argument('--output', default='results/scenario_analysis', help='输出目录')
    parser.add_argument('--file-type', choices=['costs', 'capacities'], default='costs', 
                       help='分析的文件类型 (默认: costs)')
    parser.add_argument('--verbose', '-v', action='store_true', help='详细输出')
    parser.add_argument('--config', default='config.yaml', help='配置文件路径 (默认: config.yaml)')
    
    args = parser.parse_args()
    
    # 设置日志级别
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    
    logger.info(f"开始分析场景结果，文件类型: {args.file_type}")
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
    
    # 生成场景对比图表
    logger.info("生成场景对比图表...")
    generate_scenario_plots(scenarios, output_path, args.file_type)
    
    # 生成摘要表格
    logger.info("生成摘要表格...")
    generate_summary_table(scenarios, output_path, args.file_type)
    
    # 打印数据完整性统计
    print_data_completeness_stats(scenarios, args.file_type)
    
    logger.info("分析完成！")
    logger.info(f"结果保存在: {output_path}")

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
    print(f"\n=== 数据完整性统计 ({file_type}) ===")
    
    total_scenarios = len(scenarios)
    complete_scenarios = 0
    missing_100p = 0
    missing_non_flex = 0
    missing_both = 0
    
    for scenario_code, scenario_info in scenarios.items():
        scenario_data = load_scenario_data(scenario_info, file_type)
        
        has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
        has_non_flex_data = 'non_flexible' in scenario_info and not scenario_data['non_flexible'].empty
        
        if has_100p_data and has_non_flex_data:
            complete_scenarios += 1
        elif not has_100p_data and not has_non_flex_data:
            missing_both += 1
        elif not has_100p_data:
            missing_100p += 1
        else:
            missing_non_flex += 1
    
    print(f"总场景数: {total_scenarios}")
    print(f"数据完整: {complete_scenarios} ({complete_scenarios/total_scenarios*100:.1f}%)")
    print(f"100p数据缺失: {missing_100p} ({missing_100p/total_scenarios*100:.1f}%)")
    print(f"non_flexible数据缺失: {missing_non_flex} ({missing_non_flex/total_scenarios*100:.1f}%)")
    print(f"两者都缺失: {missing_both} ({missing_both/total_scenarios*100:.1f}%)")
    
    # 按flexibility-demand-market组合显示详细统计
    print(f"\n=== 按场景组合的数据完整性 ===")
    flexibility_levels = ['L', 'M', 'H', 'N']
    demand_levels = ['L', 'M', 'H']
    market_levels = ['L', 'M', 'H']
    
    for flexibility in flexibility_levels:
        for demand in demand_levels:
            for market in market_levels:
                scenario_code = f"{flexibility}{demand}{market}"
                if scenario_code in scenarios:
                    scenario_data = load_scenario_data(scenarios[scenario_code], file_type)
                    
                    has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
                    has_non_flex_data = 'non_flexible' in scenario_data and not scenario_data['non_flexible'].empty
                    
                    if has_100p_data and has_non_flex_data:
                        status = "✓ Complete"
                    elif not has_100p_data and not has_non_flex_data:
                        status = "✗ Both missing"
                    elif not has_100p_data:
                        status = "✗ 100p missing"
                    else:
                        status = "✗ Non_flex missing"
                    
                    print(f"  {scenario_code} (F:{flexibility}, D:{demand}, M:{market}): {status}")

if __name__ == "__main__":
    main()
