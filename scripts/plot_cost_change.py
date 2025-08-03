# SPDX-FileCopyrightText: : 2022 The PyPSA-China Authors
#
# SPDX-License-Identifier: MIT

"""
对比两个不同版本的结果summary的cost和capacity差异
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

logger = logging.getLogger(__name__)

# 在这里设置要对比的版本号
VERSION1 = "0701.1H.9"  # 基准版本
VERSION2 = "0701.1H.10"  # 对比版本

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def find_summary_files(summary_dir):
    """
    递归查找summary目录下的所有CSV文件
    
    Parameters:
    -----------
    summary_dir : str
        summary目录路径
        
    Returns:
    --------
    list
        找到的CSV文件路径列表
    """
    summary_path = Path(summary_dir)
    csv_files = []
    
    if summary_path.exists():
        # 递归查找所有CSV文件
        csv_files = list(summary_path.rglob("*.csv"))
        logger.info(f"在{summary_dir}中找到{len(csv_files)}个CSV文件")
    else:
        logger.warning(f"Summary目录不存在: {summary_dir}")
    
    return csv_files

def load_summary_data(summary_dir):
    """
    加载summary目录下的所有CSV文件
    
    Parameters:
    -----------
    summary_dir : str
        summary目录路径
        
    Returns:
    --------
    dict
        包含所有summary数据的字典
    """
    summary_data = {}
    
    # 查找所有CSV文件
    csv_files = find_summary_files(summary_dir)
    
    if not csv_files:
        logger.error(f"在{summary_dir}中没有找到CSV文件")
        return summary_data
    
    for csv_file in csv_files:
        try:
            # 首先尝试读取文件的前几行来判断结构
            with open(csv_file, 'r') as f:
                first_lines = [f.readline().strip() for _ in range(5)]
            
            # 检查是否有多级索引（通过检查前几行是否有空列）
            has_multiindex = any(',' in line and line.split(',')[1] == '' for line in first_lines)
            
            if has_multiindex:
                # 多级索引文件 - 读取时保留所有列
                df = pd.read_csv(csv_file, header=None)
                
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
                    df = pd.read_csv(csv_file, index_col=[0, 1])
                    numeric_col = df.columns[0]
                    df[numeric_col] = pd.to_numeric(df[numeric_col], errors='coerce')
            else:
                # 普通单级索引文件
                df = pd.read_csv(csv_file, index_col=0)
                # 尝试将所有列转换为数值类型
                for col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            # 使用文件名作为键，避免重复
            key = csv_file.stem
            if key in summary_data:
                logger.warning(f"发现重复的文件名: {csv_file.name}，跳过")
                continue
                
            summary_data[key] = df
            logger.info(f"成功加载: {csv_file.name}")
            
        except Exception as e:
            logger.warning(f"加载{csv_file.name}时出错: {str(e)}")
            continue
    
    return summary_data

def compare_dataframes(df1, df2, name1, name2, comparison_type):
    """
    对比两个DataFrame的差异
    
    Parameters:
    -----------
    df1 : pd.DataFrame
        第一个DataFrame
    df2 : pd.DataFrame
        第二个DataFrame
    name1 : str
        第一个版本名称
    name2 : str
        第二个版本名称
    comparison_type : str
        比较类型（cost或capacity）
        
    Returns:
    --------
    pd.DataFrame
        包含差异信息的DataFrame
    """
    # 确保两个DataFrame有相同的索引和列
    common_index = df1.index.intersection(df2.index)
    common_columns = df1.columns.intersection(df2.columns)
    
    if len(common_index) == 0:
        logger.warning(f"{comparison_type}: 没有共同的索引")
        return pd.DataFrame()
    
    if len(common_columns) == 0:
        logger.warning(f"{comparison_type}: 没有共同的列")
        return pd.DataFrame()
    
    # 选择共同的数据
    df1_common = df1.loc[common_index, common_columns]
    df2_common = df2.loc[common_index, common_columns]
    
    # 确保数据类型为数值型
    for col in common_columns:
        df1_common[col] = pd.to_numeric(df1_common[col], errors='coerce')
        df2_common[col] = pd.to_numeric(df2_common[col], errors='coerce')
    
    # 计算差异
    absolute_diff = df2_common - df1_common
    relative_diff = (absolute_diff / df1_common) * 100  # 百分比差异
    
    # 创建结果DataFrame，保留原始索引
    results = pd.DataFrame(index=common_index)
    
    for col in common_columns:
        # 绝对差异
        abs_col = f"{col}_absolute_diff"
        results[abs_col] = absolute_diff[col]
        
        # 相对差异（百分比）
        rel_col = f"{col}_relative_diff_percent"
        results[rel_col] = relative_diff[col]
        
        # 原始值
        orig1_col = f"{col}_{name1}"
        orig2_col = f"{col}_{name2}"
        results[orig1_col] = df1_common[col]
        results[orig2_col] = df2_common[col]
    
    return results

def generate_comparison_report(data1, data2, name1, name2, output_dir):
    """
    生成对比报告
    
    Parameters:
    -----------
    data1 : dict
        第一个版本的数据
    data2 : dict
        第二个版本的数据
    name1 : str
        第一个版本名称
    name2 : str
        第二个版本名称
    output_dir : str
        输出目录
    """
    # 创建输出目录
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 对比costs和capacities
    comparison_files = ['costs', 'capacities']
    
    for file_type in comparison_files:
        if file_type in data1 and file_type in data2:
            logger.info(f"正在对比 {file_type}...")
            
            comparison_result = compare_dataframes(
                data1[file_type], 
                data2[file_type], 
                name1, 
                name2, 
                file_type
            )
            
            if not comparison_result.empty:
                # 保存对比结果
                output_file = output_path / f"{file_type}_comparison_{name1}_vs_{name2}.csv"
                comparison_result.to_csv(output_file)
                logger.info(f"对比结果已保存到: {output_file}")
                
                # 生成统计摘要
                generate_summary_statistics(comparison_result, file_type, name1, name2, output_path)
            else:
                logger.warning(f"{file_type} 对比结果为空")
        else:
            logger.warning(f"缺少 {file_type} 文件，跳过对比")

def generate_summary_statistics(comparison_df, file_type, name1, name2, output_path):
    """
    生成统计摘要
    
    Parameters:
    -----------
    comparison_df : pd.DataFrame
        对比结果DataFrame
    file_type : str
        文件类型（costs或capacities）
    name1 : str
        第一个版本名称
    name2 : str
        第二个版本名称
    output_path : Path
        输出路径
    """
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    # 找出绝对差异列
    abs_diff_cols = [col for col in comparison_df.columns if col.endswith('_absolute_diff')]
    
    if not abs_diff_cols:
        logger.warning(f"没有找到绝对差异列用于 {file_type}")
        return
    
    # 创建统计摘要
    summary_stats = []
    
    for col in abs_diff_cols:
        original_col = col.replace('_absolute_diff', '')
        
        # 获取对应的原始值列
        orig1_col = f"{original_col}_{name1}"
        orig2_col = f"{original_col}_{name2}"
        
        if orig1_col in comparison_df.columns and orig2_col in comparison_df.columns:
            # 计算统计信息
            abs_diff = comparison_df[col]
            rel_diff_col = col.replace('_absolute_diff', '_relative_diff_percent')
            rel_diff = comparison_df[rel_diff_col] if rel_diff_col in comparison_df.columns else pd.Series()
            
            # 转换为人民币
            orig1_sum_cny = comparison_df[orig1_col].sum() * EUR_TO_CNY
            orig2_sum_cny = comparison_df[orig2_col].sum() * EUR_TO_CNY
            abs_diff_sum_cny = abs_diff.sum() * EUR_TO_CNY
            abs_diff_max_cny = abs_diff.max() * EUR_TO_CNY
            abs_diff_min_cny = abs_diff.min() * EUR_TO_CNY
            abs_diff_std_cny = abs_diff.std() * EUR_TO_CNY
            
            stats = {
                '指标': original_col,
                f'{name1}_总值(CNY)': orig1_sum_cny,
                f'{name2}_总值(CNY)': orig2_sum_cny,
                '绝对差异总和(CNY)': abs_diff_sum_cny,
                '相对差异平均值(%)': rel_diff.mean() if not rel_diff.empty else np.nan,
                '最大绝对差异(CNY)': abs_diff_max_cny,
                '最小绝对差异(CNY)': abs_diff_min_cny,
                '标准差(CNY)': abs_diff_std_cny,
                '变化幅度(%)': ((orig2_sum_cny - orig1_sum_cny) / orig1_sum_cny * 100) if orig1_sum_cny != 0 else np.nan
            }
            
            summary_stats.append(stats)
    
    if summary_stats:
        summary_df = pd.DataFrame(summary_stats)
        summary_file = output_path / f"{file_type}_summary_{name1}_vs_{name2}.csv"
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"统计摘要已保存到: {summary_file}")
        
        # 打印主要变化
        print(f"\n=== {file_type.upper()} 主要变化 ({name1} vs {name2}) ===")
        for _, row in summary_df.iterrows():
            if abs(row['绝对差异总和(CNY)']) > 0:
                change_direction = "增加" if row['绝对差异总和(CNY)'] > 0 else "减少"
                print(f"{row['指标']}: {change_direction} {abs(row['绝对差异总和(CNY)'])/1e9:.3f}B CNY "
                      f"({row['变化幅度(%)']:.2f}%)")

def generate_yearly_comparison_plots(data1, data2, name1, name2, output_dir):
    """
    生成按年份的对比图表
    """
    # 创建输出目录
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # 从目录结构中获取所有年份数据
    years_data = collect_all_years_data(name1, name2, 'costs')
    
    if years_data:
        # 只生成成本变化量图表，不再生成total cost by carrier图表
        generate_cost_change_plot(years_data, name1, name2, plots_dir)
    else:
        logger.warning("没有找到年份信息用于成本分析")

def collect_all_years_data(name1, name2, file_type):
    """
    从所有年份目录中收集数据
    
    Parameters:
    -----------
    name1 : str
        第一个版本名称
    name2 : str
        第二个版本名称
    file_type : str
        文件类型
        
    Returns:
    --------
    dict
        按年份组织的数据
    """
    years_data = {}
    
    # 定义年份列表
    years = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
    
    for year in years:
        # 构建目录路径
        dir_pattern = f"postnetwork-ll-current+Neighbor-linear2050-{year}"
        
        # 检查两个版本的目录是否存在
        version1_dir = Path(f"results/version-{name1}/summary/postnetworks/positive/{dir_pattern}")
        version2_dir = Path(f"results/version-{name2}/summary/postnetworks/positive/{dir_pattern}")
        
        if version1_dir.exists() and version2_dir.exists():
            # 读取costs文件
            costs_file1 = version1_dir / f"{file_type}.csv"
            costs_file2 = version2_dir / f"{file_type}.csv"
            
            if costs_file1.exists() and costs_file2.exists():
                try:
                    # 读取数据
                    df1 = load_single_csv_file(costs_file1)
                    df2 = load_single_csv_file(costs_file2)
                    
                    if df1 is not None and df2 is not None:
                        years_data[year] = {
                            name1: df1,
                            name2: df2
                        }
                        logger.info(f"成功加载 {year} 年的 {file_type} 数据")
                    else:
                        logger.warning(f"无法加载 {year} 年的 {file_type} 数据")
                        
                except Exception as e:
                    logger.warning(f"加载 {year} 年的 {file_type} 数据时出错: {str(e)}")
            else:
                logger.warning(f"在 {year} 年目录中未找到 {file_type}.csv 文件")
        else:
            logger.warning(f"年份 {year} 的目录不存在")
    
    return years_data

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

def generate_cost_change_plot(years_data, name1, name2, plots_dir):
    """
    生成成本变化量图表
    
    Parameters:
    -----------
    years_data : dict
        按年份组织的数据
    name1 : str
        第一个版本名称
    name2 : str
        第二个版本名称
    plots_dir : Path
        图表输出目录
    """
    # 提取年份
    years = sorted(years_data.keys())
    
    # 定义资源类型映射 - 扩展映射以处理更多carrier类型
    carrier_mapping = {
        'coal': 'Coal',
        'coal power plant': 'Coal',  # 统一coal相关
        'coal cc': 'Coal CC',  # 碳捕获煤电
        'gas': 'Natural Gas',
        'gas boiler': 'Gas Boiler',  # 燃气锅炉
        'hydro_inflow': 'Hydro',
        'hydroelectricity': 'Hydro',
        'nuclear': 'Nuclear',
        'offwind': 'Offshore Wind',
        'onwind': 'Onshore Wind',
        'solar': 'Solar PV',
        'solar thermal': 'Solar Thermal',
        'biomass': 'Biomass',
        'battery': 'Battery',
        'battery discharger': 'Battery',  # 统一battery相关
        'PHS': 'Pumped Hydro',
        'water tanks': 'Heat Storage',
        'heat pump': 'Heat Pump',
        'resistive heater': 'Electric Heater',
        'CHP coal': 'Coal CHP',
        'CHP gas': 'Gas CHP',
        'OCGT gas': 'Gas OCGT',
        'coal boiler': 'Coal Boiler',
        'AC': 'AC Transmission',
        'stations': 'Substations',
        'biogas': 'Biogas',
        'CO2 capture': 'CO2 Capture',  # 碳捕获
        'H2': 'Hydrogen',  # 氢气
        'H2 CHP': 'H2 CHP',  # 氢气热电联产
        'Sabatier': 'Sabatier',  # Sabatier反应
        'CO2': 'CO2'  # 二氧化碳
    }
    
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
        
        # 其他分类（如果需要）
        ('capital', 'CO2 capture'): 'carbon capture',
        ('marginal', 'CO2 capture'): 'carbon capture',
        ('capital', 'Sabatier'): 'synthetic fuels',
        ('marginal', 'Sabatier'): 'synthetic fuels',
        ('capital', 'CO2'): 'carbon management',
        ('marginal', 'CO2'): 'carbon management',
    }
    
    # 按成本分类组织数据
    category_changes = {}
    net_changes = []
    
    for year in years:
        if name1 in years_data[year] and name2 in years_data[year]:
            df1 = years_data[year][name1]
            df2 = years_data[year][name2]
            
            year_net_change = 0
            
            # 计算每个成本分类的变化量
            for idx in df1.index:
                if len(idx) >= 3:
                    component_type, cost_type, carrier = idx[0], idx[1], idx[2]
                    
                    # 使用新的分类映射
                    category_key = (cost_type, carrier)
                    category_name = cost_category_mapping.get(category_key, f"{cost_type} - {carrier}")
                    
                    if category_name not in category_changes:
                        category_changes[category_name] = {'years': [], 'changes': []}
                    
                    # 获取两个版本的对应值
                    v1_value = df1.loc[idx].iloc[0]
                    v2_value = 0
                    
                    # 在第二个版本中查找对应值
                    for idx2 in df2.index:
                        if len(idx2) >= 3 and idx2[0] == component_type and idx2[1] == cost_type and idx2[2] == carrier:
                            v2_value = df2.loc[idx2].iloc[0]
                            break
                    
                    # 处理NaN值
                    if pd.isna(v1_value):
                        v1_value = 0
                    if pd.isna(v2_value):
                        v2_value = 0
                    
                    # 计算变化量（v1 - v2，节约为正，增加为负）
                    change = v1_value - v2_value
                    
                    if year not in category_changes[category_name]['years']:
                        category_changes[category_name]['years'].append(year)
                        category_changes[category_name]['changes'].append(0)
                    
                    year_idx = category_changes[category_name]['years'].index(year)
                    category_changes[category_name]['changes'][year_idx] += change
                    year_net_change += change
            
            net_changes.append(year_net_change)
    
    # 打印成本分类调试信息
    print(f"\n=== 成本分类调试信息 ===")
    print("原始成本分类名称:")
    for category in category_changes.keys():
        print(f"  - {category}")
    
    # 过滤掉nan分类和变化量为0的分类，以及不需要展示的分类
    filtered_categories = {}
    exclude_categories = {
        'nan - nan',
        'marginal - renewable',
        'marginal - heat pump', 
        'marginal - resistive heater'
    }
    
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    for category, data in category_changes.items():
        if category == 'nan' or pd.isna(category) or category in exclude_categories:
            continue
        if any(abs(change) > 1e-6 for change in data['changes']):  # 使用小的阈值避免浮点数精度问题
            # 将欧元转换为人民币
            data['changes'] = [change * EUR_TO_CNY for change in data['changes']]
            filtered_categories[category] = data
    
    # 检查synthetic fuels, marginal - renewable, marginal - heat pump, marginal - resistive heater的数量，如果很少就不展示
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
    for category in exclude_categories:
        if category in filtered_categories:
            del filtered_categories[category]
            print(f"过滤掉 {category}")
    
    print(f"\n过滤后的分类数量: {len(filtered_categories)}")
    print("过滤后的分类名称:")
    for category in filtered_categories.keys():
        print(f"  - {category}")
    
    # 添加调试信息：显示每个分类的变化量总和（人民币）
    print(f"\n=== 各分类变化量调试信息（人民币）===")
    for category, data in filtered_categories.items():
        total_change = sum(data['changes'])
        max_change = max(data['changes']) if data['changes'] else 0
        min_change = min(data['changes']) if data['changes'] else 0
        print(f"{category}: 总和={total_change/1e9:.3f}B CNY, 最大值={max_change/1e9:.3f}B CNY, 最小值={min_change/1e9:.3f}B CNY")
    
    if not filtered_categories:
        logger.warning("没有找到有效的成本变化数据")
        return
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(14, 8))
    
    x = np.arange(len(years))
    
    # 准备堆叠数据
    positive_changes = []
    negative_changes = []
    positive_labels = []
    negative_labels = []
    
    # 分离正负变化量
    for category, data in filtered_categories.items():
        changes = data['changes']
        years_for_category = data['years']
        
        # 检查数据长度是否匹配
        if len(changes) == len(x):
            positive_changes.append([max(0, change) for change in changes])
            negative_changes.append([min(0, change) for change in changes])
            positive_labels.append(category)
            negative_labels.append(category)
        else:
            # 处理数据长度不匹配的情况
            logger.warning(f"分类 {category} 数据长度不匹配: {len(changes)} vs {len(x)}")
            logger.info(f"分类 {category} 有数据的年份: {years_for_category}")
            
            # 创建完整长度的数组，缺失年份用0填充
            full_changes = []
            for year in years:
                if year in years_for_category:
                    year_idx = years_for_category.index(year)
                    full_changes.append(changes[year_idx])
                else:
                    full_changes.append(0)  # 缺失年份用0填充
            
            positive_changes.append([max(0, change) for change in full_changes])
            negative_changes.append([min(0, change) for change in full_changes])
            positive_labels.append(category)
            negative_labels.append(category)
    
    # 绘制堆叠柱状图 - 使用更丰富的颜色映射避免重复
    # 使用tab20颜色映射，提供更多不同的颜色
    colors = plt.cm.tab20(np.linspace(0, 1, len(filtered_categories)))
    
    # 为每个carrier分配唯一的颜色
    carrier_colors = {}
    for i, carrier in enumerate(filtered_categories.keys()):
        carrier_colors[carrier] = colors[i]
    
    # 用于跟踪已经添加到legend的carrier
    added_to_legend = set()
    
    # 绘制正值堆叠
    bottom_positive = np.zeros(len(x))
    for i, (changes, carrier) in enumerate(zip(positive_changes, positive_labels)):
        if any(change > 0 for change in changes):
            # 只在第一次遇到该carrier时添加到legend
            label = carrier if carrier not in added_to_legend else ""
            bars = ax.bar(x, changes, bottom=bottom_positive, label=label, color=carrier_colors[carrier], alpha=0.8)
            bottom_positive += np.array(changes)
            added_to_legend.add(carrier)
    
    # 绘制负值堆叠
    bottom_negative = np.zeros(len(x))
    for i, (changes, carrier) in enumerate(zip(negative_changes, negative_labels)):
        if any(change < 0 for change in changes):
            # 负值部分也添加到legend，但只在第一次遇到该carrier时添加
            label = carrier if carrier not in added_to_legend else ""
            bars = ax.bar(x, changes, bottom=bottom_negative, label=label, color=carrier_colors[carrier], alpha=0.8)
            bottom_negative += np.array(changes)
            added_to_legend.add(carrier)
    
    # 打印净变化量数值用于调试
    print(f"\n=== 净变化量数值调试信息（人民币）===")
    for i, (year, net_change) in enumerate(zip(years, net_changes)):
        net_change_cny = net_change * EUR_TO_CNY  # 转换为人民币
        print(f"{year}年: {net_change_cny/1e9:.3f}B CNY (原始值: {net_change:.0f} EUR)")
    
    # 绘制净变化量的黑线（加粗并确保可见）
    net_changes_cny = [change * EUR_TO_CNY for change in net_changes]  # 转换为人民币
    ax.plot(x, net_changes_cny, 'k-', linewidth=2, label='Net Change', marker='o', markersize=10, zorder=20)
    
    # 添加净变化量的数值标签 - 降低阈值显示所有数值
    for i, net_change in enumerate(net_changes_cny):
        if abs(net_change) > 1e3:  # 大幅降低阈值，显示更多数值标签
            ax.annotate(f'{net_change/1e9:.1f}B',
                        xy=(i, net_change),
                        xytext=(0, 10 if net_change > 0 else -20),
                        textcoords="offset points",
                        ha='center', va='bottom' if net_change > 0 else 'top', 
                        fontsize=10, weight='bold', color='black',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    

    
    ax.set_xlabel('Year')
    ax.set_ylabel('Cost Change (Billion CNY)')
    ax.set_title(f'Cost Change by Category ({name1} → {name2})')
    ax.set_xticks(x)
    ax.set_xticklabels(years, rotation=45)
    
    # 设置y轴刻度标签为十亿人民币单位
    y_ticks = ax.get_yticks()
    y_tick_labels = [f'{tick/1e9:.1f}B' for tick in y_ticks]
    ax.set_yticklabels(y_tick_labels)
    
    # 设置y轴范围，最高标尺为120B CNY
    ax.set_ylim(bottom=None, top=120e9)
    
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    
    # 保存图表
    plot_file = plots_dir / f"cost_change_by_category_{name1}_vs_{name2}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Cost change plot saved to: {plot_file}")
    
    # 保存数据（人民币）
    data_rows = []
    for category, data in filtered_categories.items():
        for i, year in enumerate(data['years']):
            data_rows.append({
                'Cost Category': category,
                'Year': year,
                'Cost Change (CNY)': data['changes'][i],
                'Cost Change (Billion CNY)': data['changes'][i] / 1e9
            })
    
    # 添加净变化量数据（人民币）
    for i, year in enumerate(years):
        if i < len(net_changes):
            net_change_cny = net_changes[i] * EUR_TO_CNY
            data_rows.append({
                'Cost Category': 'Net Change',
                'Year': year,
                'Cost Change (CNY)': net_change_cny,
                'Cost Change (Billion CNY)': net_change_cny / 1e9
            })
    
    data_df = pd.DataFrame(data_rows)
    data_file = plots_dir / f"cost_change_data_{name1}_vs_{name2}.csv"
    data_df.to_csv(data_file, index=False)
    logger.info(f"Cost change data saved to: {data_file}")
    
    plt.close()

def collect_all_years_data(name1, name2, file_type):
    """
    从所有年份目录中收集数据
    
    Parameters:
    -----------
    name1 : str
        第一个版本名称
    name2 : str
        第二个版本名称
    file_type : str
        文件类型
        
    Returns:
    --------
    dict
        按年份组织的数据
    """
    years_data = {}
    
    # 定义年份列表
    years = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
    
    for year in years:
        # 构建目录路径
        dir_pattern = f"postnetwork-ll-current+Neighbor-linear2050-{year}"
        
        # 检查两个版本的目录是否存在
        version1_dir = Path(f"results/version-{name1}/summary/postnetworks/positive/{dir_pattern}")
        version2_dir = Path(f"results/version-{name2}/summary/postnetworks/positive/{dir_pattern}")
        
        if version1_dir.exists() and version2_dir.exists():
            # 读取costs文件
            costs_file1 = version1_dir / f"{file_type}.csv"
            costs_file2 = version2_dir / f"{file_type}.csv"
            
            if costs_file1.exists() and costs_file2.exists():
                try:
                    # 读取数据
                    df1 = load_single_csv_file(costs_file1)
                    df2 = load_single_csv_file(costs_file2)
                    
                    if df1 is not None and df2 is not None:
                        years_data[year] = {
                            name1: df1,
                            name2: df2
                        }
                        logger.info(f"成功加载 {year} 年的 {file_type} 数据")
                    else:
                        logger.warning(f"无法加载 {year} 年的 {file_type} 数据")
                        
                except Exception as e:
                    logger.warning(f"加载 {year} 年的 {file_type} 数据时出错: {str(e)}")
            else:
                logger.warning(f"在 {year} 年目录中未找到 {file_type}.csv 文件")
        else:
            logger.warning(f"年份 {year} 的目录不存在")
    
    return years_data

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

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='对比两个不同版本的结果summary')
    parser.add_argument('--name1', help='第一个版本显示名称 (默认使用VERSION1)')
    parser.add_argument('--name2', help='第二个版本显示名称 (默认使用VERSION2)')
    parser.add_argument('--output', default='results/comparison_results', help='输出目录')
    parser.add_argument('--verbose', '-v', action='store_true', help='详细输出')
    parser.add_argument('--results-dir', default='results', help='结果目录路径 (默认: results)')
    
    args = parser.parse_args()
    
    # 设置日志级别
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    
    # 使用脚本中设置的版本号
    version1 = VERSION1
    version2 = VERSION2
    
    # 构建summary目录路径 - 使用默认的summary目录结构
    summary_dir1 = f"{args.results_dir}/version-{version1}/summary"
    summary_dir2 = f"{args.results_dir}/version-{version2}/summary"
    
    # 设置版本名称
    name1 = args.name1 if args.name1 else version1
    name2 = args.name2 if args.name2 else version2
    
    logger.info(f"开始对比版本 {name1} 和 {name2} 的结果")
    logger.info(f"Summary目录1: {summary_dir1}")
    logger.info(f"Summary目录2: {summary_dir2}")
    
    # 检查目录是否存在
    summary_path1 = Path(summary_dir1)
    summary_path2 = Path(summary_dir2)
    
    if not summary_path1.exists():
        logger.error(f"第一个版本的summary目录不存在: {summary_dir1}")
        logger.error("请确保该版本已经运行完成并生成了summary文件")
        logger.info("期望的目录结构: results/version-{version}/summary/")
        return
    
    if not summary_path2.exists():
        logger.error(f"第二个版本的summary目录不存在: {summary_dir2}")
        logger.error("请确保该版本已经运行完成并生成了summary文件")
        logger.info("期望的目录结构: results/version-{version}/summary/")
        return
    
    # 检查是否有CSV文件
    csv_files1 = find_summary_files(summary_dir1)
    csv_files2 = find_summary_files(summary_dir2)
    
    if not csv_files1:
        logger.error(f"在第一个版本的summary目录中没有找到CSV文件: {summary_dir1}")
        return
    
    if not csv_files2:
        logger.error(f"在第二个版本的summary目录中没有找到CSV文件: {summary_dir2}")
        return
    
    # 加载数据
    logger.info(f"加载版本 {name1} 的数据...")
    data1 = load_summary_data(summary_dir1)
    
    logger.info(f"加载版本 {name2} 的数据...")
    data2 = load_summary_data(summary_dir2)
    
    if not data1:
        logger.error(f"无法加载版本 {name1} 的数据")
        return
    
    if not data2:
        logger.error(f"无法加载版本 {name2} 的数据")
        return
    
    # 生成对比报告
    generate_comparison_report(data1, data2, name1, name2, args.output)
    
    # 生成按年份的对比图表
    generate_yearly_comparison_plots(data1, data2, name1, name2, Path(args.output))
    
    logger.info("对比完成！")

if __name__ == "__main__":
    main() 