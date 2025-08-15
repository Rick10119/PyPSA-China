#!/usr/bin/env python3
"""
简化版场景对比可视化脚本
专门用于画指定情景（如HHH）的对比图
"""

import logging
import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import yaml

logger = logging.getLogger(__name__)

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def load_config(config_path):
    """加载配置文件"""
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        logger.error(f"加载配置文件 {config_path} 时出错: {str(e)}")
        return None

def get_base_version_from_config(config_path):
    """从配置文件中获取基准版本号"""
    config = load_config(config_path)
    if config is None:
        logger.warning(f"无法加载配置文件 {config_path}，使用默认版本号")
        return None
    
    version = config.get('version', '')
    return version

def find_single_scenario_results(results_dir, base_version, scenario_code):
    """
    查找指定场景的结果目录
    
    Parameters:
    -----------
    results_dir : str
        结果目录路径
    base_version : str
        基准版本号
    scenario_code : str
        场景代码，如 'HHH'
        
    Returns:
    --------
    dict
        场景信息
    """
    results_path = Path(results_dir)
    if not results_path.exists():
        logger.error(f"结果目录不存在: {results_dir}")
        return {}
    
    if len(scenario_code) != 3:
        logger.error(f"场景代码必须是3位字符，如 'HHH'，当前为: {scenario_code}")
        return {}
    
    flexibility, demand, market = scenario_code[0], scenario_code[1], scenario_code[2]
    year = '2050'
    
    # 创建100p版本信息
    version_name_100p = f"{base_version}-{scenario_code}-{year}-100p"
    version_dir_100p = results_path / f"version-{version_name_100p}"
    
    # 创建non_flexible版本信息（使用中等水平的flex和demand，但保持相同的market）
    baseline_flexibility = 'M'
    baseline_demand = 'M'
    baseline_market = market
    baseline_scenario_code = f"{baseline_flexibility}{baseline_demand}{baseline_market}"
    
    version_name_non_flex = f"{base_version}-{baseline_scenario_code}-{year}-non_flexible"
    version_dir_non_flex = results_path / f"version-{version_name_non_flex}"
    
    scenario_info = {
        '100p': {
            'version_name': version_name_100p,
            'version_dir': version_dir_100p,
            'year': year,
            'flexibility': flexibility,
            'demand': demand,
            'market': market,
            'config_type': '100p'
        },
        'non_flexible': {
            'version_name': version_name_non_flex,
            'version_dir': version_dir_non_flex,
            'year': year,
            'flexibility': baseline_flexibility,
            'demand': baseline_demand,
            'market': baseline_market,
            'config_type': 'non_flexible'
        }
    }
    
    logger.info(f"Scenario {scenario_code} information:")
    logger.info(f"  100p version: {version_name_100p}")
    logger.info(f"  100p directory: {version_dir_100p}")
    logger.info(f"  non_flexible version: {version_name_non_flex}")
    logger.info(f"  non_flexible directory: {version_dir_non_flex}")
    
    return scenario_info

def load_single_csv_file(file_path):
    """加载单个CSV文件"""
    try:
        # 首先尝试读取文件的前几行来判断结构
        with open(file_path, 'r') as f:
            first_lines = [f.readline().strip() for _ in range(5)]
        
        # 检查是否有多级索引
        has_multiindex = any(',' in line and line.split(',')[1] == '' for line in first_lines)
        
        if has_multiindex:
            # 多级索引文件
            df = pd.read_csv(file_path, header=None)
            
            if len(df.columns) >= 4:
                df.set_index([0, 1, 2], inplace=True)
                df.columns = [df.columns[0]]
                df[df.columns[0]] = pd.to_numeric(df[df.columns[0]], errors='coerce')
            else:
                df = pd.read_csv(file_path, index_col=[0, 1])
                numeric_col = df.columns[0]
                df[numeric_col] = pd.to_numeric(df[numeric_col], errors='coerce')
        else:
            # 普通单级索引文件
            df = pd.read_csv(file_path, index_col=0)
            for col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        return df
        
    except Exception as e:
        logger.warning(f"Error loading {file_path}: {str(e)}")
        return None

def load_scenario_data(scenario_info, file_type='costs'):
    """加载单个场景的数据"""
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
                            logger.info(f"Successfully loaded {config_type} {file_type} data, rows: {len(df)}")
                        else:
                            logger.warning(f"Cannot load {config_type} {file_type} data")
                            data[config_type] = pd.DataFrame()
                    except Exception as e:
                        logger.warning(f"Error loading {config_type} {file_type} data: {str(e)}")
                        data[config_type] = pd.DataFrame()
                else:
                    logger.warning(f"File does not exist: {file_path}")
                    data[config_type] = pd.DataFrame()
            else:
                logger.warning(f"Year directory does not exist: {year_dir}")
                data[config_type] = pd.DataFrame()
        else:
            logger.warning(f"Summary directory does not exist: {summary_dir}")
            data[config_type] = pd.DataFrame()
    
    return data

def calculate_cost_difference(costs_100p, costs_non_flex):
    """计算100p和non_flexible之间的成本差异"""
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
    filtered_total_change = 0
    
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

def plot_single_scenario_comparison(scenario_info, output_dir, file_type='costs', scenario_code='HHH'):
    """生成单个场景的对比图表"""
    # 创建输出目录
    plots_dir = output_dir / "single_scenario_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Plot output directory: {plots_dir}")
    
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    # 加载场景数据
    scenario_data = load_scenario_data(scenario_info, file_type)
    
    # 检查是否有数据
    has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
    has_non_flex_data = 'non_flexible' in scenario_data and not scenario_data['non_flexible'].empty
    
    if not has_100p_data or not has_non_flex_data:
        logger.error("Incomplete data, cannot generate plots")
        return
    
    # 计算成本差异
    cost_diff = calculate_cost_difference(scenario_data['100p'], scenario_data['non_flexible'])
    
    if not cost_diff:
        logger.warning("No valid cost difference data")
        return
    
    # 准备绘图数据
    plot_data = []
    for category, change in cost_diff.items():
        if category != 'Total Change':
            value_cny = change * EUR_TO_CNY
            plot_data.append({
                'Category': category,
                'Value (CNY)': value_cny,
                'Value (Billion CNY)': value_cny / 1e9
            })
    
    if not plot_data:
        logger.warning("No valid plot data")
        return
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(12, 8))
    fig.suptitle(f'Scenario {scenario_code} Analysis: {file_type.capitalize()}', fontsize=16, y=0.98)
    
    # 准备堆叠图数据
    df_plot = pd.DataFrame(plot_data)
    
    # 分离正负变化量
    positive_changes = []
    negative_changes = []
    positive_labels = []
    negative_labels = []
    
    # 过滤掉nan数据，按绝对值排序，取前10个最大的变化
    df_plot = df_plot.dropna()  # 删除包含nan的行
    # 额外过滤：删除Category中包含'nan'或'NaN'的行
    df_plot = df_plot[~df_plot['Category'].str.contains('nan', case=False, na=False)]
    df_plot['Abs_Value'] = df_plot['Value (CNY)'].abs()
    top_data = df_plot.nlargest(10, 'Abs_Value')
    
    for _, row in top_data.iterrows():
        value = row['Value (CNY)']
        category = row['Category']
        
        # 调整方向：成本减少（负值）显示在上方，成本增加（正值）显示在下方
        if value < 0:  # 成本减少，显示在上方
            positive_changes.append([abs(value)])  # 取绝对值
            positive_labels.append(category)
        else:  # 成本增加，显示在下方
            negative_changes.append([value])
            negative_labels.append(category)
    
    # 创建单个柱子
    x_pos = 0
    width = 0.8
    
    # 从配置文件读取成本分类颜色
    config = load_config('config.yaml')
    category_colors = config.get('cost_category_colors', {}) if config else {}
    
    # 为每个分类分配颜色，如果不在预定义中则使用默认颜色
    for i, category in enumerate(top_data['Category']):
        if category not in category_colors:
            # 使用默认颜色映射
            category_colors[category] = plt.cm.tab20(i % 20)
    
    # 用于跟踪已经添加到legend的分类
    added_to_legend = set()
    
    # 绘制正值堆叠（成本减少，在横轴上面）
    bottom_positive = 0
    for i, (changes, category) in enumerate(zip(positive_changes, positive_labels)):
        value = changes[0]
        # 只在第一次遇到该分类时添加到legend
        label = category if category not in added_to_legend else ""
        bar = ax.bar(x_pos, value, width, bottom=bottom_positive, 
                    label=label, color=category_colors[category], alpha=0.8)
        bottom_positive += value
        added_to_legend.add(category)
    
    # 绘制负值堆叠（成本增加，在横轴下面）
    bottom_negative = 0
    for i, (changes, category) in enumerate(zip(negative_changes, negative_labels)):
        value = changes[0]
        # 负值部分也添加到legend，但只在第一次遇到该分类时添加
        label = category if category not in added_to_legend else ""
        bar = ax.bar(x_pos, -value, width, bottom=bottom_negative,  # 使用负值
                    label=label, color=category_colors[category], alpha=0.8)
        bottom_negative += value
        added_to_legend.add(category)
    
    # 设置图表
    ax.set_xlabel('Scenario', fontsize=12)
    ax.set_ylabel('Cost Change (Billion CNY)', fontsize=12)
    ax.set_title('Cost Change Stacked Bar Chart', fontsize=14, fontweight='bold')
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1.5)
    ax.grid(True, alpha=0.3, axis='y')
    
    # 设置x轴
    ax.set_xticks([x_pos])
    ax.set_xticklabels([scenario_code], fontsize=12)
    
    # 设置y轴标签为十亿人民币单位
    y_ticks = ax.get_yticks()
    y_tick_labels = [f'{tick/1e9:.1f}B' for tick in y_ticks]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels, fontsize=10)
    
    # 添加图例
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    
    plt.tight_layout()
    
    plt.show()
    
    # 保存图表
    plot_file = plots_dir / f"scenario_{scenario_code}_{file_type}_comparison.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Scenario comparison plot saved to: {plot_file}")
    
    plt.close()
    
    # 保存数据到CSV
    csv_file = plots_dir / f"scenario_{scenario_code}_{file_type}_data.csv"
    df_plot.to_csv(csv_file, index=False)
    logger.info(f"Scenario data saved to: {csv_file}")

def generate_summary_table(scenario_info, output_dir, file_type='costs', scenario_code='HHH'):
    """生成场景摘要表格"""
    # 创建输出目录
    tables_dir = output_dir / "summary_tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    # 加载场景数据
    scenario_data = load_scenario_data(scenario_info, file_type)
    
    # 检查是否有数据
    has_100p_data = '100p' in scenario_data and not scenario_data['100p'].empty
    has_non_flex_data = 'non_flexible' in scenario_data and not scenario_data['non_flexible'].empty
    
    if not has_100p_data or not has_non_flex_data:
        logger.error("数据不完整，无法生成摘要表格")
        return
    
    # 计算成本差异
    cost_diff = calculate_cost_difference(scenario_data['100p'], scenario_data['non_flexible'])
    
    if not cost_diff or 'Total Change' not in cost_diff:
        logger.warning("没有有效的成本差异数据")
        return
    
    # 准备表格数据
    table_data = []
    total_change = cost_diff['Total Change'] * EUR_TO_CNY
    
    # 添加总变化行
    table_data.append({
        'Category': 'Total Change',
        'Change (Billion EUR)': cost_diff['Total Change'] / 1e9,
        'Change (Billion CNY)': total_change / 1e9,
        'Description': '总成本变化'
    })
    
    # 添加各分类变化
    for category, change in cost_diff.items():
        if category != 'Total Change':
            change_cny = change * EUR_TO_CNY
            table_data.append({
                'Category': category,
                'Change (Billion EUR)': change / 1e9,
                'Change (Billion CNY)': change_cny / 1e9,
                'Description': f'{category} 成本变化'
            })
    
    # 创建DataFrame并排序
    df = pd.DataFrame(table_data)
    df = df.sort_values('Change (Billion CNY)', key=abs, ascending=False)
    
    # 保存表格
    table_file = tables_dir / f"scenario_{scenario_code}_{file_type}_summary.csv"
    df.to_csv(table_file, index=False)
    logger.info(f"Summary table saved to: {table_file}")
    
    # 打印表格
    print(f"\n=== Scenario {scenario_code} {file_type.capitalize()} Summary Table ===")
    print(df.to_string(index=False))
    
    # 打印关键信息
    print(f"\nKey Information:")
    print(f"Scenario Code: {scenario_code}")
    print(f"Total Cost Change: {total_change/1e9:.2f}B CNY ({cost_diff['Total Change']/1e9:.2f}B EUR)")
    if total_change > 0:
        print(f"Result: 100p configuration saves {total_change/1e9:.2f}B CNY compared to non_flexible configuration")
    else:
        print(f"Result: 100p configuration increases {abs(total_change)/1e9:.2f}B CNY compared to non_flexible configuration")

def main():
    """主函数"""
    # ===== 在这里修改您要分析的情景 =====
    SCENARIO_CODE = "HHH"  # 修改这里：HHH, LLL, MMM, HHL, HHM, 等等
    FILE_TYPE = "costs"    # 修改这里：costs 或 capacities
    # =====================================
    
    # 其他参数（通常不需要修改）
    RESULTS_DIR = "results"
    OUTPUT_DIR = "results/single_scenario_analysis"
    CONFIG_FILE = "config.yaml"
    
    # 设置日志级别
    log_level = logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    
    logger.info(f"Starting analysis for scenario {SCENARIO_CODE}, file type: {FILE_TYPE}")
    logger.info(f"Results directory: {RESULTS_DIR}")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info(f"Config file: {CONFIG_FILE}")
    
    # 从配置文件读取基准版本号
    base_version = get_base_version_from_config(CONFIG_FILE)
    if base_version is None:
        logger.error("Cannot get base version, please check config file")
        return
    
    # 查找指定场景结果
    scenario_info = find_single_scenario_results(RESULTS_DIR, base_version, SCENARIO_CODE)
    
    if not scenario_info:
        logger.error(f"Scenario {SCENARIO_CODE} results not found")
        return
    
    # 创建输出目录
    output_path = Path(OUTPUT_DIR)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 生成场景对比图表
    logger.info("Generating scenario comparison plots...")
    plot_single_scenario_comparison(scenario_info, output_path, FILE_TYPE, SCENARIO_CODE)
    
    # 生成摘要表格
    logger.info("Generating summary table...")
    generate_summary_table(scenario_info, output_path, FILE_TYPE, SCENARIO_CODE)
    
    logger.info("Analysis completed!")
    logger.info(f"Results saved to: {output_path}")

if __name__ == "__main__":
    main()
