"""
对比两个场景的就业人数情况脚本
支持_20p和_non_flexible后缀的CSV文件输入，生成上下对比的图表
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import argparse
import yaml
from pathlib import Path

def set_plot_style():
    """
    设置绘图样式
    """
    # 设置英文字体
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'sans-serif']
    plt.rcParams['axes.unicode_minus'] = False
    
    plt.style.use(['classic', 'seaborn-v0_8-whitegrid',
                   {'axes.grid': False, 'grid.linestyle': '--', 'grid.color': u'0.6',
                    'hatch.color': 'white',
                    'patch.linewidth': 0.5,
                    'font.size': 16,
                    'legend.fontsize': 'large',
                    'lines.linewidth': 1.5,
                    'pdf.fonttype': 42,
                    }])

def load_csv_data(csv_file):
    """
    从CSV文件加载容量因子数据，全部使用平均因子
    
    Parameters:
    -----------
    csv_file : str
        CSV文件路径
    
    Returns:
    --------
    tuple
        (capacity_factors, load_factors) 两个DataFrame
    """
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"CSV file not found: {csv_file}")
    
    # 读取CSV文件
    df = pd.read_csv(csv_file, index_col='Month')
    
    # 分离平均容量因子和负荷因子数据
    avg_capacity_cols = [col for col in df.columns if 'Capacity_Factor_Avg' in col]
    load_cols = [col for col in df.columns if 'Load_Factor' in col]
    
    # 优先使用平均容量因子数据
    if avg_capacity_cols:
        # 使用平均容量因子数据
        capacity_factors = df[avg_capacity_cols]
        capacity_factors.columns = [col.replace('_Capacity_Factor_Avg', '') for col in capacity_factors.columns]
        print(f"使用月度平均容量因子数据计算就业: {os.path.basename(csv_file)}")
        
    else:
        # 兼容旧格式：只有Capacity_Factor列
        capacity_cols = [col for col in df.columns if 'Capacity_Factor' in col and 'Max' not in col and 'Avg' not in col]
        capacity_factors = df[capacity_cols] if capacity_cols else pd.DataFrame()
        capacity_factors.columns = [col.replace('_Capacity_Factor', '') for col in capacity_factors.columns]
        print(f"使用旧格式容量因子数据计算就业: {os.path.basename(csv_file)}")
    
    load_factors = df[load_cols] if load_cols else pd.DataFrame()
    # 重命名负荷因子列
    load_factors.columns = [col.replace('_Load_Factor', '') for col in load_factors.columns]
    
    return capacity_factors, load_factors

def load_employment_config(config_file=None):
    """
    加载就业参数配置文件
    
    Parameters:
    -----------
    config_file : str, optional
        配置文件路径，默认为 employment_config.yaml
    
    Returns:
    --------
    tuple
        (employment_params, config) 就业参数字典和配置字典
    """
    if config_file is None:
        config_file = os.path.join(os.path.dirname(__file__), 'employment_config.yaml')
    
    if not os.path.exists(config_file):
        print(f"Warning: Config file {config_file} not found, using default parameters")
        return get_default_employment_parameters(), {}
    
    try:
        with open(config_file, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        # 转换配置格式
        employment_params = {}
        for tech, params in config['industries'].items():
            # 处理不同的单位系统
            if 'employment_per_GW' in params:
                # 新单位系统：GW和千人/GW
                employment_params[tech] = {
                    'installed_capacity': params['installed_capacity'],  # 直接使用GW
                    'employment_per_gw': params['employment_per_GW'],     # 直接使用千人/GW
                    'display_name': params['display_name'],
                    'unit_system': 'GW'
                }
            else:
                # 旧单位系统：MW和人/MW
                employment_params[tech] = {
                    'installed_capacity': params['installed_capacity'],
                    'employment_per_mw': params['employment_per_mw'],
                    'display_name': params['display_name'],
                    'unit_system': 'MW'
                }
        
        return employment_params, config
    except Exception as e:
        print(f"Warning: Failed to read config file: {e}, using default parameters")
        return get_default_employment_parameters(), {}

def get_default_employment_parameters():
    """
    获取默认的各行业参数
    
    Returns:
    --------
    dict
        包含各行业参数的字典
    """
    # 默认参数值
    # 单位：装机容量(MW)，单位就业人数(人/MW)
    employment_params = {
        'Aluminum': {
            'installed_capacity': 18.8,  # MW - 铝冶炼厂装机容量
            'employment_per_mw': 15.4,      # 人/MW - 每MW装机对应的就业人数
            'display_name': 'Aluminum Smelter'
        },
        'Coal': {
            'installed_capacity': 328.0,  # MW - 煤电装机容量
            'employment_per_mw': 0.8,      # 人/MW - 每MW装机对应的就业人数
            'display_name': 'Coal Power'
        },
        'Gas': {
            'installed_capacity': 100.80,  # MW - 天然气发电装机容量
            'employment_per_mw': 0.2,      # 人/MW - 每MW装机对应的就业人数
            'display_name': 'Gas Power'
        }
    }
    
    return employment_params

def get_scenario_specific_parameters(scenario_type):
    """
    获取特定场景的参数
    
    Parameters:
    -----------
    scenario_type : str
        场景类型 ('20p' 或 'non_flexible')
    
    Returns:
    --------
    dict
        包含各行业参数的字典
    """
    if scenario_type == 'non_flexible':
        # non_flexible场景的特殊参数
        employment_params = {
            'Aluminum': {
                'installed_capacity': 11.7,  # MW - 铝冶炼厂装机容量
                'employment_per_mw': 15.4,      # 人/MW - 每MW装机对应的就业人数
                'display_name': 'Aluminum Smelter'
            },
            'Coal': {
                'installed_capacity': 343.7,  # MW - 煤电装机容量
                'employment_per_mw': 0.8,      # 人/MW - 每MW装机对应的就业人数
                'display_name': 'Coal Power'
            },
            'Gas': {
                'installed_capacity': 100.80,  # MW - 天然气发电装机容量
                'employment_per_mw': 0.2,      # 人/MW - 每MW装机对应的就业人数
                'display_name': 'Gas Power'
            }
        }
    else:
        # 20p场景使用默认参数
        employment_params = get_default_employment_parameters()
    
    return employment_params

def calculate_monthly_employment(capacity_factors, employment_params):
    """
    基于平均容量因子计算逐月就业人数
    
    Parameters:
    -----------
    capacity_factors : pd.DataFrame
        平均容量因子数据，行为月份，列为技术类型
    employment_params : dict
        各行业的装机容量和单位就业人数参数
    
    Returns:
    --------
    pd.DataFrame
        逐月就业人数数据
    """
    employment_data = {}
    
    for tech, params in employment_params.items():
        if tech in capacity_factors.columns:
            # 计算就业人数：平均容量因子 × 装机容量 × 单位就业人数
            if params.get('unit_system') == 'GW':
                # 使用GW单位系统
                monthly_employment = (capacity_factors[tech] * 
                                    params['installed_capacity'] * 
                                    params['employment_per_gw'])
            else:
                # 使用MW单位系统
                monthly_employment = (capacity_factors[tech] * 
                                    params['installed_capacity'] * 
                                    params['employment_per_mw'])
            employment_data[params['display_name']] = monthly_employment
        else:
            print(f"Warning: No capacity factor data found for {tech}")
    
    return pd.DataFrame(employment_data, index=capacity_factors.index)

def plot_employment_comparison(employment_data_20p, employment_data_non_flexible, 
                             output_file=None, config=None):
    """
    绘制两个场景的就业人数对比图（上下布局）
    
    Parameters:
    -----------
    employment_data_20p : pd.DataFrame
        20p场景的就业人数数据
    employment_data_non_flexible : pd.DataFrame
        non_flexible场景的就业人数数据
    output_file : str, optional
        输出图片文件路径
    config : dict, optional
        配置字典，包含颜色等设置
    """
    if employment_data_20p.empty or employment_data_non_flexible.empty:
        print("Warning: No employment data to plot")
        return
    
    # 设置绘图样式
    set_plot_style()
    
    # 从配置获取图表设置
    if config and 'plot' in config:
        plot_config = config['plot']
        figsize = plot_config.get('figure_size', [12, 9])
        dpi = plot_config.get('dpi', 150)
        font_size = plot_config.get('font_size', 24)
        legend_font_size = plot_config.get('legend_font_size', 24)
        title_font_size = plot_config.get('title_font_size', 24)
        axis_font_size = plot_config.get('axis_font_size', 24)
    else:
        figsize = [12, 9]
        dpi = 150
        font_size = 24
        legend_font_size = 24
        title_font_size = 24
        axis_font_size = 24
    
    # 创建上下两个子图
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=[12, 9], sharex=False)
    
    # 从配置获取颜色设置
    if config and 'colors' in config:
        colors = config['colors']
    else:
        # 默认颜色
        colors = {
            'Aluminum Smelter': '#FF69B4',    # Hot pink
            'Coal Power': '#000000',          # Black
            'Gas Power': '#FF0000'            # Red
        }
    
    # 图例标签映射（无论是否使用配置文件都定义）
    legend_labels = {
        'Aluminum Smelter': 'Aluminum smelters',
        'Coal Power': 'Coal power plants',
        'Gas Power': 'Gas power plants'
    }
    
    # 绘制20p场景（上图）
    plot_single_scenario(ax1, employment_data_20p, colors, "Maintaining 36% Overcapacity", 
                        font_size, legend_font_size, axis_font_size, show_legend=False,
                        legend_labels=legend_labels)
    
    # 绘制non_flexible场景（下图）
    plot_single_scenario(ax2, employment_data_non_flexible, colors, "Decommissioning All Overcapacity", 
                        font_size, legend_font_size, axis_font_size, show_legend=False,
                        legend_labels=legend_labels)
    
    # 设置统一的y轴范围
    ax1.set_ylim(0, 400)
    ax2.set_ylim(0, 400)
    
    # 设置子图标题
    ax1.set_title('Maintaining 36% Overcapacity', fontsize=title_font_size, pad=20)
    ax2.set_title('Decommissioning All Overcapacity', fontsize=title_font_size, pad=20)
    
    # 创建统一的图例（放在外面下方，不带框）
    # 获取第一个子图的图例元素
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(1, -0.05), 
              ncol=len(handles), fontsize=legend_font_size, frameon=False)
    
    # 调整布局，为底部图例留出空间
    plt.tight_layout()
    plt.subplots_adjust(right=2)
    
    # 保存图表
    if output_file is None:
        output_file = "employment_scenario_comparison.png"
    
    fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.show()
    plt.close()
    
    print(f"Employment comparison plot saved to: {output_file}")

def plot_single_scenario(ax, employment_data, colors, scenario_name, 
                        font_size, legend_font_size, axis_font_size, show_legend=True, 
                        legend_labels=None):
    """
    绘制单个场景的就业数据
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        绘图轴
    employment_data : pd.DataFrame
        就业人数数据
    colors : dict
        颜色字典
    scenario_name : str
        场景名称
    font_size : int
        字体大小
    legend_font_size : int
        图例字体大小
    axis_font_size : int
        坐标轴字体大小
    show_legend : bool
        是否显示图例
    """
    # 重新排序列，将Aluminum Smelter放在最上面
    column_order = []
    if 'Aluminum Smelter' in employment_data.columns:
        column_order.append('Aluminum Smelter')
    for col in employment_data.columns:
        if col != 'Aluminum Smelter':
            column_order.append(col)
    
    # 按新顺序重新排列数据
    employment_data_ordered = employment_data[column_order]
    
    # 计算累积值用于叠加面积图
    cumulative_data = employment_data_ordered.cumsum(axis=1)
    
    # 绘制叠加面积图
    for i, industry in enumerate(employment_data_ordered.columns):
        months = employment_data_ordered.index
        # 使用图例标签映射，如果没有映射则使用原名称
        label = legend_labels.get(industry, industry) if legend_labels else industry
        
        if i == 0:
            # 第一个行业：从0开始绘制
            values = cumulative_data[industry].values
            ax.fill_between(months, 0, values, color=colors.get(industry, '#808080'), 
                           alpha=0.7, label=label)
        else:
            # 后续行业：从上一个行业的累积值开始绘制
            prev_industry = employment_data_ordered.columns[i-1]
            prev_values = cumulative_data[prev_industry].values
            values = cumulative_data[industry].values
            ax.fill_between(months, prev_values, values, color=colors.get(industry, '#808080'), 
                           alpha=0.7, label=label)
    
    # 设置图表属性
    ax.set_ylabel('Needed workers', fontsize=axis_font_size)
    ax.set_xlim(1.0, 12.0)
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                        'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], fontsize=font_size)
    
    # 设置y轴刻度标签，添加"k"后缀
    ax.set_yticks(range(0, 401, 100))
    ax.set_yticklabels([f'{i}k' for i in range(0, 401, 100)], fontsize=font_size)
    
    # 确保x轴刻度标签可见
    ax.tick_params(axis='x', labelsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    
    ax.grid(True, alpha=0.3)
    
    # 根据参数决定是否显示图例
    if show_legend:
        ax.legend(loc='best', fontsize=legend_font_size)

def calculate_scenario_differences(employment_data_20p, employment_data_non_flexible):
    """
    计算两个场景的差异统计
    
    Parameters:
    -----------
    employment_data_20p : pd.DataFrame
        20p场景的就业人数数据
    employment_data_non_flexible : pd.DataFrame
        non_flexible场景的就业人数数据
    
    Returns:
    --------
    dict
        包含差异统计的字典
    """
    differences = {}
    
    # 确保两个数据框有相同的列
    common_columns = set(employment_data_20p.columns) & set(employment_data_non_flexible.columns)
    
    for industry in common_columns:
        data_20p = employment_data_20p[industry]
        data_non_flexible = employment_data_non_flexible[industry]
        
        # 计算差异
        diff = data_20p - data_non_flexible
        diff_percent = (diff / data_non_flexible) * 100
        
        differences[industry] = {
            'absolute_diff': diff,
            'percent_diff': diff_percent,
            'avg_20p': data_20p.mean(),
            'avg_non_flexible': data_non_flexible.mean(),
            'total_20p': data_20p.sum(),
            'total_non_flexible': data_non_flexible.sum(),
            'max_diff': diff.max(),
            'min_diff': diff.min(),
            'avg_diff': diff.mean(),
            'avg_percent_diff': diff_percent.mean()
        }
    
    return differences

def print_comparison_statistics(differences, employment_params_20p, employment_params_non_flexible):
    """
    打印对比统计信息
    
    Parameters:
    -----------
    differences : dict
        包含差异统计的字典
    employment_params_20p : dict
        20p场景的就业参数
    employment_params_non_flexible : dict
        non_flexible场景的就业参数
    """
    print(f"\nEmployment Scenario Comparison Statistics")
    print("=" * 80)
    
    # 打印场景参数信息
    print(f"\nScenario Parameters:")
    print(f"20p Scenario:")
    for tech, params in employment_params_20p.items():
        print(f"  {params['display_name']}: {params['installed_capacity']} MW")
    
    print(f"Non_flexible Scenario:")
    for tech, params in employment_params_non_flexible.items():
        print(f"  {params['display_name']}: {params['installed_capacity']} MW")
    
    print(f"\nEmployment Statistics:")
    for industry, stats in differences.items():
        print(f"\n{industry}:")
        print(f"  Average Employment - 20p: {stats['avg_20p']:.1f}k")
        print(f"  Average Employment - Non-flexible: {stats['avg_non_flexible']:.1f}k")
        print(f"  Average Difference: {stats['avg_diff']:.1f}k ({stats['avg_percent_diff']:.1f}%)")
        print(f"  Total Annual - 20p: {stats['total_20p']:.1f}k")
        print(f"  Total Annual - Non-flexible: {stats['total_non_flexible']:.1f}k")
        print(f"  Max Monthly Difference: {stats['max_diff']:.1f}k")
        print(f"  Min Monthly Difference: {stats['min_diff']:.1f}k")

def save_comparison_data(employment_data_20p, employment_data_non_flexible, 
                        differences, output_file=None):
    """
    保存对比数据到CSV文件
    
    Parameters:
    -----------
    employment_data_20p : pd.DataFrame
        20p场景的就业人数数据
    employment_data_non_flexible : pd.DataFrame
        non_flexible场景的就业人数数据
    differences : dict
        差异统计字典
    output_file : str, optional
        输出CSV文件路径
    """
    if output_file is None:
        output_file = "employment_scenario_comparison.csv"
    
    # 创建对比数据框
    comparison_data = {}
    
    # 添加20p场景数据
    for industry in employment_data_20p.columns:
        comparison_data[f"{industry}_20p"] = employment_data_20p[industry]
    
    # 添加non_flexible场景数据
    for industry in employment_data_non_flexible.columns:
        comparison_data[f"{industry}_non_flexible"] = employment_data_non_flexible[industry]
    
    # 添加差异数据
    for industry, stats in differences.items():
        comparison_data[f"{industry}_difference"] = stats['absolute_diff']
        comparison_data[f"{industry}_percent_diff"] = stats['percent_diff']
    
    comparison_df = pd.DataFrame(comparison_data, index=employment_data_20p.index)
    comparison_df.to_csv(output_file, encoding='utf-8-sig')
    print(f"Comparison data saved to: {output_file}")

def find_scenario_files(base_dir="results/monthly_capacity_factors"):
    """
    查找包含_20p和_non_flexible后缀的CSV文件
    
    Parameters:
    -----------
    base_dir : str
        基础目录路径
    
    Returns:
    --------
    tuple
        (file_20p, file_non_flexible) 两个文件路径
    """
    if not os.path.exists(base_dir):
        raise FileNotFoundError(f"Directory not found: {base_dir}")
    
    csv_files = [f for f in os.listdir(base_dir) if f.endswith('.csv')]
    
    file_20p = None
    file_non_flexible = None
    
    for csv_file in csv_files:
        if '_20p.csv' in csv_file:
            file_20p = os.path.join(base_dir, csv_file)
        elif '_non_flexible.csv' in csv_file:
            file_non_flexible = os.path.join(base_dir, csv_file)
    
    if file_20p is None:
        raise FileNotFoundError("No _20p.csv file found")
    if file_non_flexible is None:
        raise FileNotFoundError("No _non_flexible.csv file found")
    
    return file_20p, file_non_flexible

def compare_employment_scenarios(file_20p=None, file_non_flexible=None, 
                               output_dir=None, config_file=None):
    """
    对比两个场景的就业人数情况的主函数
    
    Parameters:
    -----------
    file_20p : str, optional
        20p场景的CSV文件路径
    file_non_flexible : str, optional
        non_flexible场景的CSV文件路径
    output_dir : str, optional
        输出目录
    config_file : str, optional
        配置文件路径
    """
    # 如果没有指定文件，自动查找
    if file_20p is None or file_non_flexible is None:
        try:
            file_20p, file_non_flexible = find_scenario_files()
            print(f"Found 20p file: {os.path.basename(file_20p)}")
            print(f"Found non_flexible file: {os.path.basename(file_non_flexible)}")
        except FileNotFoundError as e:
            print(f"Error: {e}")
            return
    
    # 加载数据
    print("\nLoading data...")
    capacity_factors_20p, _ = load_csv_data(file_20p)
    capacity_factors_non_flexible, _ = load_csv_data(file_non_flexible)
    
    if capacity_factors_20p.empty or capacity_factors_non_flexible.empty:
        print("Warning: No capacity factor data found in CSV files")
        return
    
    # 获取就业参数和配置
    _, config = load_employment_config(config_file)
    
    # 获取场景特定的参数
    employment_params_20p = get_scenario_specific_parameters('20p')
    employment_params_non_flexible = get_scenario_specific_parameters('non_flexible')
    
    # 计算就业人数
    print("Calculating employment...")
    print("20p scenario parameters:")
    for tech, params in employment_params_20p.items():
        print(f"  {params['display_name']}: {params['installed_capacity']} MW")
    
    print("Non_flexible scenario parameters:")
    for tech, params in employment_params_non_flexible.items():
        print(f"  {params['display_name']}: {params['installed_capacity']} MW")
    
    employment_data_20p = calculate_monthly_employment(capacity_factors_20p, employment_params_20p)
    employment_data_non_flexible = calculate_monthly_employment(capacity_factors_non_flexible, employment_params_non_flexible)
    
    if employment_data_20p.empty or employment_data_non_flexible.empty:
        print("Warning: No employment data calculated")
        return
    
    # 设置输出目录
    if output_dir is None:
        if config and 'output' in config:
            output_dir = config['output'].get('directory', 'results/employment_analysis')
        else:
            output_dir = "results/employment_analysis"
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 生成输出文件名
    plot_file = os.path.join(output_dir, "employment_scenario_comparison.png")
    csv_file = os.path.join(output_dir, "employment_scenario_comparison.csv")
    
    # 计算差异统计
    print("Calculating differences...")
    differences = calculate_scenario_differences(employment_data_20p, employment_data_non_flexible)
    
    # 绘制对比图表
    print("Creating comparison plot...")
    plot_employment_comparison(employment_data_20p, employment_data_non_flexible, 
                             plot_file, config)
    
    # 保存对比数据
    save_comparison_data(employment_data_20p, employment_data_non_flexible, 
                        differences, csv_file)
    
    # 打印统计信息
    print_comparison_statistics(differences, employment_params_20p, employment_params_non_flexible)

def main():
    """
    主函数，处理命令行参数
    """
    parser = argparse.ArgumentParser(description='Compare employment between 20p and non_flexible scenarios')
    parser.add_argument('--file_20p', help='20p scenario CSV file path')
    parser.add_argument('--file_non_flexible', help='Non_flexible scenario CSV file path')
    parser.add_argument('-o', '--output', help='Output directory')
    parser.add_argument('-c', '--config', help='Config file path')
    
    args = parser.parse_args()
    
    try:
        compare_employment_scenarios(args.file_20p, args.file_non_flexible, 
                                   args.output, args.config)
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
