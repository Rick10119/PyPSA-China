"""
基于平均容量因子数据计算smelter、coal、gas行业逐月就业人数的脚本
计算公式：就业人数 = 平均容量因子 × 装机容量 × 单位装机就业人数
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import argparse
import yaml

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
                    'font.size': 18,
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
        print("使用月度平均容量因子数据计算就业")
        
    else:
        # 兼容旧格式：只有Capacity_Factor列
        capacity_cols = [col for col in df.columns if 'Capacity_Factor' in col and 'Max' not in col and 'Avg' not in col]
        capacity_factors = df[capacity_cols] if capacity_cols else pd.DataFrame()
        capacity_factors.columns = [col.replace('_Capacity_Factor', '') for col in capacity_factors.columns]
        print("使用旧格式容量因子数据计算就业")
    
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
    dict
        包含各行业参数的字典
    """
    if config_file is None:
        config_file = os.path.join(os.path.dirname(__file__), 'employment_config.yaml')
    
    if not os.path.exists(config_file):
        print(f"Warning: Config file {config_file} not found, using default parameters")
        return get_default_employment_parameters()
    
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

def plot_employment_trends(employment_data, output_file=None, title_suffix="", config=None):
    """
    绘制就业人数趋势图
    
    Parameters:
    -----------
    employment_data : pd.DataFrame
        就业人数数据
    output_file : str, optional
        输出图片文件路径
    title_suffix : str, optional
        图表标题后缀
    config : dict, optional
        配置字典，包含颜色等设置
    """
    if employment_data.empty:
        print("Warning: No employment data to plot")
        return
    
    # 设置绘图样式
    set_plot_style()
    
    # 从配置获取图表设置
    if config and 'plot' in config:
        plot_config = config['plot']
        figsize = plot_config.get('figure_size', [12, 8])
        dpi = plot_config.get('dpi', 150)
        font_size = plot_config.get('font_size', 18)
        legend_font_size = plot_config.get('legend_font_size', 20)
        title_font_size = plot_config.get('title_font_size', 24)
        axis_font_size = plot_config.get('axis_font_size', 24)
    else:
        figsize = [12, 8]
        dpi = 150
        font_size = 18
        legend_font_size = 20
        title_font_size = 24
        axis_font_size = 24
    
    # 创建图表
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
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
    
    # 绘制各行业就业人数趋势（叠加面积图）
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
        if i == 0:
            # 第一个行业：从0开始绘制
            values = cumulative_data[industry].values
            ax.fill_between(months, 0, values, color=colors.get(industry, '#808080'), 
                           alpha=0.7, label=industry)
        else:
            # 后续行业：从上一个行业的累积值开始绘制
            prev_industry = employment_data_ordered.columns[i-1]
            prev_values = cumulative_data[prev_industry].values
            values = cumulative_data[industry].values
            ax.fill_between(months, prev_values, values, color=colors.get(industry, '#808080'), 
                           alpha=0.7, label=industry)
    
    # 添加边界线使图形更清晰
    for industry in employment_data_ordered.columns:
        months = employment_data_ordered.index
        values = employment_data_ordered[industry].values
        color = colors.get(industry, '#808080')
        ax.plot(months, values, color=color, linewidth=1, alpha=0.8)
    
    # 设置图表属性
    ax.set_ylabel('Cumulative Employment (thousands)', fontsize=axis_font_size)
    ax.set_xlabel('Month', fontsize=axis_font_size)
    ax.set_title(f'Monthly Employment by Industry (Stacked Area){title_suffix}', fontsize=title_font_size)
    ax.set_xlim(1.0, 12.0)
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                        'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], fontsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=legend_font_size)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图表
    if output_file is None:
        output_file = f"employment_trends{title_suffix.replace(' ', '_')}.png"
    
    fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.show()
    plt.close()
    
    print(f"Employment trend plot saved to: {output_file}")

def save_employment_data(employment_data, output_file=None, title_suffix=""):
    """
    保存就业数据到CSV文件
    
    Parameters:
    -----------
    employment_data : pd.DataFrame
        就业人数数据
    output_file : str, optional
        输出CSV文件路径
    title_suffix : str, optional
        文件名后缀
    """
    if output_file is None:
        output_file = f"monthly_employment{title_suffix.replace(' ', '_')}.csv"
    
    employment_data.to_csv(output_file, encoding='utf-8-sig')
    print(f"Employment data saved to: {output_file}")

def print_employment_statistics(employment_data, employment_params):
    """
    打印就业统计信息
    
    Parameters:
    -----------
    employment_data : pd.DataFrame
        就业人数数据
    employment_params : dict
        各行业参数
    """
    print(f"\nMonthly Employment Statistics by Industry (Based on Average Capacity Factors)")
    print("=" * 80)
    
    for industry in employment_data.columns:
        data = employment_data[industry]
        avg_employment = data.mean()
        max_employment = data.max()
        min_employment = data.min()
        total_annual = data.sum()
        
        # 统一显示千人单位
        print(f"{industry:15s}: Avg={avg_employment:.1f}k, Max={max_employment:.1f}k, "
              f"Min={min_employment:.1f}k, Annual Total={total_annual:.1f}k")
    
    print(f"\nIndustry Parameters")
    print("=" * 40)
    for tech, params in employment_params.items():
        if tech in ['Aluminum', 'Coal', 'Gas']:
            unit_system = params.get('unit_system', 'MW')
            if unit_system == 'GW':
                # 显示GW单位
                print(f"{params['display_name']:15s}: Capacity={params['installed_capacity']:.1f}GW, "
                      f"Employment={params['employment_per_gw']:.1f}k/GW")
            else:
                # 显示MW单位
                print(f"{params['display_name']:15s}: Capacity={params['installed_capacity']:.0f}MW, "
                      f"Employment={params['employment_per_mw']:.1f}/MW")

def calculate_employment_from_csv(csv_file, output_dir=None, title_suffix="", config_file=None):
    """
    从CSV文件基于平均容量因子计算就业人数的主函数
    
    Parameters:
    -----------
    csv_file : str
        CSV文件路径
    output_dir : str, optional
        输出目录
    title_suffix : str, optional
        标题后缀
    config_file : str, optional
        配置文件路径
    """
    # 加载数据
    capacity_factors, load_factors = load_csv_data(csv_file)
    
    if capacity_factors.empty:
        print("Warning: No capacity factor data found in CSV file")
        return
    
    # 获取就业参数和配置
    employment_params, config = load_employment_config(config_file)
    
    # 计算就业人数
    employment_data = calculate_monthly_employment(capacity_factors, employment_params)
    
    if employment_data.empty:
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
    base_name = os.path.splitext(os.path.basename(csv_file))[0]
    plot_file = os.path.join(output_dir, f"{base_name}_employment_plot.png")
    csv_file_out = os.path.join(output_dir, f"{base_name}_employment_data.csv")
    
    # 绘制图表
    plot_employment_trends(employment_data, plot_file, title_suffix, config)
    
    # 保存数据
    save_employment_data(employment_data, csv_file_out, title_suffix)
    
    # 打印统计信息
    print_employment_statistics(employment_data, employment_params)

def main():
    """
    主函数，处理命令行参数
    """
    parser = argparse.ArgumentParser(description='Calculate employment from average capacity factor data')
    parser.add_argument('csv_file', help='CSV file path')
    parser.add_argument('-o', '--output', help='Output directory')
    parser.add_argument('-t', '--title', default='', help='Title suffix')
    parser.add_argument('-c', '--config', help='Config file path')
    
    args = parser.parse_args()
    
    try:
        calculate_employment_from_csv(args.csv_file, args.output, args.title, args.config)
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    # 如果直接运行脚本，自动查找并处理所有可用的CSV文件
    if len(os.sys.argv) == 1:
        # 查找所有可用的CSV文件
        csv_dir = "results/monthly_capacity_factors"
        if os.path.exists(csv_dir):
            csv_files = [f for f in os.listdir(csv_dir) if f.endswith('.csv')]
            csv_files.sort()  # 按文件名排序
            
            if csv_files:
                print(f"Found {len(csv_files)} CSV files:")
                for csv_file in csv_files:
                    print(f"  - {csv_file}")
                
                print("\nStarting file processing...")
                for csv_file in csv_files:
                    csv_path = os.path.join(csv_dir, csv_file)
                    print(f"\nProcessing file: {csv_file}")
                    
                    # 从文件名提取信息作为标题后缀
                    base_name = csv_file.replace('.csv', '')
                    if '_' in base_name:
                        parts = base_name.split('_')
                        if len(parts) >= 3:
                            year = parts[2] if parts[2].isdigit() else ""
                            province = parts[3] if len(parts) > 3 else ""
                            title_suffix = f" - {province} {year}" if province and year else f" - {base_name}"
                        else:
                            title_suffix = f" - {base_name}"
                    else:
                        title_suffix = f" - {base_name}"
                    
                    calculate_employment_from_csv(csv_path, title_suffix=title_suffix)
            else:
                print(f"No CSV files found in {csv_dir} directory")
                print("Please run plot_capacity_factors.py first to generate CSV files")
        else:
            print(f"Directory {csv_dir} does not exist")
            print("Please run plot_capacity_factors.py first to generate CSV files")
    else:
        exit(main())
