"""
从CSV文件读取月度容量因子数据并绘制图表的脚本
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import argparse

def set_plot_style():
    """
    设置绘图样式
    """
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial Unicode MS', 'Arial']
    plt.rcParams['axes.unicode_minus'] = False
    
    plt.style.use(['classic', 'seaborn-v0_8-whitegrid',
                   {'axes.grid': False, 'grid.linestyle': '--', 'grid.color': u'0.6',
                    'hatch.color': 'white',
                    'patch.linewidth': 0.5,
                    'font.size': 20,
                    'legend.fontsize': 'large',
                    'lines.linewidth': 1.5,
                    'pdf.fonttype': 42,
                    }])

def load_csv_data(csv_file):
    """
    从CSV文件加载数据
    
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
        raise FileNotFoundError(f"CSV文件不存在: {csv_file}")
    
    # 读取CSV文件
    df = pd.read_csv(csv_file, index_col='Month')
    
    # 分离容量因子和负荷因子数据
    capacity_cols = [col for col in df.columns if 'Capacity_Factor' in col]
    load_cols = [col for col in df.columns if 'Load_Factor' in col]
    
    capacity_factors = df[capacity_cols] if capacity_cols else pd.DataFrame()
    load_factors = df[load_cols] if load_cols else pd.DataFrame()
    
    # 重命名列，去掉后缀
    capacity_factors.columns = [col.replace('_Capacity_Factor', '') for col in capacity_factors.columns]
    load_factors.columns = [col.replace('_Load_Factor', '') for col in load_factors.columns]
    
    return capacity_factors, load_factors

def plot_capacity_factors_from_csv(csv_file, output_file=None, title_suffix=""):
    """
    从CSV文件绘制容量因子图表
    
    Parameters:
    -----------
    csv_file : str
        CSV文件路径
    output_file : str, optional
        输出图片文件路径
    title_suffix : str, optional
        图表标题后缀
    """
    # 加载数据
    capacity_factors, load_factors = load_csv_data(csv_file)
    
    if capacity_factors.empty and load_factors.empty:
        print("警告: CSV文件中没有找到容量因子或负荷因子数据")
        return
    
    # 设置绘图样式
    set_plot_style()
    
    # 创建图表
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # 定义颜色
    colors = {
        'Hydro': '#000080',      # Navy
        'Nuclear': '#800080',    # Purple
        'Coal': '#000000',       # Black
        'Gas': '#FF0000',        # Red
        'Wind': '#00BFFF',       # Deep sky blue
        'Solar': '#FFD700',      # Gold
        'Aluminum': '#FF69B4',   # Hot pink
        'Other': '#808080'       # Gray
    }
    
    load_colors = {
        'Electricity Load': '#1f77b4',    # Blue
        'Heating Load': '#ff7f0e',        # Orange
        'Aluminum Load': '#2ca02c'        # Green
    }
    
    # 绘制容量因子
    for tech in capacity_factors.columns:
        months = capacity_factors.index
        values = capacity_factors[tech].values
        color = colors.get(tech, '#000000')
        ax.plot(months, values, 'o-', color=color, 
                linewidth=2, markersize=6, label=tech)
    
    # 绘制负荷因子
    for load_type in load_factors.columns:
        months = load_factors.index
        values = load_factors[load_type].values
        color = load_colors.get(load_type, '#000000')
        ax.plot(months, values, 's--', color=color, 
                linewidth=2, markersize=6, label=load_type)
    
    # 设置图表属性
    ax.set_ylabel('Capacity/Load Factor (p.u.)', fontsize=20)
    ax.set_xlabel('Month', fontsize=20)
    ax.set_title(f'Monthly Capacity Factors & Load Factors{title_suffix}', 
                 fontsize=14, fontweight='bold')
    ax.set_xlim(1, 12)
    ax.set_ylim(0, 1.0)
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                         'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', ncol=2)
    
    # 添加数值标签
    for tech in capacity_factors.columns:
        months = capacity_factors.index
        values = capacity_factors[tech].values
        for month, value in zip(months, values):
            ax.annotate(f'{value:.2f}', (month, value), 
                       textcoords="offset points", xytext=(0,10), 
                       ha='center', fontsize=8)
    
    for load_type in load_factors.columns:
        months = load_factors.index
        values = load_factors[load_type].values
        for month, value in zip(months, values):
            ax.annotate(f'{value:.2f}', (month, value), 
                       textcoords="offset points", xytext=(0,10), 
                       ha='center', fontsize=8)
    
    plt.tight_layout()
    
    # 保存图表
    if output_file is None:
        output_file = csv_file.replace('.csv', '_plot.png')
    
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"图表已保存到: {output_file}")
    
    # 打印统计信息
    print(f"\n容量因子月度统计{title_suffix}")
    print("=" * 50)
    for tech in capacity_factors.columns:
        data = capacity_factors[tech]
        avg_cf = data.mean()
        max_cf = data.max()
        min_cf = data.min()
        print(f"{tech:15s}: 平均={avg_cf:.3f}, 最大={max_cf:.3f}, 最小={min_cf:.3f}")
    
    if not load_factors.empty:
        print(f"\n负荷因子月度统计{title_suffix}")
        print("=" * 50)
        for load_type in load_factors.columns:
            data = load_factors[load_type]
            avg_load = data.mean()
            max_load = data.max()
            min_load = data.min()
            print(f"{load_type:15s}: 平均={avg_load:.3f}, 最大={max_load:.3f}, 最小={min_load:.3f}")

def main():
    """
    主函数，处理命令行参数
    """
    parser = argparse.ArgumentParser(description='从CSV文件绘制容量因子图表')
    parser.add_argument('csv_file', help='CSV文件路径')
    parser.add_argument('-o', '--output', help='输出图片文件路径')
    parser.add_argument('-t', '--title', default='', help='图表标题后缀')
    
    args = parser.parse_args()
    
    try:
        plot_capacity_factors_from_csv(args.csv_file, args.output, args.title)
    except Exception as e:
        print(f"错误: {e}")
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
                print(f"找到 {len(csv_files)} 个CSV文件:")
                for csv_file in csv_files:
                    print(f"  - {csv_file}")
                
                print("\n开始处理文件...")
                for csv_file in csv_files:
                    csv_path = os.path.join(csv_dir, csv_file)
                    print(f"\n处理文件: {csv_file}")
                    
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
                    
                    plot_capacity_factors_from_csv(csv_path, title_suffix=title_suffix)
            else:
                print(f"在 {csv_dir} 目录中没有找到CSV文件")
                print("请先运行 plot_capacity_factors.py 生成CSV文件")
        else:
            print(f"目录 {csv_dir} 不存在")
            print("请先运行 plot_capacity_factors.py 生成CSV文件")
    else:
        exit(main())
