#!/usr/bin/env python3
"""
绘制最优点的箱线图
从 optimal_points_distribution_data_latest.csv 文件读取数据并绘制箱线图
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 设置matplotlib中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def load_optimal_points_data(csv_file='results/optimal_points_analysis/optimal_points_distribution_data_latest.csv'):
    """
    从CSV文件加载最优点数据
    
    Parameters:
    -----------
    csv_file : str
        CSV文件路径
        
    Returns:
    --------
    pd.DataFrame
        加载的数据
    """
    try:
        df = pd.read_csv(csv_file)
        logger.info(f"成功加载数据: {len(df)} 个数据点")
        logger.info(f"年份: {sorted(df['year'].unique())}")
        logger.info(f"市场: {sorted(df['market'].unique())}")
        logger.info(f"灵活性: {sorted(df['flexibility'].unique())}")
        return df
    except Exception as e:
        logger.error(f"加载数据失败: {str(e)}")
        return None

def plot_capacity_boxplot(df, output_dir='results/optimal_points_analysis'):
    """
    绘制产能箱线图
    
    Parameters:
    -----------
    df : pd.DataFrame
        数据框
    output_dir : str
        输出目录
    """
    # 按年份分组
    years = sorted(df['year'].unique())
    
    # 准备数据
    capacity_data = []
    year_labels = []
    
    for year in years:
        year_data = df[df['year'] == year]
        capacities = year_data['capacity'].values
        capacity_data.append(capacities)
        year_labels.append(f'{year}')
    
    # 创建图形
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # 设置颜色
    colors = plt.cm.viridis(np.linspace(0, 1, len(years)))
    year_colors = dict(zip(years, colors))
    
    # 绘制箱线图
    bp = ax.boxplot(capacity_data, labels=year_labels, patch_artist=True, 
                   boxprops=dict(alpha=0.7), medianprops=dict(color='black', linewidth=2))
    
    # 为每个箱子设置颜色
    for patch, year in zip(bp['boxes'], years):
        patch.set_facecolor(year_colors[year])
    
    # 设置标题和标签
    ax.set_title('Optimal Aluminum Smelting Capacity Distribution by Year', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Year', fontsize=14, fontweight='bold')
    ax.set_ylabel('Aluminum Smelting Capacity (10,000 tons/year)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # 格式化y轴显示整数
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x)}'))
    
    # 添加年度需求线
    demand_by_year = {
        2030: 29.0241717,
        2040: 15.0817033,
        2050: 11.6668363,
    }
    
    demand_colors = {2030: 'red', 2040: 'orange', 2050: 'purple'}
    
    for i, year in enumerate(years, 1):
        if year in demand_by_year:
            demand = demand_by_year[year]
            ax.axhline(y=demand, xmin=(i-1)/len(years), xmax=i/len(years), 
                      color=year_colors[year], linestyle='--', linewidth=3, alpha=0.8)
    
    # 添加图例
    legend_elements = []
    for year in years:
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=year_colors[year], alpha=0.7, label=f'{year}'))
    
    # 添加需求线到图例
    for year, demand in demand_by_year.items():
        if year in years:
            legend_elements.append(plt.Line2D([0], [0], color=year_colors[year], linestyle='--', linewidth=3, 
                                            label=f'{year} Demand: {demand:.0f} Mt/year'))
    
    ax.legend(handles=legend_elements, loc='lower center', fontsize=12)
    
    plt.tight_layout()
    
    # 保存图形
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    plot_file = output_path / "optimal_points_capacity_boxplot.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"产能箱线图已保存到: {plot_file}")
    
    return fig

def plot_net_value_boxplot(df, output_dir='results/optimal_points_analysis'):
    """
    绘制净价值箱线图
    
    Parameters:
    -----------
    df : pd.DataFrame
        数据框
    output_dir : str
        输出目录
    """
    # 按年份分组
    years = sorted(df['year'].unique())
    
    # 准备数据
    net_value_data = []
    year_labels = []
    
    for year in years:
        year_data = df[df['year'] == year]
        net_values = year_data['net_value'].values
        net_value_data.append(net_values)
        year_labels.append(f'{year}')
    
    # 创建图形
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.5))
    
    # 设置颜色
    colors = plt.cm.viridis(np.linspace(0, 1, len(years)))
    year_colors = dict(zip(years, colors))
    
    # 绘制箱线图
    bp = ax.boxplot(net_value_data, labels=year_labels, patch_artist=True, 
                   boxprops=dict(alpha=0.7), medianprops=dict(color='black', linewidth=2))
    
    # 为每个箱子设置颜色
    for patch, year in zip(bp['boxes'], years):
        patch.set_facecolor(year_colors[year])
    
    # 设置标题和标签
    ax.set_title('Optimal Net System Value Distribution by Year', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Year', fontsize=14, fontweight='bold')
    ax.set_ylabel('Net System Value (Billion CNY)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # 格式化y轴显示整数
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x)}'))
    
    # 添加图例
    legend_elements = []
    for year in years:
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=year_colors[year], alpha=0.7, label=f'{year}'))
    
    ax.legend(handles=legend_elements, loc='lower center', fontsize=12)
    
    plt.tight_layout()
    
    # 保存图形
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    plot_file = output_path / "optimal_points_net_value_boxplot.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"净价值箱线图已保存到: {plot_file}")
    
    # Show plot
    # plt.show()
    
    return fig

def plot_excess_ratio_boxplot(df, output_dir='results/optimal_points_analysis'):
    """
    绘制过剩比例箱线图
    
    Parameters:
    -----------
    df : pd.DataFrame
        数据框
    output_dir : str
        输出目录
    """
    # 按年份分组
    years = sorted(df['year'].unique())
    
    # 准备数据
    excess_ratio_data = []
    year_labels = []
    
    for year in years:
        year_data = df[df['year'] == year]
        excess_ratios = year_data['excess_ratio'].values
        excess_ratio_data.append(excess_ratios)
        year_labels.append(f'{year}')
    
    # 创建图形
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # 设置颜色
    colors = plt.cm.viridis(np.linspace(0, 1, len(years)))
    year_colors = dict(zip(years, colors))
    
    # 绘制箱线图
    bp = ax.boxplot(excess_ratio_data, labels=year_labels, patch_artist=True, 
                   boxprops=dict(alpha=0.7), medianprops=dict(color='black', linewidth=2))
    
    # 为每个箱子设置颜色
    for patch, year in zip(bp['boxes'], years):
        patch.set_facecolor(year_colors[year])
    
    # 设置标题和标签
    ax.set_title('Optimal Excess Ratio Distribution by Year', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Year', fontsize=14, fontweight='bold')
    ax.set_ylabel('Excess Ratio', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # 格式化y轴显示百分比
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.1%}'))
    
    # 添加图例
    legend_elements = []
    for year in years:
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=year_colors[year], alpha=0.7, label=f'{year}'))
    
    ax.legend(handles=legend_elements, loc='lower center', fontsize=12)
    
    plt.tight_layout()
    
    # 保存图形
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    plot_file = output_path / "optimal_points_excess_ratio_boxplot.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"过剩比例箱线图已保存到: {plot_file}")
    
    return fig

def plot_combined_boxplot(df, output_dir='results/optimal_points_analysis'):
    """
    绘制组合箱线图（产能和净价值）
    
    Parameters:
    -----------
    df : pd.DataFrame
        数据框
    output_dir : str
        输出目录
    """
    # 按年份分组
    years = sorted(df['year'].unique())
    
    # 准备数据
    capacity_data = []
    net_value_data = []
    year_labels = []
    
    for year in years:
        year_data = df[df['year'] == year]
        capacities = year_data['capacity'].values
        net_values = year_data['net_value'].values
        
        capacity_data.append(capacities)
        net_value_data.append(net_values)
        year_labels.append(f'{year}')
    
    # 创建图形 - 改为上下排列，压扁子图
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8.5))
    
    # 设置颜色
    colors = plt.cm.viridis(np.linspace(0, 1, len(years)))
    year_colors = dict(zip(years, colors))
    
    # 上图：产能箱线图
    bp1 = ax1.boxplot(capacity_data, labels=year_labels, patch_artist=True, 
                     boxprops=dict(alpha=0.7), medianprops=dict(color='black', linewidth=2))
    
    for patch, year in zip(bp1['boxes'], years):
        patch.set_facecolor(year_colors[year])
    
    ax1.set_title('Optimal Aluminum Smelting Capacity Distribution by Year', fontsize=14, fontweight='bold', pad=15)
    # ax1.set_xlabel('Year', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Aluminum Smelting Capacity (Mt/year)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x)}'))
    
    # 添加年度需求线 - 使用与箱子相同的颜色
    demand_by_year = {
        2030: 29.0241717,
        2040: 15.0817033,
        2050: 11.6668363,
    }
    
    for i, year in enumerate(years, 1):
        if year in demand_by_year:
            demand = demand_by_year[year]
            # 使用与箱子相同的颜色
            ax1.axhline(y=demand, xmin=(i-1)/len(years), xmax=i/len(years), 
                       color=year_colors[year], linestyle='--', linewidth=3, alpha=0.8)
    
    # 计算每年的平均产能
    year_capacity_means = {}
    for i, year in enumerate(years):
        year_data = df[df['year'] == year]
        if not year_data.empty:
            year_capacity_means[year] = year_data['capacity'].mean()
    
    # 添加连接需求值的虚线
    demand_years = [year for year in years if year in demand_by_year]
    if len(demand_years) > 1:
        demand_values = [demand_by_year[year] for year in demand_years]
        # 将年份转换为x轴位置（1, 2, 3...）
        demand_x_positions = [years.index(year) + 1 for year in demand_years]
        ax1.plot(demand_x_positions, demand_values, 'k--', linewidth=2, alpha=0.6, 
                label='Demand Trend')
    
    # 添加连接平均产能的虚线
    if len(year_capacity_means) > 1:
        capacity_years = sorted(year_capacity_means.keys())
        capacity_values = [year_capacity_means[year] for year in capacity_years]
        # 将年份转换为x轴位置（1, 2, 3...）
        capacity_x_positions = [years.index(year) + 1 for year in capacity_years]
        ax1.plot(capacity_x_positions, capacity_values, 'b--', linewidth=2, alpha=0.6, 
                label='Average Capacity Trend')
    
    # 下图：净价值箱线图
    bp2 = ax2.boxplot(net_value_data, labels=year_labels, patch_artist=True, 
                     boxprops=dict(alpha=0.7), medianprops=dict(color='black', linewidth=2))
    
    for patch, year in zip(bp2['boxes'], years):
        patch.set_facecolor(year_colors[year])
    
    ax2.set_title('Optimal Net System Value Distribution by Year', fontsize=14, fontweight='bold', pad=15)
    ax2.set_xlabel('Year', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Net System Value (Billion CNY)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x)}'))
    
    # 添加图例
    legend_elements = []
    # for year in years:
    #     legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=year_colors[year], alpha=0.7, label=f'{year}'))
    
    # 添加需求线到图例
    for year, demand in demand_by_year.items():
        if year in years:
            legend_elements.append(plt.Line2D([0], [0], color=year_colors[year], linestyle='--', linewidth=3, 
                                            label=f'{year} Demand'))
    
    # 添加趋势线到图例
    if len(demand_years) > 1:
        legend_elements.append(plt.Line2D([0], [0], color='k', linestyle='--', linewidth=2, 
                                        label='Primary Aluminum Demand'))
    
    if len(year_capacity_means) > 1:
        legend_elements.append(plt.Line2D([0], [0], color='b', linestyle='--', linewidth=2, 
                                        label='Optimal Smelting Capacity'))
    
    # 只在图形外部添加一个图例
    fig.legend(handles=legend_elements, loc='lower center',
              ncol=len(legend_elements), fontsize=10)
    
    # 调整布局，为图例留出空间
    # plt.tight_layout(rect=[0, 0.1, 1, 1])
    
    # 保存图形
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    plot_file = output_path / "optimal_points_combined_boxplot.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"组合箱线图已保存到: {plot_file}")
    
    return fig

def print_data_summary(df):
    """
    打印数据摘要统计
    
    Parameters:
    -----------
    df : pd.DataFrame
        数据框
    """
    logger.info("=== 数据摘要统计 ===")
    logger.info(f"总数据点数: {len(df)}")
    logger.info(f"年份范围: {df['year'].min()} - {df['year'].max()}")
    logger.info(f"包含年份: {sorted(df['year'].unique())}")
    logger.info(f"市场类型: {sorted(df['market'].unique())}")
    logger.info(f"灵活性类型: {sorted(df['flexibility'].unique())}")
    
    logger.info("\n=== 按年份统计 ===")
    for year in sorted(df['year'].unique()):
        year_data = df[df['year'] == year]
        logger.info(f"{year}年:")
        logger.info(f"  数据点数: {len(year_data)}")
        logger.info(f"  产能范围: {year_data['capacity'].min():.1f} - {year_data['capacity'].max():.1f} 万吨/年")
        logger.info(f"  净价值范围: {year_data['net_value'].min():.2f} - {year_data['net_value'].max():.2f} 十亿人民币")
        logger.info(f"  过剩比例范围: {year_data['excess_ratio'].min():.1%} - {year_data['excess_ratio'].max():.1%}")

def main():
    """
    主函数
    """
    # 数据文件路径
    csv_file = 'results/optimal_points_analysis/optimal_points_distribution_data_latest.csv'
    
    # 加载数据
    df = load_optimal_points_data(csv_file)
    if df is None:
        logger.error("无法加载数据，程序退出")
        return
    
    # 打印数据摘要
    print_data_summary(df)
    
    # 绘制组合箱线图（上下排列）
    logger.info("开始绘制组合箱线图...")
    plot_combined_boxplot(df)
    
    logger.info("箱线图绘制完成！")

if __name__ == "__main__":
    main()
