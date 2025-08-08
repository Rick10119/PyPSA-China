# SPDX-FileCopyrightText: : 2022 The PyPSA-China Authors
#
# SPDX-License-Identifier: MIT

"""
对比不同容量比例下电解铝导致的碳排放变化（相对于no_aluminum基准）
生成电解铝导致的月度碳排放变化图表
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

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def collect_all_emissions_data(capacity_ratios, file_type):
    """
    从所有容量比例目录中收集2050年排放数据
    
    Parameters:
    -----------
    capacity_ratios : list
        容量比例列表，如 ['no_aluminum', '55p', '60p', '70p', '80p', '90p', '100p']
    file_type : str
        文件类型
        
    Returns:
    --------
    dict
        按容量比例组织的数据
    """
    emissions_data = {}
    year = 2050
    
    for ratio in capacity_ratios:
        if ratio == 'no_aluminum':
            version_name = "0723.8H.5-no-aluminum"
        else:
            version_name = f"0723.8H.5-{ratio}"
        
        # 构建目录路径
        dir_pattern = f"postnetwork-ll-current+Neighbor-linear2050-{year}"
        version_dir = Path(f"results/version-{version_name}/summary/postnetworks/positive/{dir_pattern}")
        
        if version_dir.exists():
            # 读取emissions文件
            emissions_file = version_dir / f"{file_type}.csv"
            
            if emissions_file.exists():
                try:
                    # 读取数据
                    df = load_single_csv_file(emissions_file)
                    
                    if df is not None:
                        emissions_data[ratio] = df
                        logger.info(f"成功加载 {version_name} 的 {file_type} 数据")
                    else:
                        logger.warning(f"无法加载 {version_name} 的 {file_type} 数据")
                        
                except Exception as e:
                    logger.warning(f"加载 {version_name} 的 {file_type} 数据时出错: {str(e)}")
            else:
                logger.warning(f"在 {version_name} 目录中未找到 {file_type}.csv 文件")
        else:
            logger.warning(f"版本 {version_name} 的目录不存在")
    
    return emissions_data

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

def generate_aluminum_emissions_plot(emissions_data, output_dir):
    """
    生成电解铝导致的碳排放变化图表
    
    Parameters:
    -----------
    emissions_data : dict
        按容量比例组织的排放数据
    output_dir : Path
        输出目录
    """
    # 创建输出目录
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # 检查是否有55p基准数据
    if '55p' not in emissions_data:
        logger.error("缺少55p基准数据，无法计算变化量")
        return
    
    baseline_data = emissions_data['55p']
    
    # 计算相对于55p的变化量
    comparison_ratios = ['60p', '70p', '80p', '90p', '100p']
    
    # 月份标签
    month_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    
    # 计算每个容量比例的月度排放变化
    change_data = {}
    total_changes = {}
    
    for ratio in comparison_ratios:
        if ratio not in emissions_data:
            logger.warning(f"缺少 {ratio} 数据，跳过")
            continue
            
        comparison_data = emissions_data[ratio]
        change_data[ratio] = []
        total_changes[ratio] = 0
        
        # 计算总排放变化（煤炭+天然气）
        total_changes_monthly = []
        for month in month_labels:
            coal_indicator = f"coal_emissions_{month.lower()}"
            gas_indicator = f"gas_emissions_{month.lower()}"
            
            # 计算煤炭排放变化
            coal_baseline = 0
            coal_comparison = 0
            if coal_indicator in baseline_data.index:
                coal_baseline = baseline_data.loc[coal_indicator].iloc[0] if not pd.isna(baseline_data.loc[coal_indicator].iloc[0]) else 0
            if coal_indicator in comparison_data.index:
                coal_comparison = comparison_data.loc[coal_indicator].iloc[0] if not pd.isna(comparison_data.loc[coal_indicator].iloc[0]) else 0
            
            # 计算天然气排放变化
            gas_baseline = 0
            gas_comparison = 0
            if gas_indicator in baseline_data.index:
                gas_baseline = baseline_data.loc[gas_indicator].iloc[0] if not pd.isna(baseline_data.loc[gas_indicator].iloc[0]) else 0
            if gas_indicator in comparison_data.index:
                gas_comparison = comparison_data.loc[gas_indicator].iloc[0] if not pd.isna(comparison_data.loc[gas_indicator].iloc[0]) else 0
            
            # 计算总排放变化
            total_change = (coal_comparison + gas_comparison) - (coal_baseline + gas_baseline)
            total_changes_monthly.append(total_change)
            total_changes[ratio] += total_change
        
        change_data[ratio] = total_changes_monthly
    
    # 创建图表
    fig, ax1 = plt.subplots(figsize=(15, 8))
    
    # 设置颜色
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    # 月度总排放变化
    x = np.arange(len(month_labels))
    bar_width = 0.15
    opacity = 0.8
    
    for i, ratio in enumerate(comparison_ratios):
        if ratio in change_data:
            bars = ax1.bar(x + i * bar_width, change_data[ratio], bar_width,
                          color=colors[i], alpha=opacity, 
                          label=f"{ratio.replace('p', '%')} vs no_aluminum")
    
    ax1.set_xlabel('Month')
    ax1.set_ylabel('Total Emissions Change (tonnes CO2)')
    ax1.set_title('Monthly Total Emissions Change Due to Aluminum Production')
    ax1.set_xticks(x + bar_width * 2.5)
    ax1.set_xticklabels(month_labels)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    
    # 保存图表
    plot_file = plots_dir / f"aluminum_total_emissions_change_by_month_2050.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Aluminum total emissions change plot saved to: {plot_file}")
    
    # 创建年度总排放变化图表
    fig2, ax3 = plt.subplots(figsize=(12, 8))
    
    # 准备数据
    x = np.arange(len(comparison_ratios))
    x_labels = [f"{ratio.replace('p', '')}%" for ratio in comparison_ratios]
    
    # 绘制总排放变化
    total_emissions = [total_changes[ratio] for ratio in comparison_ratios if ratio in total_changes]
    
    bar_width = 0.6
    opacity = 0.8
    
    bars = ax3.bar(x, total_emissions, bar_width, 
                   color='#2ca02c', alpha=opacity, label='Total Emissions')
    
    # 添加数值标签
    for i, total in enumerate(total_emissions):
        if abs(total) > 1e3:
            ax3.annotate(f'{total/1e6:.1f}M',
                        xy=(i, total),
                        xytext=(0, 10 if total > 0 else -20),
                        textcoords="offset points",
                        ha='center', va='bottom' if total > 0 else 'top', 
                        fontsize=10, weight='bold')
    
    ax3.set_xlabel('Capacity Ratio (vs 55%)')
    ax3.set_ylabel('Annual Emissions Change (tonnes CO2)')
    ax3.set_title('Annual Total Emissions Change Due to Aluminum Production (2050)')
    ax3.set_xticks(x)
    ax3.set_xticklabels(x_labels)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    plt.tight_layout()
    
    # 保存总排放变化图表
    plot_file2 = plots_dir / f"aluminum_annual_emissions_change_2050.png"
    plt.savefig(plot_file2, dpi=300, bbox_inches='tight')
    logger.info(f"Annual aluminum emissions change plot saved to: {plot_file2}")
    
    # 保存数据
    data_rows = []
    for ratio in comparison_ratios:
        if ratio in change_data:
            for i, month in enumerate(month_labels):
                data_rows.append({
                    'Capacity Ratio': ratio.replace('p', '%'),
                    'Month': month,
                    'Total Emissions Change (tonnes CO2)': change_data[ratio][i]
                })
    
    data_df = pd.DataFrame(data_rows)
    data_file = plots_dir / f"aluminum_total_emissions_change_data_2050.csv"
    data_df.to_csv(data_file, index=False)
    logger.info(f"Aluminum total emissions change data saved to: {data_file}")
    
    # 保存年度汇总数据
    annual_data_rows = []
    for ratio in comparison_ratios:
        if ratio in total_changes:
            annual_data_rows.append({
                'Capacity Ratio': ratio.replace('p', '%'),
                'Annual Total Emissions Change (tonnes CO2)': total_changes[ratio]
            })
    
    annual_data_df = pd.DataFrame(annual_data_rows)
    annual_data_file = plots_dir / f"aluminum_annual_total_emissions_change_data_2050.csv"
    annual_data_df.to_csv(annual_data_file, index=False)
    logger.info(f"Annual aluminum total emissions change data saved to: {annual_data_file}")
    
    plt.close()
    plt.close(fig2)
    
    # 打印调试信息
    print(f"\n=== 电解铝导致的年度碳排放变化（相对于55%）===")
    for ratio in comparison_ratios:
        if ratio in total_changes:
            total_change = total_changes[ratio]
            print(f"{ratio.replace('p', '%')} vs 55%:")
            print(f"  总排放增加: {total_change/1e6:.3f}M tonnes CO2")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='生成电解铝导致的碳排放变化图表（相对于55%基准）')
    parser.add_argument('--output', default='results/comparison_results', help='输出目录')
    parser.add_argument('--verbose', '-v', action='store_true', help='详细输出')
    parser.add_argument('--results-dir', default='results', help='结果目录路径 (默认: results)')
    
    args = parser.parse_args()
    
    # 设置日志级别
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    
    # 定义容量比例（包括基准版本）
    capacity_ratios = ['55p', '60p', '70p', '80p', '90p', '100p']
    
    logger.info(f"开始收集容量比例 {capacity_ratios} 的2050年排放数据")
    
    # 收集所有容量比例的排放数据
    emissions_data = collect_all_emissions_data(capacity_ratios, 'emissions')
    
    if not emissions_data:
        logger.error("没有找到任何容量比例的排放数据")
        return
    
    logger.info(f"成功收集到 {len(emissions_data)} 个容量比例的排放数据")
    
    # 生成电解铝排放变化图表
    generate_aluminum_emissions_plot(emissions_data, Path(args.output))
    
    logger.info("电解铝排放变化图表生成完成！")

if __name__ == "__main__":
    main() 