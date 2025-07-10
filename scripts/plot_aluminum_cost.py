# SPDX-FileCopyrightText: : 2022 The PyPSA-China Authors
#
# SPDX-License-Identifier: MIT

"""
对比两个不同版本的结果summary的铝冶炼电费成本差异
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
VERSION1 = "0701.1H.7"  # 基准版本
VERSION2 = "0701.1H.8"  # 对比版本

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

def collect_aluminum_statistics_data(name1, name2):
    """
    从所有年份目录中收集铝冶炼统计数据
    
    Parameters:
    -----------
    name1 : str
        第一个版本名称（基准版本，使用平均边际电价计算）
    name2 : str
        第二个版本名称（对比版本，读取铝冶炼数据）
        
    Returns:
    --------
    dict
        按年份组织的铝冶炼统计数据
    """
    aluminum_data = {}
    
    # 定义年份列表
    years = [2020, 2025, 2030, 2035, 2040, 2045, 2050]
    
    for year in years:
        # 构建目录路径
        dir_pattern = f"postnetwork-ll-current+Neighbor-linear2050-{year}"
        
        # 检查对比版本的目录是否存在
        version2_dir = Path(f"results/version-{name2}/summary/postnetworks/positive/{dir_pattern}")
        
        # 只处理对比版本目录存在的情况
        if version2_dir.exists():
            # 读取对比版本的aluminum_statistics文件
            aluminum_file2 = version2_dir / "aluminum_statistics.csv"
            
            # 只处理对比版本文件存在的情况
            if aluminum_file2.exists():
                try:
                    # 读取对比版本数据
                    df2 = pd.read_csv(aluminum_file2, index_col=0)
                    
                    # 确保数据类型为数值型
                    for col in df2.columns:
                        df2[col] = pd.to_numeric(df2[col], errors='coerce')
                    
                    print(f"\n==== {year}年 aluminum_statistics index ====")
                    print(df2.index.tolist())
                    if 'total_electricity_cost' in df2.index:
                        print(f"{year}年 total_electricity_cost 原始值: {df2.loc['total_electricity_cost'].values}")
                    else:
                        print(f"{year}年 没有total_electricity_cost，所有index: {df2.index.tolist()}")
                    
                    # 获取对比版本的铝冶炼用电量
                    if 'total_electricity_consumption' in df2.index:
                        aluminum_consumption = df2.loc['total_electricity_consumption'].iloc[0] if not pd.isna(df2.loc['total_electricity_consumption'].iloc[0]) else 0
                    else:
                        aluminum_consumption = 0
                    
                    # 如果对比版本的铝冶炼用电量为0，跳过该年份
                    if aluminum_consumption == 0:
                        logger.info(f"跳过 {year} 年：对比版本铝冶炼用电量为0")
                        continue
                    
                    # 获取基准版本的平均边际电价
                    version1_dir = Path(f"results/version-{name1}/summary/postnetworks/positive/{dir_pattern}")
                    if version1_dir.exists():
                        # 尝试读取prices.csv文件获取平均边际电价
                        prices_file1 = version1_dir / "prices.csv"
                        if prices_file1.exists():
                            try:
                                df_prices1 = pd.read_csv(prices_file1, index_col=0)
                                # 获取AC电价的平均值
                                if 'AC' in df_prices1.index:
                                    avg_price_v1 = df_prices1.loc['AC'].iloc[0] if not pd.isna(df_prices1.loc['AC'].iloc[0]) else 0
                                else:
                                    avg_price_v1 = 0
                            except Exception as e:
                                logger.warning(f"读取 {year} 年基准版本电价数据时出错: {str(e)}")
                                avg_price_v1 = 0
                        else:
                            logger.warning(f"在 {year} 年基准版本目录中未找到 prices.csv 文件")
                            avg_price_v1 = 0
                    else:
                        logger.warning(f"基准版本 {year} 年目录不存在")
                        avg_price_v1 = 0
                    
                    # 获取对比版本的平均电价
                    if 'average_electricity_price' in df2.index:
                        avg_price_v2 = df2.loc['average_electricity_price'].iloc[0] if not pd.isna(df2.loc['average_electricity_price'].iloc[0]) else 0
                    else:
                        avg_price_v2 = 0
                    
                    # 计算基准版本的电费成本（用电量 × 平均边际电价）
                    cost_v1 = aluminum_consumption * avg_price_v1
                    cost_v2 = df2.loc['total_electricity_cost'].iloc[0] if 'total_electricity_cost' in df2.index and not pd.isna(df2.loc['total_electricity_cost'].iloc[0]) else 0
                    
                    # 创建基准版本的模拟数据
                    df1_data = {
                        'total_electricity_cost': cost_v1,
                        'total_electricity_consumption': aluminum_consumption,
                        'average_electricity_price': avg_price_v1
                    }
                    df1 = pd.DataFrame([df1_data])
                    
                    aluminum_data[year] = {
                        name1: df1,
                        name2: df2
                    }
                    logger.info(f"成功加载 {year} 年的铝冶炼统计数据")
                    logger.info(f"  基准版本：用电量={aluminum_consumption/1e6:.2f}TWh, 平均电价={avg_price_v1:.2f}EUR/MWh, 电费成本={cost_v1/1e9:.3f}BEUR")
                    logger.info(f"  对比版本：用电量={aluminum_consumption/1e6:.2f}TWh, 平均电价={avg_price_v2:.2f}EUR/MWh, 电费成本={cost_v2/1e9:.3f}BEUR")
                    
                except Exception as e:
                    logger.warning(f"加载 {year} 年的铝冶炼统计数据时出错: {str(e)}")
            else:
                logger.info(f"跳过 {year} 年：对比版本缺少aluminum_statistics.csv文件")
        else:
            logger.info(f"跳过 {year} 年：对比版本目录不存在")
    
    return aluminum_data

def generate_detailed_aluminum_cost_plot(aluminum_data, name1, name2, plots_dir):
    """
    生成详细的铝冶炼电费成本对比图表
    
    Parameters:
    -----------
    aluminum_data : dict
        按年份组织的铝冶炼统计数据
    name1 : str
        第一个版本名称
    name2 : str
        第二个版本名称
    plots_dir : Path
        图表输出目录
    """
    # 提取年份
    years = sorted(aluminum_data.keys())
    
    if not years:
        logger.warning("没有找到任何有效的铝冶炼统计数据")
        return
    
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    # 准备数据
    total_electricity_costs_v1 = []
    total_electricity_costs_v2 = []
    average_prices_v1 = []
    average_prices_v2 = []
    electricity_consumptions_v1 = []
    electricity_consumptions_v2 = []
    valid_years = []
    
    for year in years:
        if name1 in aluminum_data[year] and name2 in aluminum_data[year]:
            df1 = aluminum_data[year][name1]
            df2 = aluminum_data[year][name2]
            
            print(f"\n==== {year}年 详细对比 ====")
            print(f"df1 index: {df1.index.tolist()}")
            print(f"df2 index: {df2.index.tolist()}")
            if 'total_electricity_cost' in df2.index:
                print(f"cost_v2: {df2.loc['total_electricity_cost'].values}")
            else:
                print(f"cost_v2: 没有该行")
            
            # 获取铝冶炼电费成本数据
            # df1是单行DataFrame，df2是index为统计项的DataFrame
            if 'total_electricity_cost' in df2.index:
                cost_v2 = df2.loc['total_electricity_cost'].values[0] if not pd.isna(df2.loc['total_electricity_cost'].values[0]) else 0
            else:
                cost_v2 = 0
            if 'total_electricity_cost' in df1.columns:
                cost_v1 = df1['total_electricity_cost'].iloc[0] if not pd.isna(df1['total_electricity_cost'].iloc[0]) else 0
            else:
                cost_v1 = 0

            # 获取平均电价数据
            if 'average_electricity_price' in df2.index:
                price_v2 = df2.loc['average_electricity_price'].values[0] if not pd.isna(df2.loc['average_electricity_price'].values[0]) else 0
            else:
                price_v2 = 0
            if 'average_electricity_price' in df1.columns:
                price_v1 = df1['average_electricity_price'].iloc[0] if not pd.isna(df1['average_electricity_price'].iloc[0]) else 0
            else:
                price_v1 = 0

            # 获取用电量数据
            if 'total_electricity_consumption' in df2.index:
                consumption_v2 = df2.loc['total_electricity_consumption'].values[0] if not pd.isna(df2.loc['total_electricity_consumption'].values[0]) else 0
            else:
                consumption_v2 = 0
            if 'total_electricity_consumption' in df1.columns:
                consumption_v1 = df1['total_electricity_consumption'].iloc[0] if not pd.isna(df1['total_electricity_consumption'].iloc[0]) else 0
            else:
                consumption_v1 = 0

            # 检查对比版本（name2）的数值是否为0，如果是则跳过
            if cost_v2 == 0:
                logger.info(f"跳过 {year} 年：对比版本电费成本为0")
                continue
            
            # 添加到有效数据列表
            valid_years.append(year)
            total_electricity_costs_v1.append(cost_v1 * EUR_TO_CNY)  # 转换为人民币
            total_electricity_costs_v2.append(cost_v2 * EUR_TO_CNY)
            average_prices_v1.append(price_v1)
            average_prices_v2.append(price_v2)
            electricity_consumptions_v1.append(consumption_v1)
            electricity_consumptions_v2.append(consumption_v2)
    
    if not valid_years:
        logger.warning("没有找到有效的对比数据")
        return
    
    logger.info(f"将生成 {len(valid_years)} 年的对比图表：{valid_years}")
    
    # 创建子图
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    x = np.arange(len(valid_years))
    width = 0.35
    
    # 1. 铝冶炼电费成本对比
    bars1 = ax1.bar(x - width/2, total_electricity_costs_v1, width, label=name1, alpha=0.8, color='skyblue')
    bars2 = ax1.bar(x + width/2, total_electricity_costs_v2, width, label=name2, alpha=0.8, color='lightcoral')
    
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Total Electricity Cost (Billion CNY)')
    ax1.set_title('Aluminum Smelting Electricity Cost Comparison')
    ax1.set_xticks(x)
    ax1.set_xticklabels(valid_years, rotation=45)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 添加数值标签
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax1.annotate(f'{height/1e9:.1f}B',
                            xy=(bar.get_x() + bar.get_width()/2, height),
                            xytext=(0, 3),
                            textcoords="offset points",
                            ha='center', va='bottom', fontsize=8)
    
    # 2. 平均电价对比
    ax2.plot(x, average_prices_v1, 'o-', label=name1, linewidth=2, markersize=8, color='blue')
    ax2.plot(x, average_prices_v2, 's-', label=name2, linewidth=2, markersize=8, color='red')
    
    ax2.set_xlabel('Year')
    ax2.set_ylabel('Average Electricity Price (EUR/MWh)')
    ax2.set_title('Average Electricity Price Comparison')
    ax2.set_xticks(x)
    ax2.set_xticklabels(valid_years, rotation=45)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 添加数值标签
    for i, (price1, price2) in enumerate(zip(average_prices_v1, average_prices_v2)):
        if price1 > 0:
            ax2.annotate(f'{price1:.1f}',
                        xy=(i, price1),
                        xytext=(0, 5),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=8)
        if price2 > 0:
            ax2.annotate(f'{price2:.1f}',
                        xy=(i, price2),
                        xytext=(0, -15),
                        textcoords="offset points",
                        ha='center', va='top', fontsize=8)
    
    # 3. 用电量对比
    bars3 = ax3.bar(x - width/2, electricity_consumptions_v1, width, label=name1, alpha=0.8, color='lightgreen')
    bars4 = ax3.bar(x + width/2, electricity_consumptions_v2, width, label=name2, alpha=0.8, color='orange')
    
    ax3.set_xlabel('Year')
    ax3.set_ylabel('Electricity Consumption (TWh)')
    ax3.set_title('Aluminum Smelting Electricity Consumption')
    ax3.set_xticks(x)
    ax3.set_xticklabels(valid_years, rotation=45)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 添加数值标签
    for bars in [bars3, bars4]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax3.annotate(f'{height/1e6:.1f}T',
                            xy=(bar.get_x() + bar.get_width()/2, height),
                            xytext=(0, 3),
                            textcoords="offset points",
                            ha='center', va='bottom', fontsize=8)
    
    # 4. 成本变化量（节约为正，增加为负）
    cost_changes = [cost1 - cost2 for cost1, cost2 in zip(total_electricity_costs_v1, total_electricity_costs_v2)]
    colors = ['green' if change > 0 else 'red' for change in cost_changes]
    
    bars5 = ax4.bar(x, cost_changes, color=colors, alpha=0.7)
    ax4.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    
    ax4.set_xlabel('Year')
    ax4.set_ylabel('Cost Change (Billion CNY)')
    ax4.set_title(f'Electricity Cost Change ({name1} - {name2})')
    ax4.set_xticks(x)
    ax4.set_xticklabels(valid_years, rotation=45)
    ax4.grid(True, alpha=0.3)
    
    # 添加数值标签
    for i, (bar, change) in enumerate(zip(bars5, cost_changes)):
        if abs(change) > 1e6:  # 只显示大于1M的变化
            ax4.annotate(f'{change/1e9:.1f}B',
                        xy=(bar.get_x() + bar.get_width()/2, change),
                        xytext=(0, 10 if change > 0 else -20),
                        textcoords="offset points",
                        ha='center', va='bottom' if change > 0 else 'top', 
                        fontsize=8, weight='bold',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    plt.tight_layout()
    
    # 保存图表
    plot_file = plots_dir / f"detailed_aluminum_electricity_cost_{name1}_vs_{name2}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Detailed aluminum electricity cost plot saved to: {plot_file}")
    
    # 保存数据
    data_rows = []
    for i, year in enumerate(valid_years):
        data_rows.append({
            'Year': year,
            f'{name1}_Total_Cost_CNY': total_electricity_costs_v1[i],
            f'{name2}_Total_Cost_CNY': total_electricity_costs_v2[i],
            f'{name1}_Average_Price_EUR': average_prices_v1[i],
            f'{name2}_Average_Price_EUR': average_prices_v2[i],
            f'{name1}_Consumption_TWh': electricity_consumptions_v1[i] / 1e6,
            f'{name2}_Consumption_TWh': electricity_consumptions_v2[i] / 1e6,
            'Cost_Change_CNY': cost_changes[i]
        })
    
    data_df = pd.DataFrame(data_rows)
    data_file = plots_dir / f"detailed_aluminum_electricity_cost_data_{name1}_vs_{name2}.csv"
    data_df.to_csv(data_file, index=False)
    logger.info(f"Detailed aluminum electricity cost data saved to: {data_file}")
    
    plt.close()

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='对比两个不同版本的铝冶炼电费成本')
    parser.add_argument('--name1', help='第一个版本显示名称 (默认使用VERSION1)')
    parser.add_argument('--name2', help='第二个版本显示名称 (默认使用VERSION2)')
    parser.add_argument('--output', default='results/aluminum_comparison', help='输出目录')
    parser.add_argument('--verbose', '-v', action='store_true', help='详细输出')
    parser.add_argument('--results-dir', default='results', help='结果目录路径 (默认: results)')
    
    args = parser.parse_args()
    
    # 设置日志级别
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    
    # 使用脚本中设置的版本号
    version1 = VERSION1
    version2 = VERSION2
    
    # 设置版本名称
    name1 = args.name1 if args.name1 else version1
    name2 = args.name2 if args.name2 else version2
    
    logger.info(f"开始对比版本 {name1} 和 {name2} 的铝冶炼电费成本")
    
    # 创建输出目录
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 收集铝冶炼统计数据并生成对比图表
    aluminum_data = collect_aluminum_statistics_data(name1, name2)
    
    if aluminum_data:
        # 创建图表目录
        plots_dir = output_path / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        
        # 生成铝冶炼电费成本详细对比图表
        generate_detailed_aluminum_cost_plot(aluminum_data, name1, name2, plots_dir)
        logger.info("铝冶炼电费成本对比完成！")
    else:
        logger.warning("没有找到铝冶炼统计数据")
        logger.info("请确保运行了包含铝冶炼功能的版本，并生成了aluminum_statistics.csv文件")

if __name__ == "__main__":
    main() 