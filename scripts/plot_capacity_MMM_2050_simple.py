#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
简化版MMM-2050情景成本分析图表绘制脚本
直接读取已生成的mmm_2050_M_detailed_data.csv文件来绘制图表
横轴：不同电解铝容量比例（5p-100p）
纵轴：成本节约（十亿人民币）和碳排放减少（百万吨CO2）
显示电力系统成本节约、电解铝运行成本变化、净成本节约和碳排放减少
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def load_detailed_data(csv_path):
    """
    加载详细数据CSV文件
    
    Parameters:
    -----------
    csv_path : str or Path
        CSV文件路径
        
    Returns:
    --------
    pd.DataFrame
        详细数据
    """
    try:
        df = pd.read_csv(csv_path)
        logger.info(f"成功加载数据文件: {csv_path}")
        logger.info(f"数据形状: {df.shape}")
        return df
    except Exception as e:
        logger.error(f"加载数据文件失败: {str(e)}")
        return None

def plot_mmm_2050_from_csv(csv_path, output_dir=None):
    """
    从CSV文件绘制MMM-2050情景分析图表
    
    Parameters:
    -----------
    csv_path : str or Path
        CSV文件路径
    output_dir : str or Path, optional
        输出目录，如果为None则使用CSV文件所在目录
    """
    # 加载数据
    df = load_detailed_data(csv_path)
    if df is None:
        return
    
    # 设置输出目录
    if output_dir is None:
        output_dir = Path(csv_path).parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建图表
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    # 准备数据
    x = df['Capacity_Value_Mt'].values
    power_savings = df['Power_Cost_Changes_Billion_CNY'].values
    aluminum_changes = df['Aluminum_Cost_Changes_Billion_CNY'].values
    net_savings = df['Net_Cost_Savings_Billion_CNY'].values
    emissions = df['Emissions_Changes_Million_Tonnes_CO2'].values
    
    # 创建双y轴图表
    ax2 = ax.twinx()
    
    # 设置柱子宽度
    bar_width = 150
    
    # 为电力系统成本和电解铝成本创建稍微错开的位置
    x_power = [pos - bar_width/6 for pos in x]  # 电力系统成本稍微向左
    x_aluminum = [pos + bar_width/6 for pos in x]  # 电解铝成本稍微向右
    
    # 绘制电力系统成本节约（底部）
    bars1 = ax.bar(x_power, power_savings, bar_width*0.8, color='#1f77b4', alpha=0.8, 
                   label='Power System Cost Savings')
    
    # 绘制电解铝运行成本变化（从电力系统成本顶部往上堆叠，稍微错开）
    bars2 = ax.bar(x_aluminum, aluminum_changes, bar_width*0.8, bottom=power_savings, 
                   color='#ff7f0e', alpha=0.8, label='Aluminum Operation Cost Increase')
    
    # 绘制净值曲线（黑色线，使用原始x轴位置）
    ax.plot(x, net_savings, 'k-', linewidth=3, label='Net Cost Savings', marker='o', markersize=8, zorder=20)
    
    # 找出净值最大的位置
    max_saving_index = np.argmax(net_savings)
    max_saving_value = net_savings[max_saving_index]
    max_saving_capacity = x[max_saving_index]
    
    # 用星号标出净值最大处
    ax.plot(max_saving_capacity, max_saving_value, 'r*', markersize=15, 
            label=f'Highest Net Savings: {max_saving_value:.1f}B CNY', zorder=30)
    
    # 为净值最大点添加数值标签
    ax.annotate(f'{max_saving_value:.1f}B',
                xy=(max_saving_capacity, max_saving_value),
                xytext=(0, 20),
                textcoords="offset points",
                ha='center', va='bottom', 
                fontsize=20, weight='bold', color='red',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # 绘制碳排放变化（右y轴，使用原始x轴位置）
    line1 = ax2.plot(x, emissions, linewidth=2, marker='o', 
                     markersize=6, label='Emissions Reduction', color='red')
    
    # 设置标签
    ax.set_xlabel('Aluminum Smelting Capacity (Mt)', fontsize=20)
    ax.set_ylabel('Cost Savings (Billion CNY)', fontsize=20, color='blue')
    ax2.set_ylabel('Emissions Reduction (Million Tonnes CO2)', fontsize=20, color='red')
    
    # 添加零线
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
    ax2.axhline(y=0, color='red', linestyle='--', alpha=0.5, linewidth=1)
    
    # 添加网格
    ax.grid(True, alpha=0.3, axis='y')
    
    # 设置x轴刻度和标签
    ax.set_xticks(x)
    ax.set_xticklabels([f'{cap/100:.0f}' for cap in x], fontsize=20)
    
    # 设置y轴标签
    y_ticks = ax.get_yticks()
    y_tick_labels = [f'{tick:.0f}' for tick in y_ticks]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels, fontsize=20, color='blue')
    
    # 设置右y轴标签
    y2_ticks = ax2.get_yticks()
    y2_tick_labels = [f'{tick:.0f}' for tick in y2_ticks]
    ax2.set_yticks(y2_ticks)
    ax2.set_yticklabels(y2_tick_labels, fontsize=20, color='red')
    
    # 创建图例
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor='#1f77b4', alpha=0.8, label='Power System Cost Savings'),
        plt.Rectangle((0,0),1,1, facecolor='#ff7f0e', alpha=0.8, label='Aluminum Operation Cost Increase'),
        plt.Line2D([0], [0], color='black', linewidth=3, marker='o', markersize=8, label='Net Cost Savings'),
        plt.Line2D([0], [0], marker='*', color='red', markersize=15, linestyle='', label='Highest Net Savings'),
        plt.Line2D([0], [0], color='red', linewidth=2, marker='o', markersize=6, label='Carbon Emissions Reduction')
    ]
    
    # 添加图例
    ax.legend(handles=legend_elements, loc='lower left', fontsize=20)
    
    plt.tight_layout()
    
    # 保存图表
    plot_file = output_dir / "mmm_2050_analysis.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"MMM-2050情景分析图表已保存到: {plot_file}")
    
    # 打印关键信息
    logger.info(f"Highest Net Cost Savings: {max_saving_value:.1f} Billion CNY")
    logger.info(f"Corresponding Capacity: {max_saving_capacity/100:.0f} Mt")
    logger.info(f"Corresponding Capacity Ratio: {df.iloc[max_saving_index]['Capacity_Ratio']}")
    
    return fig, ax, ax2

def main():
    """主函数"""
    # 默认CSV文件路径
    default_csv_path = "results/mmm_2050_analysis/mmm_2050_M_detailed_data.csv"
    
    # 检查文件是否存在
    csv_path = Path(default_csv_path)
    if not csv_path.exists():
        logger.error(f"CSV文件不存在: {csv_path}")
        logger.info("请确保mmm_2050_M_detailed_data.csv文件存在")
        return
    
    logger.info(f"开始绘制MMM-2050情景分析图表")
    logger.info(f"数据文件: {csv_path}")
    
    # 绘制图表
    try:
        fig, ax, ax2 = plot_mmm_2050_from_csv(csv_path)
        logger.info("图表绘制完成！")
        
        # 显示图表（如果在交互环境中）
        # plt.show()
        
    except Exception as e:
        logger.error(f"绘制图表时出错: {str(e)}")
        raise

if __name__ == "__main__":
    main()
