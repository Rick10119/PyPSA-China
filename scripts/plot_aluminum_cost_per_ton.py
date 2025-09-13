#!/usr/bin/env python3
"""
电解铝成本构成柱状图可视化程序
展示2020年和2050年的电解铝成本构成变化
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
import os

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def create_aluminum_cost_bar_chart():
    """创建电解铝成本构成柱状图"""
    
    # 数据准备
    categories = ['Raw Materials', 'Labor', 'Fixed o&m', 'Restart', 'Depreciation', 'Storage', 'Electricity']
    costs_2020 = [8451.2, 150, 1900, 0.01, 300, 1, 6250]
    costs_2050_5p = [8451.2, 159.889114, 2171.893745, 769.1493874, 1156.812339, 2.4, 4529.662849]
    costs_2050_20p = [8451.2, 159.1901429, 2984.318766, 299.1122492, 1156.812339, 15.6, 2476.284987]
    costs_2050_100p = [8451.2, 160.7041192, 7326.478149, 707.9584908, 1156.812339, 25.2, 797.0111764]
    
    # 创建DataFrame
    df = pd.DataFrame({
        '类别': categories,
        '2020成本（每吨）': costs_2020,
        '2050_5p成本（每吨）': costs_2050_5p,
        '2050_20p成本（每吨）': costs_2050_20p,
        '2050_100p成本（每吨）': costs_2050_100p
    })
    
    # 计算总成本
    total_2020 = sum(costs_2020)
    total_2050_5p = sum(costs_2050_5p)
    total_2050_20p = sum(costs_2050_20p)
    total_2050_100p = sum(costs_2050_100p)
    
    print(f"2020年总成本: {total_2020:.2f} 元/吨")
    print(f"2050_5p总成本: {total_2050_5p:.2f} 元/吨")
    print(f"2050_20p总成本: {total_2050_20p:.2f} 元/吨")
    print(f"2050_100p总成本: {total_2050_100p:.2f} 元/吨")
    
    # 创建单个子图用于堆叠柱状图
    fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    
    # 定义颜色 - 与含义相关的颜色组合
    colors = [
        '#8B4513',  # 原材料 - 棕色（代表原材料/矿物）
        '#2E8B57',  # 人工 - 深绿色（代表人力/生命）
        '#4169E1',  # 固定运维 - 蓝色（代表稳定/技术）
        '#FF4500',  # 重启 - 橙红色（代表重启/动态）
        '#696969',  # 折旧 - 灰色（代表折旧/时间）
        '#9370DB',  # 存储 - 紫色（代表存储/空间）
        '#FFD700'   # 电力 - 金色（代表电力/能源）
    ]
    
    # 准备数据 - 每个场景作为一个柱子
    scenarios = ['2024\nCurrent', 
                '2050\n5% overcapacity', 
                '2050\n36% overcapacity', 
                '2050\nNo-discommissioning']
    scenario_costs = [costs_2020, costs_2050_5p, costs_2050_20p, costs_2050_100p]
    scenario_totals = [total_2020, total_2050_5p, total_2050_20p, total_2050_100p]
    
    # 设置x轴位置
    x = np.arange(len(scenarios))
    width = 0.6
    
    # 创建堆叠柱状图
    bottom = np.zeros(len(scenarios))
    
    for i, (category, color) in enumerate(zip(categories, colors)):
        # 获取每个场景中该类别的成本
        category_costs = [costs[i] for costs in scenario_costs]
        
        # 绘制堆叠柱状图
        bars = ax.bar(x, category_costs, width, bottom=bottom, 
                     label=category, color=color, alpha=0.8, 
                     edgecolor='black', linewidth=0.5)
        
        # 添加数值标签（只对较大的值显示）
        for j, (bar, cost) in enumerate(zip(bars, category_costs)):
            if cost > 1 and cost < 8000:  # 只显示大于100的值
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2, 
                           bar.get_y() + height/2, 
                           f'{cost:.0f}', ha='center', va='center', 
                           fontsize=16, fontweight='bold', color='black')
        
        # 更新底部位置
        bottom += category_costs
    
    # 设置图表属性
    # ax.set_xlabel('Scenarios', fontsize=20, fontweight='bold')
    ax.set_ylabel('Levelized cost (CNY/tonne)', fontsize=20, fontweight='bold')
    # ax.set_title('Aluminum Cost Composition Comparison\n(Stacked Bar Chart)', 
    #             fontsize=16, fontweight='bold', pad=20)
    
    # 设置x轴标签
    ax.set_xticks(x)
    ax.set_xticklabels(scenarios, fontsize=20)
    
    # 设置y轴标签
    ax.set_yticks(np.arange(8000, 20000, 2000))
    ax.set_yticklabels([f'{int(tick)}' for tick in np.arange(8000, 20000, 2000)], fontsize=20)
    
    # 添加总成本标签在柱子顶部
    for i, total in enumerate(scenario_totals):
        ax.text(i, total + 200, f'Total: {total:.0f}', 
               ha='center', va='bottom', fontsize=20, fontweight='bold')
    
    # 添加图例 - 放在图中间，分为两行，顺序相反
    handles, labels = ax.get_legend_handles_labels()
    # 反转顺序
    handles = handles[::-1]
    labels = labels[::-1]
    ax.legend(handles, labels, bbox_to_anchor=(0.5, -0.1), loc='upper center', 
              ncol=4, fontsize=20, frameon=False)
    
    # 添加网格
    ax.grid(True, alpha=0.3, axis='y')
    
    # 设置y轴范围
    ax.set_ylim(8000, 20000)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    output_path = 'results/aluminum_cost_composition_2020_2050_stacked_bar.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Image saved to: {output_path}")
    
    # 显示图片
    plt.show()
    
    return df

def create_detailed_comparison_table():
    """创建详细的成本对比表格"""
    
    # 数据准备
    categories = ['Raw Materials', 'Labor', 'Fixed o&m', 'Restart', 'Depreciation', 'Storage', 'Electricity']
    costs_2020 = [8451.2, 150, 1900, 0, 300, 1, 6250]
    costs_2050_5p = [8451.2, 159.889114, 2171.893745, 769.1493874, 1156.812339, 1, 4529.662849]
    costs_2050_20p = [8451.2, 159.1901429, 2984.318766, 299.1122492, 1156.812339, 15.6, 2476.284987]
    costs_2050_100p = [8451.2, 160.7041192, 7326.478149, 707.9584908, 1156.812339, 15.6, 797.0111764]
    
    # 创建DataFrame
    df = pd.DataFrame({
        '类别': categories,
        '2020成本（元/吨）': costs_2020,
        '2050_5p成本（元/吨）': costs_2050_5p,
        '2050_20p成本（元/吨）': costs_2050_20p,
        '2050_100p成本（元/吨）': costs_2050_100p
    })
    
    # 计算总成本
    total_2020 = df['2020成本（元/吨）'].sum()
    total_2050_5p = df['2050_5p成本（元/吨）'].sum()
    total_2050_20p = df['2050_20p成本（元/吨）'].sum()
    total_2050_100p = df['2050_100p成本（元/吨）'].sum()
    
    print("\n=== 电解铝成本构成详细对比 ===")
    print(df.to_string(index=False, float_format='%.2f'))
    print(f"\n总成本变化:")
    print(f"2020年总成本: {total_2020:.2f} 元/吨")
    print(f"2050_5p总成本: {total_2050_5p:.2f} 元/吨 (变化: {total_2050_5p - total_2020:.2f} 元/吨)")
    print(f"2050_20p总成本: {total_2050_20p:.2f} 元/吨 (变化: {total_2050_20p - total_2020:.2f} 元/吨)")
    print(f"2050_100p总成本: {total_2050_100p:.2f} 元/吨 (变化: {total_2050_100p - total_2020:.2f} 元/吨)")
    
    return df

def main():
    """主函数"""
    print("Creating aluminum cost composition bar chart...")
    
    # 创建柱状图
    df = create_aluminum_cost_bar_chart()
    
    # 创建详细对比表格
    create_detailed_comparison_table()
    
    print("\nVisualization completed!")

if __name__ == "__main__":
    main()
