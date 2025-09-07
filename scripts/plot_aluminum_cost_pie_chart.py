#!/usr/bin/env python3
"""
电解铝成本构成饼状图可视化程序
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

def create_aluminum_cost_pie_chart():
    """创建电解铝成本构成饼状图"""
    
    # 数据准备
    categories = ['Raw Materials', 'Labor', 'Fixed o&m', 'Restart', 'Depreciation', 'Storage', 'Electricity']
    costs_2020 = [8451.2, 150, 1900, 0.01, 300, 1, 6250]
    costs_2050_non_flex = [8451.2, 159.889114, 2171.893745, 769.1493874, 1156.812339, 1, 4529.662849]
    costs_2050_20p = [8451.2, 159.1901429, 2984.318766, 299.1122492, 1156.812339, 15.6, 2476.284987]
    costs_2050_100p = [8451.2, 160.7041192, 7326.478149, 707.9584908, 1156.812339, 15.6, 797.0111764]
    
    # 创建DataFrame
    df = pd.DataFrame({
        '类别': categories,
        '2020成本（每吨）': costs_2020,
        '2050_non_flex成本（每吨）': costs_2050_non_flex,
        '2050_20p成本（每吨）': costs_2050_20p,
        '2050_100p成本（每吨）': costs_2050_100p
    })
    
    # 计算总成本
    total_2020 = sum(costs_2020)
    total_2050_non_flex = sum(costs_2050_non_flex)
    total_2050_20p = sum(costs_2050_20p)
    total_2050_100p = sum(costs_2050_100p)
    
    print(f"2020年总成本: {total_2020:.2f} 元/吨")
    print(f"2050_non_flex总成本: {total_2050_non_flex:.2f} 元/吨")
    print(f"2050_20p总成本: {total_2050_20p:.2f} 元/吨")
    print(f"2050_100p总成本: {total_2050_100p:.2f} 元/吨")
    
    # 创建子图 - 2x2布局
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 12))
    
    # 定义颜色
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7', '#DDA0DD', '#98D8C8']
    
    # 2020年饼状图
    # 过滤掉成本为0的项目
    non_zero_2020 = [(cat, cost) for cat, cost in zip(categories, costs_2020) if cost > 0]
    labels_2020 = [item[0] for item in non_zero_2020]
    sizes_2020 = [item[1] for item in non_zero_2020]
    colors_2020 = colors[:len(non_zero_2020)]
    
    wedges1, texts1, autotexts1 = ax1.pie(sizes_2020, labels=labels_2020, colors=colors_2020, 
                                         autopct='%1.1f%%', startangle=90, textprops={'fontsize': 12})
    ax1.set_title('2020 Aluminum Cost Composition\n(Total Cost: {:.0f} CNY/ton)'.format(total_2020), 
                  fontsize=14, fontweight='bold', pad=15)
    
    # 2050_non_flex饼状图
    wedges2, texts2, autotexts2 = ax2.pie(costs_2050_non_flex, labels=categories, colors=colors, 
                                         autopct='%1.1f%%', startangle=90, textprops={'fontsize': 12})
    ax2.set_title('2050_non_flex Aluminum Cost Composition\n(Total Cost: {:.0f} CNY/ton)'.format(total_2050_non_flex), 
                  fontsize=14, fontweight='bold', pad=15)
    
    # 2050_20p饼状图
    wedges3, texts3, autotexts3 = ax3.pie(costs_2050_20p, labels=categories, colors=colors, 
                                         autopct='%1.1f%%', startangle=90, textprops={'fontsize': 12})
    ax3.set_title('2050_20p Aluminum Cost Composition\n(Total Cost: {:.0f} CNY/ton)'.format(total_2050_20p), 
                  fontsize=14, fontweight='bold', pad=15)
    
    # 2050_100p饼状图
    wedges4, texts4, autotexts4 = ax4.pie(costs_2050_100p, labels=categories, colors=colors, 
                                         autopct='%1.1f%%', startangle=90, textprops={'fontsize': 12})
    ax4.set_title('2050_100p Aluminum Cost Composition\n(Total Cost: {:.0f} CNY/ton)'.format(total_2050_100p), 
                  fontsize=14, fontweight='bold', pad=15)
    
    # 调整布局
    plt.tight_layout()
    
    # 添加总标题
    fig.suptitle('Aluminum Cost Composition Comparison (2020 vs 2050 Scenarios)', fontsize=18, fontweight='bold', y=0.95)
    
    # 保存图片
    output_path = 'results/aluminum_cost_composition_2020_2050.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Image saved to: {output_path}")
    
    # 显示图片
    # plt.show()
    
    return df

def create_detailed_comparison_table():
    """创建详细的成本对比表格"""
    
    # 数据准备
    categories = ['Raw Materials', 'Labor', 'Fixed o&m', 'Restart', 'Depreciation', 'Storage', 'Electricity']
    costs_2020 = [8451.2, 150, 1900, 0, 300, 1, 6250]
    costs_2050_non_flex = [8451.2, 159.889114, 2171.893745, 769.1493874, 1156.812339, 1, 4529.662849]
    costs_2050_20p = [8451.2, 159.1901429, 2984.318766, 299.1122492, 1156.812339, 15.6, 2476.284987]
    costs_2050_100p = [8451.2, 160.7041192, 7326.478149, 707.9584908, 1156.812339, 15.6, 797.0111764]
    
    # 创建DataFrame
    df = pd.DataFrame({
        '类别': categories,
        '2020成本（元/吨）': costs_2020,
        '2050_non_flex成本（元/吨）': costs_2050_non_flex,
        '2050_20p成本（元/吨）': costs_2050_20p,
        '2050_100p成本（元/吨）': costs_2050_100p
    })
    
    # 计算总成本
    total_2020 = df['2020成本（元/吨）'].sum()
    total_2050_non_flex = df['2050_non_flex成本（元/吨）'].sum()
    total_2050_20p = df['2050_20p成本（元/吨）'].sum()
    total_2050_100p = df['2050_100p成本（元/吨）'].sum()
    
    print("\n=== 电解铝成本构成详细对比 ===")
    print(df.to_string(index=False, float_format='%.2f'))
    print(f"\n总成本变化:")
    print(f"2020年总成本: {total_2020:.2f} 元/吨")
    print(f"2050_non_flex总成本: {total_2050_non_flex:.2f} 元/吨 (变化: {total_2050_non_flex - total_2020:.2f} 元/吨)")
    print(f"2050_20p总成本: {total_2050_20p:.2f} 元/吨 (变化: {total_2050_20p - total_2020:.2f} 元/吨)")
    print(f"2050_100p总成本: {total_2050_100p:.2f} 元/吨 (变化: {total_2050_100p - total_2020:.2f} 元/吨)")
    
    return df

def main():
    """主函数"""
    print("Creating aluminum cost composition pie chart...")
    
    # 创建饼状图
    df = create_aluminum_cost_pie_chart()
    
    # 创建详细对比表格
    create_detailed_comparison_table()
    
    print("\nVisualization completed!")

if __name__ == "__main__":
    main()
