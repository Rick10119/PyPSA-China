#!/usr/bin/env python3
"""
生成分省铝需求和装机数据（单位：吨/小时）

输出文件：
1. aluminum_demand_by_province.csv - 分省需求（按年份和情景），单位：吨/小时
2. aluminum_capacity_by_province.csv - 分省装机，单位：吨/小时
"""

import pandas as pd
import json
import numpy as np
from pathlib import Path

# 文件路径
BASE_DIR = Path(__file__).parent.parent
DEMAND_JSON = BASE_DIR / "data" / "aluminum_demand" / "aluminum_demand_all_scenarios.json"
CAPACITY_CSV = BASE_DIR / "data" / "p_nom" / "al_smelter_p_max.csv"
OUTPUT_DIR = BASE_DIR / "data" / "aluminum_demand"

# 转换常数
HOURS_PER_YEAR = 8760  # 一年的小时数
TONS_PER_10KT = 10000  # 1万吨 = 10000吨


def load_demand_data():
    """加载需求数据"""
    with open(DEMAND_JSON, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data['primary_aluminum_demand']


def load_capacity_data():
    """加载装机容量数据"""
    df = pd.read_csv(CAPACITY_CSV)
    df = df.set_index('Province')
    # 过滤掉产量为0或很小的省份
    df = df[df['p_nom'] > 0.01]
    return df


def calculate_production_ratio(capacity_df):
    """计算各省份的生产比例"""
    total_capacity = capacity_df['p_nom'].sum()
    production_ratio = capacity_df['p_nom'] / total_capacity
    return production_ratio


def convert_demand_to_tons_per_hour(demand_10kt):
    """
    将需求从10kt转换为吨/小时
    
    参数:
        demand_10kt: 需求（10kt，万吨）
    
    返回:
        需求（吨/小时）
    """
    demand_tons = demand_10kt * TONS_PER_10KT  # 转换为吨
    demand_tons_per_hour = demand_tons / HOURS_PER_YEAR  # 转换为吨/小时
    return demand_tons_per_hour


def convert_capacity_to_tons_per_hour(capacity_10kt_per_year):
    """
    将装机容量从10kt/year转换为吨/小时
    
    参数:
        capacity_10kt_per_year: 年产量（10kt/year，万吨/年）
    
    返回:
        产能（吨/小时）
    """
    capacity_tons_per_year = capacity_10kt_per_year * TONS_PER_10KT  # 转换为吨/年
    capacity_tons_per_hour = capacity_tons_per_year / HOURS_PER_YEAR  # 转换为吨/小时
    return capacity_tons_per_hour


def generate_demand_by_province():
    """生成分省需求数据"""
    # 加载数据
    demand_data = load_demand_data()
    capacity_df = load_capacity_data()
    production_ratio = calculate_production_ratio(capacity_df)
    
    # 准备输出数据
    results = []
    
    # 遍历所有情景和年份
    for scenario in ['low', 'mid', 'high']:
        if scenario not in demand_data:
            continue
        
        # 添加2025年数据（所有情景都是4000万吨）
        demand_data[scenario]['2025'] = 4000.0
        
        for year, demand_10kt in demand_data[scenario].items():
            # 转换为吨/小时
            national_demand_tons_per_hour = convert_demand_to_tons_per_hour(demand_10kt)
            
            # 按省份分配
            for province in production_ratio.index:
                province_demand = national_demand_tons_per_hour * production_ratio[province]
                results.append({
                    'Province': province,
                    'Year': year,
                    'Scenario': scenario,
                    'Demand_ton_per_h': province_demand
                })
    
    # 转换为DataFrame
    df = pd.DataFrame(results)
    
    # 重新排列列顺序
    df = df[['Province', 'Year', 'Scenario', 'Demand_ton_per_h']]
    
    return df


def generate_capacity_by_province():
    """生成分省装机数据"""
    # 加载数据
    capacity_df = load_capacity_data()
    
    # 转换为吨/小时
    capacity_tons_per_hour = convert_capacity_to_tons_per_hour(capacity_df['p_nom'])
    
    # 创建输出DataFrame，重置索引避免冲突
    df = pd.DataFrame({
        'Province': capacity_df.index,
        'Capacity_10kt_per_year': capacity_df['p_nom'].values,
        'Capacity_ton_per_h': capacity_tons_per_hour.values
    })
    
    # 按省份排序
    df = df.sort_values('Province').reset_index(drop=True)
    
    return df


def main():
    """主函数"""
    # 创建输出目录
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("正在生成分省铝需求数据...")
    demand_df = generate_demand_by_province()
    demand_output = OUTPUT_DIR / "aluminum_demand_by_province.csv"
    demand_df.to_csv(demand_output, index=False, encoding='utf-8-sig')
    print(f"✓ 需求数据已保存到: {demand_output}")
    print(f"  包含 {len(demand_df)} 条记录")
    print(f"  省份数: {demand_df['Province'].nunique()}")
    print(f"  年份: {sorted(demand_df['Year'].unique())}")
    print(f"  情景: {sorted(demand_df['Scenario'].unique())}")
    
    print("\n正在生成分省装机数据...")
    capacity_df = generate_capacity_by_province()
    capacity_output = OUTPUT_DIR / "aluminum_capacity_by_province.csv"
    capacity_df.to_csv(capacity_output, index=False, encoding='utf-8-sig')
    print(f"✓ 装机数据已保存到: {capacity_output}")
    print(f"  包含 {len(capacity_df)} 个省份")
    print(f"  总装机: {capacity_df['Capacity_ton_per_h'].sum():.6f} 吨/小时")
    
    # 显示一些统计信息
    print("\n=== 需求数据预览 ===")
    print(demand_df.head(10))
    
    print("\n=== 装机数据预览 ===")
    print(capacity_df.head(10))
    
    # 按年份汇总需求
    print("\n=== 全国需求汇总（吨/小时）===")
    demand_summary = demand_df.groupby(['Year', 'Scenario'])['Demand_ton_per_h'].sum()
    print(demand_summary)
    
    print("\n完成！")


if __name__ == "__main__":
    main()

