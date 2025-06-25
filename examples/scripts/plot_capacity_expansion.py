# plot_capacity_expansion.py
# 容量扩展规划绘图模块

import matplotlib.pyplot as plt
import pandas as pd
import pypsa
import os

def plot_results(n, p_min_pu, save_path="examples/results"):
    """绘制电解槽用电情况结果"""
    # 获取电解槽用电数据
    smelter_p = n.links_t.p0['smelter']  # 电解槽每小时用电量
    
    # 创建图形
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # 绘制小时级数据
    ax1.plot(smelter_p, label='Hourly Usage')
    ax1.set_title('Aluminum Smelter Hourly Electricity Usage')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Power (MW)')
    ax1.legend()
    
    # 计算并绘制月度数据
    monthly_usage = smelter_p.resample('M').mean()
    ax2.bar(range(len(monthly_usage)), monthly_usage, label='Monthly Average')
    ax2.set_title('Aluminum Smelter Monthly Average Electricity Usage')
    ax2.set_xlabel('Month')
    ax2.set_ylabel('Power (MW)')
    ax2.set_xticks(range(len(monthly_usage)))
    ax2.set_xticklabels([d.strftime('%Y-%m') for d in monthly_usage.index])
    ax2.tick_params(axis='x', rotation=45)
    ax2.legend()
    
    plt.tight_layout()
    
    # 确保保存目录存在
    os.makedirs(save_path, exist_ok=True)
    
    # 保存图像
    plt.savefig(f"{save_path}/aluminum_smelter_usage_{p_min_pu}.png")
    # 显示图像
    plt.show()
    # 关闭图形，释放内存
    plt.close()

def analyze_ramp_constraints(n, config, save_path="examples/results"):
    """分析爬坡约束对电解槽运行的影响"""
    
    if 'p0' not in n.links_t:
        print("警告: 找不到链接的时间序列结果")
        return
    
    # 获取电解槽输入功率时间序列
    p = n.links_t.p0['smelter']
    
    # 计算每个时间步的功率变化
    p_diff = p.diff()
    
    # 计算最大爬坡率（上升和下降）
    max_ramp_up = p_diff[p_diff > 0].max()
    max_ramp_down = abs(p_diff[p_diff < 0].min())
    
    # 计算相对于额定功率的百分比
    p_nom = n.links.at['smelter', 'p_nom']
    max_ramp_up_pu = max_ramp_up / p_nom
    max_ramp_down_pu = max_ramp_down / p_nom
    
    print(f"\n爬坡约束分析:")
    print(f"最大上升爬坡率: {max_ramp_up:.2f} MW ({max_ramp_up_pu:.3f} p.u.)")
    print(f"最大下降爬坡率: {max_ramp_down:.2f} MW ({max_ramp_down_pu:.3f} p.u.)")
    print(f"设定的爬坡限制: {1/config['al_start_up_time']:.3f} p.u.")
    
    # 绘制功率变化图
    plt.figure(figsize=(12, 6))
    plt.plot(p_diff, 'b-', label='Power Change')
    plt.axhline(y=p_nom / config['al_start_up_time'], color='r', linestyle='--', label='Ramp Up Limit')
    plt.axhline(y=-p_nom / config['al_start_up_time'], color='r', linestyle='--', label='Ramp Down Limit')
    plt.title('Aluminum Smelter Power Changes')
    plt.xlabel('Time')
    plt.ylabel('Power Change (MW)')
    plt.legend()
    plt.grid(True)
    
    # 确保保存目录存在
    os.makedirs(save_path, exist_ok=True)
    
    plt.savefig(f'{save_path}/smelter_ramp_analysis.png', dpi=300)
    plt.show()
    plt.close()

def plot_capacity_comparison(results, config, save_path="examples/results"):
    """绘制不同最小功率比例下的容量对比图"""
    
    # 获取第一年的结果
    year_results = results[config["years"][0]]
    
    # 准备数据
    excess_rates = []
    wind_capacities = []
    solar_capacities = []
    ocgt_capacities = []
    battery_capacities = []
    h2_capacities = []
    system_costs = []
    co2_emissions = []
    
    for rate in sorted(year_results.keys()):
        result = year_results[rate]
        excess_rates.append(rate * 100)  # 转换为百分比
        
        # 发电机容量
        wind_cap = (result['generator_capacities']['onwind'] + 
                   result['generator_capacities']['offwind'])
        wind_capacities.append(wind_cap)
        solar_capacities.append(result['generator_capacities']['solar'])
        ocgt_capacities.append(result['generator_capacities']['OCGT'])
        
        # 储能容量
        battery_capacities.append(result['storage_capacities']['battery storage'])
        h2_capacities.append(result['storage_capacities']['hydrogen storage underground'])
        
        # 成本和排放
        system_costs.append(result['system_cost'])
        co2_emissions.append(result['co2_emissions'])
    
    # 创建子图
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # 发电机容量对比
    ax1.plot(excess_rates, wind_capacities, 'o-', label='Wind', linewidth=2, markersize=8)
    ax1.plot(excess_rates, solar_capacities, 's-', label='Solar', linewidth=2, markersize=8)
    ax1.plot(excess_rates, ocgt_capacities, '^-', label='OCGT', linewidth=2, markersize=8)
    ax1.set_xlabel('Minimum Power Ratio (%)')
    ax1.set_ylabel('Capacity (GW)')
    ax1.set_title('Generator Capacities vs Minimum Power Ratio')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 储能容量对比
    ax2.plot(excess_rates, battery_capacities, 'o-', label='Battery', linewidth=2, markersize=8)
    ax2.plot(excess_rates, h2_capacities, 's-', label='H2 Storage', linewidth=2, markersize=8)
    ax2.set_xlabel('Minimum Power Ratio (%)')
    ax2.set_ylabel('Capacity (GW)')
    ax2.set_title('Storage Capacities vs Minimum Power Ratio')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 系统成本
    ax3.plot(excess_rates, system_costs, 'o-', color='red', linewidth=2, markersize=8)
    ax3.set_xlabel('Minimum Power Ratio (%)')
    ax3.set_ylabel('System Cost (B€)')
    ax3.set_title('System Cost vs Minimum Power Ratio')
    ax3.grid(True, alpha=0.3)
    
    # CO2排放
    ax4.plot(excess_rates, co2_emissions, 'o-', color='green', linewidth=2, markersize=8)
    ax4.set_xlabel('Minimum Power Ratio (%)')
    ax4.set_ylabel('CO2 Emissions (Mt)')
    ax4.set_title('CO2 Emissions vs Minimum Power Ratio')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # 确保保存目录存在
    os.makedirs(save_path, exist_ok=True)
    
    plt.savefig(f'{save_path}/capacity_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

def plot_network_summary(n, p_min_pu, save_path="examples/results"):
    """绘制网络概览图"""
    
    # 创建网络图
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # 绘制网络拓扑
    n.plot(ax=ax, title=f'Network Overview (Min Power: {p_min_pu*100:.0f}%)')
    
    # 确保保存目录存在
    os.makedirs(save_path, exist_ok=True)
    
    plt.savefig(f'{save_path}/network_overview_{p_min_pu}.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

def plot_time_series_analysis(n, p_min_pu, save_path="examples/results"):
    """绘制时间序列分析图"""
    
    # 创建子图
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # 1. 电解槽用电量
    if 'p0' in n.links_t:
        smelter_p = n.links_t.p0['smelter']
        ax1.plot(smelter_p, label='Smelter Usage')
        ax1.set_title('Aluminum Smelter Power Usage')
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Power (MW)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
    
    # 2. 发电机出力
    if 'p' in n.generators_t:
        for gen in n.generators.index:
            if gen in n.generators_t.p.columns:
                ax2.plot(n.generators_t.p[gen], label=gen)
        ax2.set_title('Generator Power Output')
        ax2.set_xlabel('Time')
        ax2.set_ylabel('Power (MW)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    
    # 3. 储能状态
    if 'state_of_charge' in n.storage_units_t:
        for storage in n.storage_units.index:
            if storage in n.storage_units_t.state_of_charge.columns:
                ax3.plot(n.storage_units_t.state_of_charge[storage], label=storage)
        ax3.set_title('Storage State of Charge')
        ax3.set_xlabel('Time')
        ax3.set_ylabel('State of Charge (MWh)')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    
    # 4. 负荷曲线
    if 'p' in n.loads_t:
        for load in n.loads.index:
            if load in n.loads_t.p.columns:
                ax4.plot(n.loads_t.p[load], label=load)
        ax4.set_title('Load Profiles')
        ax4.set_xlabel('Time')
        ax4.set_ylabel('Power (MW)')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # 确保保存目录存在
    os.makedirs(save_path, exist_ok=True)
    
    plt.savefig(f'{save_path}/time_series_analysis_{p_min_pu}.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

def create_summary_plots(results, config, save_path="examples/results"):
    """创建所有汇总图表"""
    
    # 确保保存目录存在
    os.makedirs(save_path, exist_ok=True)
    
    # 创建容量对比图
    plot_capacity_comparison(results, config, save_path)
    
    print(f"所有图表已保存到 {save_path} 目录") 