# plot_capacity_expansion.py
# 容量扩展规划绘图模块

import matplotlib.pyplot as plt
import pandas as pd
import pypsa
import os

"""
容量扩展规划绘图模块

本模块包含所有与电解铝迭代优化算法相关的绘图函数：

基础绘图函数：
- plot_aluminum_usage(): 绘制电解铝用能模式
- plot_nodal_prices(): 绘制节点电价时间序列
- save_and_show_plot(): 统一处理图片保存和显示

迭代优化绘图函数：
- plot_iteration_results(): 绘制每次迭代的结果
- plot_iterative_results(): 绘制最终迭代结果

传统绘图函数：
- plot_results(): 绘制电解槽用电情况结果
- analyze_ramp_constraints(): 分析爬坡约束
- plot_capacity_comparison(): 绘制容量对比图
- plot_network_summary(): 绘制网络概览图
- plot_time_series_analysis(): 绘制时间序列分析图
- create_summary_plots(): 创建所有汇总图表

模块化设计：
- 所有绘图函数集中在一个模块中，便于维护和复用
- 支持不同的绘图需求（迭代过程、结果分析、网络展示等）
- 统一的图片保存和显示机制
"""

def plot_results(n, p_min_pu, save_path="examples/results"):
    """绘制电解槽用电情况结果（仅2月数据）"""
    # 获取电解槽用电数据
    smelter_p = n.links_t.p0['smelter']  # 电解槽每小时用电量
    
    # 筛选4月数据
    month_data = smelter_p[smelter_p.index.month == 2]
    
    if len(month_data) == 0:
        print("警告: 没有找到4月的数据")
        return
    
    # 创建图形
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # 绘制4月小时级数据
    ax1.plot(month_data, label='April Hourly Usage', color='blue', linewidth=1)
    ax1.set_title('Aluminum Smelter April Hourly Electricity Usage')
    ax1.set_xlabel('Time (April)')
    ax1.set_ylabel('Power (MW)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 计算并绘制4月日平均数据
    daily_usage = month_data.resample('D').mean()
    ax2.bar(range(len(daily_usage)), daily_usage, label='April Daily Average', color='orange', alpha=0.7)
    ax2.set_title('Aluminum Smelter April Daily Average Electricity Usage')
    ax2.set_xlabel('Day of April')
    ax2.set_ylabel('Power (MW)')
    ax2.set_xticks(range(len(daily_usage)))
    ax2.set_xticklabels([d.strftime('%m-%d') for d in daily_usage.index])
    ax2.tick_params(axis='x', rotation=45)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # 确保保存目录存在
    os.makedirs(save_path, exist_ok=True)
    
    # 保存图像
    plt.savefig(f"{save_path}/aluminum_smelter_usage_april_{p_min_pu}.png")
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
    
    # plt.savefig(f'{save_path}/time_series_analysis_{p_min_pu}.png', dpi=300, bbox_inches='tight')
    # plt.show()
    # plt.close()

def create_summary_plots(results, config, save_path="examples/results"):
    """创建所有汇总图表"""
    
    # 确保保存目录存在
    os.makedirs(save_path, exist_ok=True)
    
    # 创建容量对比图
    plot_capacity_comparison(results, config, save_path)
    
    print(f"所有图表已保存到 {save_path} 目录")

def plot_aluminum_usage(ax, aluminum_usage, ts, title_suffix=""):
    """
    绘制电解铝用能模式
    ax: matplotlib轴对象
    aluminum_usage: 电解铝用能数据
    ts: 时间序列数据
    title_suffix: 标题后缀
    """
    if 'smelter' in aluminum_usage.columns:
        ax.plot(aluminum_usage.index, aluminum_usage['smelter'], 
                label='Aluminum Smelter Power', linewidth=2, color='red')
        print(f"电解铝用能范围: {aluminum_usage['smelter'].min():.2f} - {aluminum_usage['smelter'].max():.2f} MW")
    else:
        print("警告：电解铝用能数据中没有'smelter'列")
        print(f"可用列: {list(aluminum_usage.columns)}")
    
    ax.plot(ts.index, ts.aluminum, label='Aluminum Demand', linewidth=2, color='blue', linestyle='--')
    ax.set_title(f'Aluminum Smelter Power Pattern{title_suffix}', fontsize=14)
    ax.set_ylabel('Power (MW)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

def plot_nodal_prices(ax, network, title_suffix=""):
    """
    绘制节点电价
    ax: matplotlib轴对象
    network: PyPSA网络对象
    title_suffix: 标题后缀
    """
    # 检查是否有边际价格数据
    if hasattr(network, 'buses_t') and hasattr(network.buses_t, 'marginal_price'):
        print(f"找到边际价格数据，列数: {len(network.buses_t.marginal_price.columns)}")
        print(f"边际价格列名: {list(network.buses_t.marginal_price.columns)}")
        
        # 获取电力节点（通过节点名称识别）
        electricity_buses = [bus for bus in network.buses.index if bus == "electricity"]
        print(f"电力节点: {electricity_buses}")
        
        if len(electricity_buses) > 0:
            price_bus = electricity_buses[0]
            if price_bus in network.buses_t.marginal_price.columns:
                prices = network.buses_t.marginal_price[price_bus]
                ax.plot(prices.index, prices, 
                        label=f'Nodal Price ({price_bus})', linewidth=2, color='green')
                print(f"绘制节点电价: {price_bus}, 价格范围: {prices.min():.2f} - {prices.max():.2f} €/MWh")
            else:
                # 如果没有找到特定节点，绘制所有电力节点的平均电价
                available_buses = [bus for bus in electricity_buses if bus in network.buses_t.marginal_price.columns]
                if len(available_buses) > 0:
                    elec_prices = network.buses_t.marginal_price[available_buses]
                    avg_price = elec_prices.mean(axis=1)
                    ax.plot(avg_price.index, avg_price, 
                            label='Average Nodal Price', linewidth=2, color='green')
                    print(f"绘制平均节点电价，使用节点: {available_buses}, 价格范围: {avg_price.min():.2f} - {avg_price.max():.2f} €/MWh")
                else:
                    # 如果没有任何电力节点，使用第一个可用的节点
                    if len(network.buses_t.marginal_price.columns) > 0:
                        first_bus = network.buses_t.marginal_price.columns[0]
                        prices = network.buses_t.marginal_price[first_bus]
                        ax.plot(prices.index, prices, 
                                label=f'Nodal Price ({first_bus})', linewidth=2, color='green')
                        print(f"绘制节点电价: {first_bus} (备选), 价格范围: {prices.min():.2f} - {prices.max():.2f} €/MWh")
                    else:
                        ax.text(0.5, 0.5, 'No Marginal Price Data Available', 
                                transform=ax.transAxes, ha='center', va='center', fontsize=12)
        else:
            # 如果没有电力节点，尝试其他方法
            if len(network.buses_t.marginal_price.columns) > 0:
                first_bus = network.buses_t.marginal_price.columns[0]
                prices = network.buses_t.marginal_price[first_bus]
                ax.plot(prices.index, prices, 
                        label=f'Nodal Price ({first_bus})', linewidth=2, color='green')
                print(f"绘制节点电价: {first_bus} (无电力节点), 价格范围: {prices.min():.2f} - {prices.max():.2f} €/MWh")
            else:
                ax.text(0.5, 0.5, 'No Marginal Price Data Available', 
                        transform=ax.transAxes, ha='center', va='center', fontsize=12)
    else:
        # 如果没有边际价格数据，显示提示
        print("未找到边际价格数据")
        ax.text(0.5, 0.5, 'Marginal Price Data Not Available', 
                transform=ax.transAxes, ha='center', va='center', fontsize=12)
    
    ax.set_title(f'Electricity Price Time Series{title_suffix}', fontsize=14)
    ax.set_ylabel('Price (€/MWh)', fontsize=12)
    ax.set_xlabel('Time', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

def save_and_show_plot(fig, filename, output_dir="examples/results"):
    """
    保存和显示图片
    fig: matplotlib图形对象
    filename: 文件名
    output_dir: 输出目录
    """
    # 调整布局
    plt.tight_layout()
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    filepath = os.path.join(output_dir, filename)
    
    # 保存图片
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    print(f"图片已保存到: {filepath}")
    
    # 显示图片
    plt.show()

def plot_iteration_results(network, aluminum_usage, ts, p_min_pu, iteration, stage):
    """
    绘制每次迭代的结果：电解槽用能和节点电价
    """
    print(f"绘制第 {iteration} 次迭代 {stage} 阶段结果...")
    
    # 创建图形
    fig, axes = plt.subplots(2, 1, figsize=(15, 10))
    
    # 绘制电解槽用能
    title_suffix = f" (Iteration {iteration}, {stage}, p_min_pu={p_min_pu})"
    plot_aluminum_usage(axes[0], aluminum_usage, ts, title_suffix)
    
    # 绘制节点电价
    plot_nodal_prices(axes[1], network, title_suffix)
    
    # 保存和显示图片
    filename = f"iteration_{iteration}_{stage}_pmin_{p_min_pu}.png"
    save_and_show_plot(fig, filename)

def plot_iterative_results(network, aluminum_usage, ts, p_min_pu, iteration):
    """
    绘制迭代优化结果：电解槽用能和节点电价
    """
    print(f"绘制第 {iteration} 次迭代结果...")
    print(f"电解铝用能数据形状: {aluminum_usage.shape}")
    print(f"电解铝用能列名: {list(aluminum_usage.columns)}")
    
    # 创建图形
    fig, axes = plt.subplots(2, 1, figsize=(15, 10))
    
    # 绘制电解槽用能
    title_suffix = f" (Iteration {iteration}, p_min_pu={p_min_pu})"
    plot_aluminum_usage(axes[0], aluminum_usage, ts, title_suffix)
    
    # 绘制节点电价
    plot_nodal_prices(axes[1], network, title_suffix)
    
    # 保存和显示图片
    filename = f"iterative_result_{iteration}_pmin_{p_min_pu}.png"
    save_and_show_plot(fig, filename)
    
    # 打印统计信息
    print("\n=== Iterative Optimization Results Statistics ===")
    if 'smelter' in aluminum_usage.columns:
        print(f"Average Aluminum Smelter Power: {aluminum_usage['smelter'].mean():.2f} MW")
        print(f"Maximum Aluminum Smelter Power: {aluminum_usage['smelter'].max():.2f} MW")
        print(f"Minimum Aluminum Smelter Power: {aluminum_usage['smelter'].min():.2f} MW")
        print(f"Aluminum Smelter Power Std Dev: {aluminum_usage['smelter'].std():.2f} MW")
    
    # 获取电价统计信息
    if hasattr(network, 'buses_t') and hasattr(network.buses_t, 'marginal_price'):
        electricity_buses = network.buses[network.buses.carrier == "AC"].index
        if len(electricity_buses) > 0:
            available_buses = [bus for bus in electricity_buses if bus in network.buses_t.marginal_price.columns]
            if len(available_buses) > 0:
                elec_prices = network.buses_t.marginal_price[available_buses]
                avg_price = elec_prices.mean(axis=1)
                print(f"Average Nodal Price: {avg_price.mean():.2f} €/MWh")
                print(f"Maximum Nodal Price: {avg_price.max():.2f} €/MWh")
                print(f"Minimum Nodal Price: {avg_price.min():.2f} €/MWh")
            else:
                # 使用第一个可用的节点
                if len(network.buses_t.marginal_price.columns) > 0:
                    first_bus = network.buses_t.marginal_price.columns[0]
                    prices = network.buses_t.marginal_price[first_bus]
                    print(f"Average Nodal Price ({first_bus}): {prices.mean():.2f} €/MWh")
                    print(f"Maximum Nodal Price ({first_bus}): {prices.max():.2f} €/MWh")
                    print(f"Minimum Nodal Price ({first_bus}): {prices.min():.2f} €/MWh")
    
    print(f"Number of Iterations: {iteration}")
    print("=" * 50) 