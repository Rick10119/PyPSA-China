import pandas as pd
import matplotlib.pyplot as plt
from config import CONFIG

def analyze_emissions(n):
    """分析系统的碳排放情况"""
    
    # 1. 获取每个发电机的排放因子
    carrier_emissions = n.carriers.co2_emissions
    
    # 2. 获取每个发电机的发电量
    if hasattr(n, 'generators_t') and hasattr(n.generators_t, 'p'):
        gen_dispatch = n.generators_t.p
        
        # 3. 计算每个发电机的小时排放量
        emissions = pd.DataFrame(index=gen_dispatch.index)
        
        for generator in gen_dispatch.columns:
            carrier = n.generators.at[generator, 'carrier']
            emissions_factor = carrier_emissions.get(carrier, 0)
            emissions[generator] = gen_dispatch[generator] * emissions_factor
        
        # 4. 计算总排放量
        total_emissions = emissions.sum().sum()  # tonCO2
        total_energy = gen_dispatch.sum().sum()  # MWh
        
        # 5. 计算平均排放强度
        average_intensity = total_emissions / total_energy if total_energy > 0 else 0  # tonCO2/MWh
        
        # 6. 计算月度排放量
        monthly_emissions = emissions.sum(axis=1).resample('M').sum()
        
        # 7. 按发电机类型计算排放量
        emissions_by_generator = emissions.sum()
        emissions_by_carrier = pd.Series(dtype=float)
        
        for generator in emissions_by_generator.index:
            carrier = n.generators.at[generator, 'carrier']
            if carrier in emissions_by_carrier:
                emissions_by_carrier[carrier] += emissions_by_generator[generator]
            else:
                emissions_by_carrier[carrier] = emissions_by_generator[generator]
        
        # 8. 打印结果
        print("\n碳排放分析:")
        print(f"总排放量: {total_emissions/1e6:.2f} 百万吨 CO2")
        print(f"平均排放强度: {average_intensity:.2f} kgCO2/MWh")
        print(f"\n按能源类型排放量 (百万吨 CO2):")
        for carrier, value in emissions_by_carrier.items():
            print(f"{carrier}: {value/1e6:.4f}")
        
        # 9. 绘制排放图表
        # plot_emissions(monthly_emissions, emissions_by_carrier)
        
        # 10. 返回结果
        return {
            'total_emissions': total_emissions,
            'average_intensity': average_intensity,
            'monthly_emissions': monthly_emissions,
            'emissions_by_carrier': emissions_by_carrier
        }
    else:
        print("警告: 无法找到发电机的发电量数据")
        return None

def plot_emissions(monthly_emissions, emissions_by_carrier):
    """绘制排放结果图表"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 绘制月度排放量
    ax1.bar(range(len(monthly_emissions)), monthly_emissions/1e6)
    ax1.set_title('Monthly CO2 Emissions')
    ax1.set_xlabel('Month')
    ax1.set_ylabel('Emissions (Million Tons CO2)')
    ax1.set_xticks(range(len(monthly_emissions)))
    ax1.set_xticklabels([d.strftime('%Y-%m') for d in monthly_emissions.index])
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(axis='y')
    
    # 绘制按能源类型的排放量
    if not emissions_by_carrier.empty:
        # 只绘制非零排放的能源类型
        non_zero = emissions_by_carrier[emissions_by_carrier > 0]
        if not non_zero.empty:
            non_zero = non_zero.sort_values(ascending=False)
            ax2.pie(non_zero/1e6, 
                   labels=non_zero.index, 
                   autopct='%1.1f%%',
                   textprops={'fontsize': 9})
            ax2.set_title('CO2 Emissions by Source')
    
    plt.tight_layout()
    plt.savefig('examples/results/emissions_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

def compare_scenarios(results, scenario_names):
    """比较不同情景的排放结果"""
    # 确保结果列表与情景名称列表长度一致
    if len(results) != len(scenario_names):
        print("错误: 结果数量与情景名称数量不一致")
        return
    
    # 提取总排放量和平均排放强度
    total_emissions = [res['total_emissions']/1e6 if res else 0 for res in results]
    avg_intensity = [res['average_intensity'] if res else 0 for res in results]
    
    # 创建比较图表
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 绘制总排放量比较
    x = range(len(scenario_names))
    ax1.bar(x, total_emissions)
    ax1.set_title('Total CO2 Emissions Comparison')
    ax1.set_xlabel('Scenario')
    ax1.set_ylabel('Emissions (Million Tons CO2)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(scenario_names)
    ax1.grid(axis='y')
    
    # 绘制平均排放强度比较
    ax2.bar(x, avg_intensity)
    ax2.set_title('Average Emission Intensity Comparison')
    ax2.set_xlabel('Scenario')
    ax2.set_ylabel('Intensity (kgCO2/MWh)')
    ax2.set_xticks(x)
    ax2.set_xticklabels(scenario_names)
    ax2.grid(axis='y')
    
    plt.tight_layout()
    plt.savefig('examples/results/emissions_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

# 在主程序中使用这些函数
def main():
    """主程序"""
    # ... 前面的代码 ...
    
    # 针对不同情景的结果存储
    emission_results = []
    scenario_names = []
    
    for year in CONFIG["years"]:
        year_results = {}
        
        # 数据处理
        costs = process_cost_data(year)
        ts = process_time_series()
        
        # 对每个过剩率进行计算
        for excess_rate in CONFIG["al_excess_rate"]:
            # 创建并优化网络
            n = create_network(costs, ts, excess_rate)
            n.optimize(solver_name="gurobi", solver_options={"OutputFlag": 0})
            
            # 分析排放结果
            scenario_name = f"Year {year}, Excess {excess_rate}%"
            scenario_names.append(scenario_name)
            emission_result = analyze_emissions(n)
            emission_results.append(emission_result)
            
            # ... 其他结果处理 ...
    
    # 比较不同情景的排放结果
    compare_scenarios(emission_results, scenario_names) 