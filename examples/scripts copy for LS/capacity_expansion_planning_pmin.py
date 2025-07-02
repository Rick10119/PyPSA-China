# capacity_expansion_planning.py

import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from config import CONFIG  # 导入配置
import os
from analyze_startups import analyze_startups
# 导入排放分析模块
from analyze_emissions import analyze_emissions, plot_emissions, compare_scenarios
plt.style.use("bmh")

def process_cost_data(year):
    """处理成本数据"""
    costs = pd.read_csv(f"data/costs/costs_{year}.csv", index_col=[0, 1])
    
    # 单位转换
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.unit = costs.unit.str.replace("/kW", "/MW")
    
    costs = costs.value.unstack().fillna(CONFIG["default_costs"])
    
    # 设置燃气相关参数
    costs.at["OCGT", "fuel"] = costs.at["gas", "fuel"]
    costs.at["CCGT", "fuel"] = costs.at["gas", "fuel"]
    costs.at["OCGT", "CO2 intensity"] = costs.at["gas", "CO2 intensity"]
    costs.at["CCGT", "CO2 intensity"] = costs.at["gas", "CO2 intensity"]
    
    # 计算成本
    def annuity(r, n):
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    
    costs["marginal_cost"] = costs["VOM"] + costs["fuel"] / costs["efficiency"]
    annuity_values = costs.apply(lambda x: annuity(x["discount rate"], x["lifetime"]), axis=1)
    costs["capital_cost"] = (annuity_values + costs["FOM"] / 100) * costs["investment"]
    
    return costs

def process_time_series():
    """处理时间序列数据"""
    ts = pd.read_csv("examples/data/time-series-lecture-2.csv", index_col=0, parse_dates=True)
    ts.load -= CONFIG["al_demand"] # 负荷需求
    ts['aluminum'] = CONFIG["al_demand"] # 铝需求，10%负荷
    ts.load *= 1e3 # 负荷需求单位GW, 转换为MW
    ts.aluminum *= 1e3 # 铝需求单位转换为MW
    return ts.resample(f"{CONFIG['resolution']}h").first()

def create_network(costs, ts, p_min_pu, excess_rate):
    """创建和配置网络"""
    n = pypsa.Network()
    
    # 计算电解槽容量
    al_p_nom = CONFIG["al_demand"] * (1 + excess_rate) * 1e3  # MW
    
    # 添加基础节点
    n.add("Bus", "electricity")
    n.add("Bus", "aluminum", carrier="AL")
    n.set_snapshots(ts.index)
    n.snapshot_weightings.loc[:, :] = CONFIG["resolution"]
    
    # 添加carriers
    carriers = ["onwind", "offwind", "solar", "OCGT", "hydrogen storage underground", "battery storage"]
    n.add(
        "Carrier",
        carriers,
        color=["dodgerblue", "aquamarine", "gold", "indianred", "magenta", "yellowgreen"],
        co2_emissions=[costs.at[c, "CO2 intensity"] for c in carriers],
    )
    
    # 添加负载
    n.add("Load", "demand", bus="electricity", p_set=ts.load)
    n.add("Load", "al demand", bus="aluminum", p_set=ts.aluminum)
    
    # 添加发电机组
    add_generators(n, costs, ts)
    
    # 添加存储设施
    add_storage(n, costs)
    
    # 添加其他组件
    add_other_components(n, al_p_nom, p_min_pu)  # 传入计算得到的电解槽容量
    
    return n

def add_generators(n, costs, ts):
    """添加发电设备"""
    # 添加OCGT
    n.add(
        "Generator",
        "OCGT",
        bus="electricity",
        carrier="OCGT",
        capital_cost=costs.at["OCGT", "capital_cost"],
        marginal_cost=costs.at["OCGT", "marginal_cost"],
        efficiency=costs.at["OCGT", "efficiency"],
        p_nom_extendable=True,
    )
    
    # 添加可再生能源发电机
    for tech in ["onwind", "offwind", "solar"]:
        n.add(
            "Generator",
            tech,
            bus="electricity",
            carrier=tech,
            p_max_pu=ts[tech],
            capital_cost=costs.at[tech, "capital_cost"],
            marginal_cost=costs.at[tech, "marginal_cost"],
            efficiency=costs.at[tech, "efficiency"],
            p_nom_extendable=True,
        )

def add_storage(n, costs):
    """添加存储设施"""
    # 添加电池存储
    n.add(
        "StorageUnit",
        "battery storage",
        bus="electricity",
        carrier="battery storage",
        max_hours=2,
        capital_cost=costs.at["battery inverter", "capital_cost"] + 2 * costs.at["battery storage", "capital_cost"],
        marginal_cost=costs.at["battery inverter", "marginal_cost"],
        efficiency_store=costs.at["battery inverter", "efficiency"],
        efficiency_dispatch=costs.at["battery inverter", "efficiency"],
        p_nom_extendable=True,
        cyclic_state_of_charge=True,
    )
    
    # 添加氢能存储
    capital_costs = (
        costs.at["electrolysis", "capital_cost"]
        + costs.at["fuel cell", "capital_cost"]
        + 168 * costs.at["hydrogen storage underground", "capital_cost"]
    )
    
    n.add(
        "StorageUnit",
        "hydrogen storage underground",
        bus="electricity",
        carrier="hydrogen storage underground",
        max_hours=168,
        capital_cost=capital_costs,
        marginal_cost=costs.at["electrolysis", "marginal_cost"],
        efficiency_store=costs.at["electrolysis", "efficiency"],
        efficiency_dispatch=costs.at["fuel cell", "efficiency"],
        p_nom_extendable=True,
        cyclic_state_of_charge=True,
    )

def add_other_components(n, al_p_nom, p_min_pu):
    """添加其他组件"""
    # 添加铝冶炼设备
    n.add(
        "Link", 
        "smelter", 
        bus0="electricity", 
        bus1="aluminum", 
        p_nom=al_p_nom,  # 使用计算得到的容量, MW
        efficiency=1, 
        # capital_cost=CONFIG["al_capital_cost"] * al_p_nom, # $  
        start_up_cost=CONFIG["al_start_up_cost"] * al_p_nom, # $
        committable=True,
        p_min_pu=p_min_pu,
    )
    
    # 添加铝存储
    n.add(
        "Store", 
        "aluminum storage", 
        bus="aluminum", 
        e_nom=CONFIG["al_storage_limit"] * al_p_nom,  
        e_cyclic=True, 
        marginal_cost_storage=CONFIG["al_marginal_cost_storage"]
    )
    
    # 添加CO2约束
    n.add(
        "GlobalConstraint",
        "CO2Limit",
        carrier_attribute="co2_emissions",
        sense="<=",
        constant=CONFIG["al_co2_limit"], # kgCO2/MW/year -> kgCO2/MW/hour
    )

def print_results_table(results):
    """将结果整理成表格形式输出"""
    # 创建表头
    print("\n=== System Results Summary ===")
    print(f"{'Excess Rate':^12} | {'Al Capacity':^12} | {'Saved Cost':^12} | {'CO2 Emissions':^12} | {'OCGT':^10} | {'Wind':^10} | {'Solar':^10} | {'Battery':^10} | {'H2 Store':^10}")
    print("-" * 120)
    
    # 获取第一年的结果（因为只有一年）
    year_results = results[CONFIG["years"][0]]
    
    # 按excess_rate排序输出结果
    for rate in sorted(year_results.keys()):
        result = year_results[rate]
        
        # 获取发电机容量
        ocgt_cap = result['generator_capacities']['OCGT']
        wind_cap = (result['generator_capacities']['onwind'] + 
                   result['generator_capacities']['offwind'])
        solar_cap = result['generator_capacities']['solar']
        
        # 获取CO2排放
        co2_emissions = result.get('co2_emissions', 0)  # Get CO2 emissions from results
        
        # 获取储能容量
        battery_cap = result['storage_capacities']['battery storage']
        h2_cap = result['storage_capacities']['hydrogen storage underground']
        
        # 格式化输出
        print(f"{rate * 10:^12.1f} | {result['al_capacity']:^12.2f} | {result['system_cost']:^12.2f} | "
              f"{co2_emissions:^12.2f} | {ocgt_cap:^10.2f} | {wind_cap:^10.2f} | {solar_cap:^10.2f} | {battery_cap:^10.2f} | {h2_cap:^10.2f}")

def plot_results(n, p_min_pu):
    """绘制结果"""
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
    # 先保存图像
    plt.savefig(f"examples/results/aluminum_smelter_usage_{p_min_pu}.png")
    # 然后显示
    plt.show()
    # 关闭图形，释放内存
    plt.close()

def analyze_ramp_constraints(n):
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
    print(f"设定的爬坡限制: {1/CONFIG['al_start_up_time']:.3f} p.u.")
    
    # 绘制功率变化图
    plt.figure(figsize=(12, 6))
    plt.plot(p_diff, 'b-', label='Power Change')
    plt.axhline(y=p_nom / CONFIG['al_start_up_time'], color='r', linestyle='--', label='Ramp Up Limit')
    plt.axhline(y=-p_nom / CONFIG['al_start_up_time'], color='r', linestyle='--', label='Ramp Down Limit')
    plt.title('Aluminum Smelter Power Changes')
    plt.xlabel('Time')
    plt.ylabel('Power Change (MW)')
    plt.legend()
    plt.grid(True)
    plt.savefig('examples/results/smelter_ramp_analysis.png', dpi=300)
    plt.show()
    plt.close()

def main():
    """主运行函数"""
    # 创建结果保存目录
    os.makedirs("examples/results", exist_ok=True)
    
    results = {}
    
    # 用于比较不同场景排放的数据
    emission_results = []
    scenario_names = []
    
    for year in CONFIG["years"]:
        year_results = {}
        
        # 数据处理
        costs = process_cost_data(year)
        ts = process_time_series()
        
        # 对每个p_min_pu进行计算
        for p_min_pu in CONFIG["al_p_min_pu"]:
            # 创建并优化网络
            n = create_network(costs, ts, p_min_pu, CONFIG["al_excess_rate"])
            # 最大求解时间2分钟
            # 设置求解器, 不显示求解信息
            n.optimize(solver_name="gurobi", solver_options={
                "OutputFlag": 0,  # 不显示求解过程
                "IntFeasTol": 1e-5  # 整数可行性容差
            })
            
            # 在优化后调用启动分析
            # analyze_startups(n)
            
            # 分析排放情况
            scenario_name = f"Min Power {p_min_pu*100:.0f}%"
            scenario_names.append(scenario_name)
            emission_result = analyze_emissions(n)
            emission_results.append(emission_result)
            
            # 记录总排放量（百万吨）用于结果表格
            co2_emissions = emission_result['total_emissions'] / 1e6 if emission_result else 0
            
            # 保存结果
            year_results[p_min_pu] = {
                'p_min_pu': p_min_pu,
                'system_cost': CONFIG["original_cost"] - n.objective * 1e-9, # Billion
                'generator_capacities': n.generators.p_nom_opt * 1e-3, # GW
                'storage_capacities': n.storage_units.p_nom_opt * 1e-3, # GW
                'al_capacity': CONFIG["al_demand"] * (1 + CONFIG["al_excess_rate"]),
                'co2_emissions': co2_emissions  # 添加CO2排放结果
            }
            
            # 绘制该过剩率下的用电情况
            # plot_results(n, p_min_pu)
        
        results[year] = year_results
    
    # 比较不同情景的排放结果
    # compare_scenarios(emission_results, scenario_names)
    
    # 输出结果表格到excel
    # 创建更详细的Excel结果
    writer = pd.ExcelWriter("examples/results/capacity_expansion_planning.xlsx", engine='openpyxl')
    
    # # 创建一个汇总结果的数据框
    # summary_data = []
    # for year in results:
    #     for p_min in results[year]:
    #         result = results[year][p_min]
    #         row = {
    #             'Year': year,
    #             'Min Power (%)': p_min * 100,
    #             'System Cost (B€)': result['system_cost'],
    #             'CO2 Emissions (Mt)': result['co2_emissions'],
    #             'Al Capacity (GW)': result['al_capacity'],
    #             'Wind (GW)': result['generator_capacities']['onwind'] + result['generator_capacities']['offwind'],
    #             'Solar (GW)': result['generator_capacities']['solar'],
    #             'OCGT (GW)': result['generator_capacities']['OCGT'],
    #             'Battery (GW)': result['storage_capacities']['battery storage'],
    #             'H2 Storage (GW)': result['storage_capacities']['hydrogen storage underground']
    #         }
    #         summary_data.append(row)
    
    # summary_df = pd.DataFrame(summary_data)
    # summary_df.to_excel(writer, sheet_name='Summary', index=False)
    
    # 保存并关闭Excel
    # writer.save()
    
    # 打印结果表格
    print_results_table(results)

    # 分析爬坡约束
    # analyze_ramp_constraints(n)

if __name__ == "__main__":
    main()

