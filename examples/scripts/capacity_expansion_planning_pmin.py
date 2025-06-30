import matplotlib.pyplot as plt
import pandas as pd
import pypsa
import sys
import os
import argparse
import yaml

# 添加scripts目录到Python路径
sys.path.append(os.path.join(os.path.dirname(__file__)))

from config import CONFIG  # 导入配置
from analyze_startups import analyze_startups
# 导入排放分析模块
from analyze_emissions import analyze_emissions, plot_emissions, compare_scenarios
# 导入绘图模块
from plot_capacity_expansion import (
    plot_results, 
    analyze_ramp_constraints, 
    plot_capacity_comparison,
    plot_network_summary,
    plot_time_series_analysis,
    create_summary_plots
)

plt.style.use("bmh")

def load_config(config_file="config.yaml"):
    """加载配置文件"""
    with open(config_file, 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)
    return config

def get_solver_options(config):
    """从配置文件中获取求解器选项"""
    solving_config = config.get('solving', {})
    solver_name = solving_config.get('solver', {}).get('name', 'gurobi')
    solver_options_name = solving_config.get('solver', {}).get('options', 'default')
    solver_options = solving_config.get('solver_options', {}).get(solver_options_name, {})
    
    return solver_name, solver_options

def process_cost_data(year):
    """处理成本数据"""
    costs = pd.read_csv(f"examples/data/costs/costs_{year}.csv", index_col=[0, 1])
    
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
    n.add("Bus", "aluminum", carrier="aluminum")
    n.set_snapshots(ts.index)
    n.snapshot_weightings.loc[:, :] = CONFIG["resolution"]
    
    # 添加carriers
    carriers = ["onwind", "offwind", "solar", "OCGT", "hydrogen storage underground", "battery storage"]
    colors = ["dodgerblue", "aquamarine", "gold", "indianred", "magenta", "yellowgreen"]
    
    for i, carrier in enumerate(carriers):
        n.add(
            "Carrier",
            carrier,
            color=colors[i],
            co2_emissions=costs.at[carrier, "CO2 intensity"],
        )
    
    # 添加铝相关的carriers
    n.add("Carrier", "aluminum")
    n.add("Carrier", "aluminum smelter")
    n.add("Carrier", "aluminum storage")
    
    # 添加负载
    n.add("Load", "demand", bus="electricity", p_set=ts.load)
    n.add("Load", "al demand", bus="aluminum", p_set=ts.aluminum)
    
    # 添加切负荷发电机（如果启用）
    if CONFIG["enable_load_shedding"]:
        # 添加切负荷carrier
        n.add("Carrier", "load", color="#dd2e23", nice_name="Load shedding")
        
        # 为电力负荷添加切负荷发电机
        n.add(
            "Generator",
            "demand load",
            bus="electricity",
            carrier="load",
            sign=1e-3,  # 调整单位，p和p_nom以kW为单位而不是MW
            marginal_cost=CONFIG["load_shedding_cost"],  # 欧元/kWh
            p_nom=1e9,  # kW，设置一个很大的容量
        )
        
        # 为铝负荷添加切负荷发电机
        n.add(
            "Generator",
            "al demand load",
            bus="aluminum",
            carrier="load",
            sign=1e-3,  # 调整单位，p和p_nom以kW为单位而不是MW
            marginal_cost=CONFIG["al_load_shedding_cost"],  # 欧元/kWh
            p_nom=1e9,  # kW，设置一个很大的容量
        )
    
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
        carrier="aluminum smelter",
        p_nom=al_p_nom,  # 使用计算得到的容量, MW
        efficiency=1, 
        # capital_cost=CONFIG["al_capital_cost"] * al_p_nom, # $  
        start_up_cost=CONFIG["al_start_up_cost"] * al_p_nom, # $
        committable=CONFIG["al_committable"],
        p_min_pu=p_min_pu if CONFIG["al_committable"] else 0,
    )
    print(CONFIG["al_start_up_cost"] * al_p_nom)
    print(CONFIG["al_start_up_cost"], al_p_nom)
    
    # 添加铝存储
    n.add(
        "Store", 
        "aluminum storage", 
        bus="aluminum", 
        carrier="aluminum storage",
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
    print(f"{'Excess Rate':^12} | {'Al Capacity':^12} | {'Saved Cost':^12} | {'CO2 Emissions':^12} | {'Load Shed':^10} | {'Al Load Shed':^12} | {'OCGT':^10} | {'Wind':^10} | {'Solar':^10} | {'Battery':^10} | {'H2 Store':^10}")
    print("-" * 140)
    
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
        
        # 获取切负荷信息
        load_shed = result.get('load_shedding', {}).get('demand', 0) if result.get('load_shedding') else 0
        al_load_shed = result.get('load_shedding', {}).get('al demand', 0) if result.get('load_shedding') else 0
        
        # 格式化输出
        print(f"{rate * 10:^12.1f} | {result['al_capacity']:^12.2f} | {result['system_cost']:^12.2f} | "
              f"{co2_emissions:^12.2f} | {load_shed:^10.2f} | {al_load_shed:^12.2f} | {ocgt_cap:^10.2f} | {wind_cap:^10.2f} | {solar_cap:^10.2f} | {battery_cap:^10.2f} | {h2_cap:^10.2f}")

def main(config_file="config.yaml"):
    """主运行函数"""
    # 加载配置文件
    config = load_config(config_file)
    
    # 获取求解器设置
    solver_name, solver_options = get_solver_options(config)
    
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
            
            # 使用优化函数
            n.optimize(solver_name=solver_name, **solver_options)
            
            # 在优化后调用启动分析
            if CONFIG["al_committable"]:
                analyze_startups(n)
            
            # 分析排放情况
            scenario_name = f"Min Power {p_min_pu*100:.0f}%"
            scenario_names.append(scenario_name)
            emission_result = analyze_emissions(n)
            emission_results.append(emission_result)
            
            # 记录总排放量（百万吨）用于结果表格
            co2_emissions = emission_result['total_emissions'] / 1e6 if emission_result else 0
            
            # 收集切负荷数据
            load_shedding_data = {}
            if CONFIG["enable_load_shedding"]:
                # 计算切负荷量（GWh）
                # 由于使用了sign=1e-3，切负荷发电机的出力需要乘以1e3来转换为MW
                demand_load_shed = n.generators_t.p['demand load'].sum() * 1e3 if 'demand load' in n.generators_t.p.columns else 0
                al_demand_load_shed = n.generators_t.p['al demand load'].sum() * 1e3 if 'al demand load' in n.generators_t.p.columns else 0
                
                load_shedding_data = {
                    'demand': demand_load_shed,
                    'al demand': al_demand_load_shed
                }
                # 转换为GWh
                load_shedding_data = {k: v * CONFIG["resolution"] / 1e3 for k, v in load_shedding_data.items()}
            
            # 安全获取目标值
            try:
                objective_value = n.objective if hasattr(n, 'objective') else 0.0
            except:
                objective_value = 0.0
            
            # 保存结果
            year_results[p_min_pu] = {
                'p_min_pu': p_min_pu,
                'system_cost': CONFIG["original_cost"] - objective_value * 1e-9, # Billion
                'generator_capacities': n.generators.p_nom_opt * 1e-3, # GW
                'storage_capacities': n.storage_units.p_nom_opt * 1e-3, # GW
                'al_capacity': CONFIG["al_demand"] * (1 + CONFIG["al_excess_rate"]),
                'co2_emissions': co2_emissions,  # 添加CO2排放结果
                'load_shedding': load_shedding_data  # 添加切负荷数据
            }
            
            # 绘制该过剩率下的用电情况
            plot_results(n, p_min_pu)
            
            # 可选：绘制网络概览图
            # plot_network_summary(n, p_min_pu)
            
            # 可选：绘制时间序列分析图
            # plot_time_series_analysis(n, p_min_pu)
        
        results[year] = year_results
    
    # 比较不同情景的排放结果
    # compare_scenarios(emission_results, scenario_names)
    
    # 创建汇总图表
    # create_summary_plots(results, CONFIG)
    
    # 打印结果表格
    print_results_table(results)
    
    # 分析爬坡约束
    # analyze_ramp_constraints(n, CONFIG)

if __name__ == "__main__":
    # 可以通过命令行参数指定config文件
    parser = argparse.ArgumentParser(description='Capacity expansion planning with aluminum smelter (Fixed Version)')
    parser.add_argument('--config', type=str, default="config.yaml", help='Path to config file (default: config.yaml)')
    args = parser.parse_args()
    
    main(config_file=args.config) 