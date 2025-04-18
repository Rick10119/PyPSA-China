# capacity_expansion_planning.py

import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from config import CONFIG  # 导入配置
import os

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

def create_network(costs, ts, excess_rate):
    """创建和配置网络"""
    n = pypsa.Network()
    
    # 计算电解槽容量
    al_p_nom = CONFIG["al_demand"] * (1 + excess_rate) * 1e3  # 转换为MW
    
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
    add_other_components(n, al_p_nom)  # 传入计算得到的电解槽容量
    
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

def add_other_components(n, al_p_nom):
    """添加其他组件"""
    # 添加铝冶炼设备
    n.add(
        "Link", 
        "smelter", 
        bus0="electricity", 
        bus1="aluminum", 
        p_nom=al_p_nom,  # 使用计算得到的容量
        efficiency=1, 
        # capital_cost=CONFIG["al_capital_cost"],
        start_up_cost=CONFIG["al_start_up_cost"],
        start_up_time=0,
        committable=True,
        p_min_pu=CONFIG["al_p_min_pu"],
    )
    
    # 添加铝存储
    n.add(
        "Store", 
        "aluminum storage", 
        bus="aluminum", 
        e_nom=float('inf'),  
        e_cyclic=True, 
        marginal_cost_storage=CONFIG["al_marginal_cost_storage"]
    )
    
    # 添加CO2约束
    n.add(
        "GlobalConstraint",
        "CO2Limit",
        carrier_attribute="co2_emissions",
        sense="<=",
        constant=0,
    )

def calculate_system_cost(n):
    """计算系统总成本"""
    tsc = n.statistics.capex()
    costs_by_carrier = (
        tsc.groupby(level=1)
        .sum()
        .div(1e6)
    )
    return costs_by_carrier.sum()

def print_results_table(results):
    """将结果整理成表格形式输出"""
    # 创建表头
    print("\n=== System Results Summary ===")
    print(f"{'Excess Rate':^12} | {'Al Capacity':^12} | {'System Cost':^12} | {'Wind':^10} | {'Solar':^10} | {'Battery':^10} | {'H2 Store':^10}")
    print("-" * 82)
    
    # 获取第一年的结果（因为只有一年）
    year_results = results[CONFIG["years"][0]]
    
    # 按excess_rate排序输出结果
    for rate in sorted(year_results.keys()):
        result = year_results[rate]
        
        # 获取发电机容量
        wind_cap = (result['generator_capacities']['onwind'] + 
                   result['generator_capacities']['offwind'])
        solar_cap = result['generator_capacities']['solar']
        
        # 获取储能容量
        battery_cap = result['storage_capacities']['battery storage']
        h2_cap = result['storage_capacities']['hydrogen storage underground']
        
        # 格式化输出
        print(f"{rate * 10:^12.1f} | {result['al_capacity']:^12.2f} | {result['system_cost']:^12.2f} | "
              f"{wind_cap:^10.2f} | {solar_cap:^10.2f} | {battery_cap:^10.2f} | {h2_cap:^10.2f}")

def plot_results(n, excess_rate):
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
    # plt.savefig(f"examples/results/aluminum_smelter_usage_{excess_rate}.png")
    # 然后显示
    # plt.show()
    # 关闭图形，释放内存
    plt.close()

def main():
    """主运行函数"""
    # 创建结果保存目录
    os.makedirs("examples/results", exist_ok=True)
    
    results = {}
    
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
            
            # 保存结果
            year_results[excess_rate] = {
                'system_cost': n.objective * 1e-9, # Billion
                'generator_capacities': n.generators.p_nom_opt * 1e-3,
                'storage_capacities': n.storage_units.p_nom_opt * 1e-3,
                'al_capacity': CONFIG["al_demand"] * (1 + excess_rate)
            }
            
            # 绘制该过剩率下的用电情况
            # plot_results(n, excess_rate)
        
        results[year] = year_results
    
    # 输出结果表格到excel
    df = pd.DataFrame(results)
    df.to_excel("examples/results/capacity_expansion_planning.xlsx")
    print_results_table(results)

if __name__ == "__main__":
    main()

