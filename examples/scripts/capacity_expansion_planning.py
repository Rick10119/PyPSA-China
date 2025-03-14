# capacity_expansion_planning.py

import matplotlib.pyplot as plt
import pandas as pd
import pypsa

plt.style.use("bmh")

# 配置参数
CONFIG = {
    "years": [2025],
    "resolution": 4,
    "al_p_nom": 10 * 1e3,
    # 铝存储成本, 计算方式：1/13.4 是每MWh生产的铝，单位转换为百万欧元/吨，0.74 是一吨铝需要的空间（立方米），0.8 是存储价格（0.8元/平方米/天），1e-6 是转换为百万，7.55 是欧元兑人民币汇率，24 是换算成小时
    "al_marginal_cost_storage": 1/13.4 * 0.74 * 0.8 * 1e-6 / 7.55 / 24,
    # 默认成本
    "default_costs": {
        "FOM": 0,
        "VOM": 0,
        "efficiency": 1,
        "fuel": 0,
        "investment": 0,
        "lifetime": 25,
        "CO2 intensity": 0,
        "discount rate": 0.07,
    }
}

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
    ts.load -= 5.5 # 负荷需求
    ts['aluminum'] = 5.5 # 铝需求，10%负荷
    ts.load *= 1e3 # 负荷需求单位转换为G瓦
    ts.aluminum *= 1e3 # 铝需求单位转换为G瓦
    return ts.resample(f"{CONFIG['resolution']}h").first()

def create_network(costs, ts):
    """创建和配置网络"""
    n = pypsa.Network()
    
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
    add_other_components(n)
    
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

def add_other_components(n):
    """添加其他组件"""
    # 添加铝冶炼设备
    n.add(
        "Link", 
        "smelter", 
        bus0="electricity", 
        bus1="aluminum", 
        p_nom=CONFIG["al_p_nom"], 
        efficiency=1, 
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

def main():
    """主运行函数"""
    results = {}
    
    for year in CONFIG["years"]:
        print(f"\nCalculating for year {year}")
        
        # 数据处理
        costs = process_cost_data(year)
        ts = process_time_series()
        
        # 创建并优化网络
        n = create_network(costs, ts)
        n.optimize(solver_name="gurobi")
        
        # 保存结果
        results[year] = {
            'system_cost': calculate_system_cost(n),
            'generator_capacities': n.generators.p_nom_opt * 1e-3,
            'storage_capacities': n.storage_units.p_nom_opt * 1e-3
        }
    
    # 输出结果
    for year in CONFIG["years"]:
        print(f"\nResults for {year}:")
        print(f"System costs (million €/a):")
        print(results[year]['system_cost'])
        print(f"\nGenerator capacities (GW):")
        print(results[year]['generator_capacities'])
        print(f"\nStorage capacities (GW):")
        print(results[year]['storage_capacities'])

if __name__ == "__main__":
    main()

