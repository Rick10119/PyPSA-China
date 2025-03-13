# capacity_expansion_planning.py

import matplotlib.pyplot as plt
import pandas as pd
import pypsa

plt.style.use("bmh")

# 处理不同年份的技术数据和成本
years = [2025, 2030, 2035, 2040]
# years = [2025]
results = {}

for year in years:
    print(f"\nCalculating for year {year}")
    
    # 读取该年份的成本数据
    costs = pd.read_csv(f"data/costs/costs_{year}.csv", index_col=[0, 1])
    
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.unit = costs.unit.str.replace("/kW", "/MW")
    
    defaults = {
        "FOM": 0,
        "VOM": 0,
        "efficiency": 1,
        "fuel": 0,
        "investment": 0,
        "lifetime": 25,
        "CO2 intensity": 0,
        "discount rate": 0.07,
    }
    costs = costs.value.unstack().fillna(defaults)
    
    costs.at["OCGT", "fuel"] = costs.at["gas", "fuel"]
    costs.at["CCGT", "fuel"] = costs.at["gas", "fuel"]
    costs.at["OCGT", "CO2 intensity"] = costs.at["gas", "CO2 intensity"]
    costs.at["CCGT", "CO2 intensity"] = costs.at["gas", "CO2 intensity"]

    def annuity(r, n):
        return r / (1.0 - 1.0 / (1.0 + r) ** n)

    costs["marginal_cost"] = costs["VOM"] + costs["fuel"] / costs["efficiency"]
    annuity = costs.apply(lambda x: annuity(x["discount rate"], x["lifetime"]), axis=1)
    costs["capital_cost"] = (annuity + costs["FOM"] / 100) * costs["investment"]

    # 加载时间序列数据
    file_path = "examples/data/time-series-lecture-2.csv"
    ts = pd.read_csv(file_path, index_col=0, parse_dates=True)  # 直接读取本地文件

    ts.load -= 5.5
    ts['aluminum'] = 5.5
    ts.load *= 1e3
    ts.aluminum *= 1e3
    resolution = 4
    ts = ts.resample(f"{resolution}h").first()

    # 模型初始化
    n = pypsa.Network()
    n.add("Bus", "electricity")
    n.add("Bus", "aluminum", carrier="AL")
    n.set_snapshots(ts.index)
    n.snapshot_weightings.loc[:, :] = resolution

    # 添加技术
    carriers = [
        "onwind",
        "offwind",
        "solar",
        "OCGT",
        "hydrogen storage underground",
        "battery storage",
    ]

    n.add(
        "Carrier",
        carriers,
        color=["dodgerblue", "aquamarine", "gold", "indianred", "magenta", "yellowgreen"],
        co2_emissions=[costs.at[c, "CO2 intensity"] for c in carriers],
    )

    n.add("Load", "demand", bus="electricity", p_set=ts.load)
    n.add("Load", "al demand", bus="aluminum", p_set=ts.aluminum)

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

    n.add(
        "Link", 
        "smelter", 
        bus0="electricity", 
        bus1="aluminum", 
        p_nom=10 * 1e3, 
        efficiency=1, 
    )

    n.add(
        "Store", 
        "aluminum storage", 
        bus="aluminum", 
        e_nom=float('inf'),  
        e_cyclic=True, 
        marginal_cost_storage=1/13.4 * 0.74 * 0.8 * 1e-6 / 7.55 / 24  # 1/13.4 是每MWh生产的铝，单位转换为百万欧元/吨，0.74 是一吨铝需要的空间（立方米），0.8 是存储价格（0.8元/平方米/天），1e-6 是转换为百万，7.55 是欧元兑人民币汇率，24 是换算成小时
    )

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

    # 添加存储单元
    n.add(
        "StorageUnit",
        "battery storage",
        bus="electricity",
        carrier="battery storage",
        max_hours=2,
        capital_cost=costs.at["battery inverter", "capital_cost"]
        + 2 * costs.at["battery storage", "capital_cost"],
        marginal_cost=costs.at["battery inverter", "marginal_cost"],
        efficiency_store=costs.at["battery inverter", "efficiency"],
        efficiency_dispatch=costs.at["battery inverter", "efficiency"],
        p_nom_extendable=True,
        cyclic_state_of_charge=True,
    )

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

    n.add(
        "GlobalConstraint",
        "CO2Limit",
        carrier_attribute="co2_emissions",
        sense="<=",
        constant=0,
    )

    # 优化网络
    n.optimize(solver_name="gurobi")
    
    # 计算系统总成本
    def system_cost(n):
        # 获取所有组件的资本支出
        tsc = n.statistics.capex()
        
        # 确保包含所有carrier的成本
        costs_by_carrier = (
            tsc.groupby(level=1)  # 按carrier分组
            .sum()  # 合计每个carrier的成本
            .div(1e6)  # 转换为百万欧元
        )
        
        return costs_by_carrier.sum()

    # 存储结果
    results[year] = {
        'system_cost': system_cost(n),
        'generator_capacities': n.generators.p_nom_opt * 1e-3,  # 转换为GW
        'storage_capacities': n.storage_units.p_nom_opt * 1e-3  # 转换为GW
    }

# 输出各年份的结果
for year in years:
    print(f"\nResults for {year}:")
    print(f"System costs (million €/a):")
    print(results[year]['system_cost'])
    print(f"\nGenerator capacities (GW):")
    print(results[year]['generator_capacities'])
    print(f"\nStorage capacities (GW):")
    print(results[year]['storage_capacities'])

