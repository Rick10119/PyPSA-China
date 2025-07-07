# capacity_expansion_planning.py

from itertools import product
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pypsa
from config import CONFIG  # 导入配置
from plot_results import plot_results  # 导入绘图函数
import os

plt.style.use("bmh")

def process_cost_data(year):
    """处理成本数据"""
    costs = pd.read_csv(f"data/costs/costs_{year}.csv", index_col=[0, 2])
    
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

    return ts.resample(f"{CONFIG['resolution']}h").first()

def process_products():
    """处理产品数据"""

    df_pp = pd.read_csv("examples/data/production_parameters_forLS_modified.csv")
    df_pm = pd.read_csv("examples/data/production_relationship_matrix_modified.csv")

    # 处理产品关系矩阵
    relationship_matrix = df_pm.iloc[:, 1:].values.tolist()
    
    # 将产品参数转换为字典
    products = {
        "materials": df_pp['materials'].tolist(),
        "hierarchy": df_pp['hierarchy'].tolist(),
        "demand": df_pp['demand'].tolist(),
        "net_demand_rate": df_pp['net_demand_rate'].tolist(),
        "price": df_pp['price'].tolist(),
        "energy_consumption": df_pp['energy_consumption'].tolist(),
        "initial_inventory": df_pp['initial_inventory'].tolist(),
        "max_inventory": df_pp['max_inventory'].tolist(),
        "excess_rate": df_pp['excess_rate'].tolist(),
        "type": df_pp['type'].tolist(),
        "relationship_matrix": relationship_matrix,
    }
    
    # 价格矩阵修正为1欧元,1欧元=7.55人民币
    products["price"] = [price * 10000 / 7.55 for price in products["price"]]

    # 计算所有产品的日总电能需求，平均到每个时段
    products["daily_total_load"] = sum(np.array(products["demand"]) * np.array(products["energy_consumption"])) # 每天总电能需求（MWh）
    products["load_t"] = products["daily_total_load"] / 24  # 每小时平均功率（MW）

    # 初始化产品总增加值字典
    products["product_total_profit"] = {}
    
    # 计算每个产品的总增加值
    for i, material in enumerate(products["materials"]):
        products["product_total_profit"][material] = 366 * products["price"][i] * products["demand"][i] * products["net_demand_rate"][i] # 产品i一年的总增加值（欧元）

    # 计算所有产品的日总增加值
    products["daily_total_profit"] = np.sum(np.array(products["price"]) * np.array(products["demand"]) * np.array(products["net_demand_rate"])) # 每天总增加值（欧元）
    # products["fit_t"] = products["daily_total_fit"] / 24  # 每小时增加值（欧元）
    products["total_profit"] = products["daily_total_profit"] * 366 * 1e-9 # 每年总增加值（十亿欧元）

    return products

def create_network(costs, ts, products, scenario):
    """创建和配置网络"""
    n = pypsa.Network()

    # 添加基础节点
    n.add("Bus", "electricity")
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
    
    # 添加发电机组
    add_generators(n, costs, ts)
    
    # 添加存储设施
    add_storage(n, costs)
    
    # 添加负荷
    add_load(n, ts, scenario, products)

    # 添加其他组件
    add_other_components(n, costs, ts)

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

def add_load(n, ts, scenario, products):
    """
    添加负荷到网络中。
    根据不同的场景添加不同类型的负荷。
    """

    # 根据senario添加不同的负荷
    if scenario == 1:
        # 场景1：刚性负荷
        
        # 添加一个总的负荷节点
        n.add("Load", "demand", bus="electricity", p_set=products["load_t"])

    elif scenario == 2:
        # 场景2：简单切负荷

        # 为每个产品添加带惩罚成本的负荷，用负出力发电机建模柔性负荷
        for i, material in enumerate(products["materials"]):
            n.add(
                "Generator",
                f"{material} industrial load",
                bus="electricity",
                p_max_pu=0,  # 最大出力为0，只能切负荷
                p_min_pu=-1,  # 需求量（负号表示消耗）
                p_nom=products["demand"][i] * products["energy_consumption"][i] / 24,  # 需求规模
                p_nom_extendable=False,  # 不允许模型优化此负荷的规模
                marginal_cost=products["price"][i] / products["energy_consumption"][i] * products["net_demand_rate"][i]  # 切负荷成本（欧元）
                # 这里需要注意，price是1单位产品的增加值，energy_consumption是1单位产品的能量消耗，二者相除即为产品i生产用能1单位时产生的边际成本
                # 为保证总增加值与考虑产业链的切负荷一致，需要再乘以净需求率
            )

    elif scenario == 3:
        # 场景3：考虑产业链协同的切负荷

        # 添加产品仓储和生产过程
        add_FlexibleLoad(n, ts, products)


def add_FlexibleLoad(n, ts, products):
    """
    向网络中添加柔性负荷（产品、仓储和生产过程）。
    """

    ## 第一部分：添加产品仓储设施
    # 添加产品节点（每种物料一个bus，对应文档bus_i）
    for material in products["materials"]:
        n.add("Bus", material)
    # 添加产品仓储设施（每个产品添加一个store，对应文档s_i）
    for i, material in enumerate(products["materials"]):
        n.add(
            "Store",
            f"{material}_store",
            bus=material,
            carrier="product",  # 仓储类型
            e_nom=products["max_inventory"][i],  # 仓储容量上限（可调）目前假设为7天需求量
            e_initial=products["initial_inventory"][i],  # 初始库存
            e_cyclic=False  # 是否周期性
        )
    
    ## 第二部分：添加生产过程的Link
    
    # 使用单位用能向量修正关系矩阵，矩阵的每一列除以单位用能向量的对应元素
    # 此时矩阵[i][i]表示输入1MWh电能可生产的产品i的数量；矩阵[i][j]表示产品i生产过程中，输入1MWh电能时消耗的产品j的数量
    # 由于link输入的是功率，而关系矩阵对角元素表示的是每单位电能生产的产品数量，所以需要乘以时间分辨率
    products["relationship_matrix"] = np.array(products["relationship_matrix"]) / np.array(products["energy_consumption"]).reshape(-1, 1) * CONFIG["resolution"]  # 每小时的能量需求

    # 添加产品生产Link（每个产品一个生产Link）
    for j, material in enumerate(products["materials"]):

        inputs = []
        efficiencies = []
        for i, input_material in enumerate(products["materials"]):
            if products["relationship_matrix"][i][j] > 0 and i != j:
                inputs.append(input_material)
                efficiencies.append(products["relationship_matrix"][i][j])
        if inputs: 
            link_kwargs = {
                "bus0": "electricity",  # 电力输入总线
                "carrier": products["hierarchy"][j],  # 根据配置中的产品层级添加carrier
                # p_nom为每个时段的功率需求，时段能量需求为日总用能除每个时段的长度（24/resolution），由于p_nom是输入功率，因此还需要再除以时间分辨率
                # 故最终p_nom为产品的日总能量需求除以24小时
                "p_nom": products["energy_consumption"][j] * products["demand"][j] / 24,  
                "p_max_pu": 0.9*products['excess_rate'][j],  # 最大出力比例
                "p_min_pu": 0.0,  # 最小出力比例
                "p_nom_extendable": False,  
            }
            # 设置输入物料的效率（消耗量为负）
            link_kwargs[f"efficiency"] = -efficiencies[0]
            for idx, input_material in enumerate(inputs):
                link_kwargs[f"bus{idx + 1}"] = input_material
                if idx + 1 != 1:
                    link_kwargs[f"efficiency{idx + 1}"] = -efficiencies[idx]
            # 添加产品输出bus
            link_kwargs[f"bus{len(inputs) + 1}"] = material
            # 产品输出效率
            link_kwargs[f"efficiency{len(inputs) + 1}"] = products["relationship_matrix"][j][j]
            link_name = f"{material}_production"
            n.add("Link", link_name, **link_kwargs)
    
    ## 第三部分：添加弹性负荷（负荷切除）
    
    # 构造需求时间序列：每个产品每天最后一个时段有需求
    AllTime = len(ts.index) # 试验次数，等于时间序列长度
    days = AllTime * CONFIG["resolution"]  // 24  # 天数
    demand_series = np.zeros(AllTime)
    for i in range(days):
        demand_series[i * 24 // CONFIG["resolution"] + 24 // CONFIG["resolution"] - 1] = -1  # 负号表示消耗
    # 为每个产品添加带惩罚成本的负荷
    for i, material in enumerate(products["materials"]):
        n.add("Generator",
            f"{material} industrial load",
            bus = material,
            p_max_pu = 0,  # 最大出力为0，只能切负荷
            p_min_pu = demand_series,  # 需求量
            p_nom = products["demand"][i] * products["net_demand_rate"][i],  # 需求规模
            p_nom_extendable = False,  # 允许模型优化此负荷的规模
            marginal_cost = products["price"][i] / CONFIG["resolution"] # 切负荷成本（欧元）
            # 这里需要注意，price_matrix是1单位产品的增加值，而实际输出的产品会乘以时间分辨率，所以需要除以时间分辨率以确保边际成本正确
      )

def add_other_components(n, costs, ts):

    # 添加CO2约束,取当前值的一半
    n.add(
        "GlobalConstraint",
        "CO2Limit",
        carrier_attribute="co2_emissions",
        sense="<=",
        constant=CONFIG["al_co2_limit"]*0.5,
    )

def print_results_table(results, products):
    """将结果整理成表格形式输出"""

    # 创建表头
    print("\n=== System Results Summary ===")
    print(f"{'Scenario':^12} | {'System Cost':^12} | {'Power System Cost':^16} | {'Load Shedding Cost':^12} | {'Wind':^10} | {'Solar':^10} | {'OCGT':^10} | {'Battery':^10} | {'H2 Store':^10}")
    print("-" * 100)
    
    # 获取第一年的结果（因为只有一年）
    year_results = results[CONFIG["years"][0]]
    
    # 按scenario排序输出结果
    for scenario in sorted(year_results.keys()):
        result = year_results[scenario]
        
        # 获取发电机容量
        wind_cap = (result['generator_capacities']['onwind'] + 
                   result['generator_capacities']['offwind'])
        solar_cap = result['generator_capacities']['solar']

        OCGT_cap = result['generator_capacities']['OCGT']

        # 获取储能容量
        battery_cap = result['storage_capacities']['battery storage']
        h2_cap = result['storage_capacities']['hydrogen storage underground']
        
        # 获取切负荷成本
        load_shedding_cost = result['load_shedding_cost']

        # 这里需要根据scenario的值来计算电力系统成本，并且修正切负荷成本
        # 具体要不要定义在这里还需研判
        if scenario == 1:
            # 电力系统成本等于总成本，无切负荷成本
            power_system_cost = result['system_cost']
        elif scenario == 2 or scenario == 3: 
            # 修正切负荷成本为正数
            load_shedding_cost = load_shedding_cost + products['total_profit'] 
            # 修正总成本
            result['system_cost'] = result['system_cost'] + products['total_profit']
            # 电力系统成本为总成本减去负荷切除成本
            power_system_cost = result['system_cost'] - load_shedding_cost
        else: # 未知情景：最好有一个警告
            pass
                
        # 格式化输出
        print(f"{scenario:^12d} | {result['system_cost']:^12.2f} | "
              f"{power_system_cost:^16.2f} | {load_shedding_cost:^16.2f} | {wind_cap:^10.2f} | {solar_cap:^10.2f} | {OCGT_cap:^10.2f} | {battery_cap:^10.2f} | {h2_cap:^10.2f}")


def main():
    """主运行函数"""
    # 创建结果保存目录
    os.makedirs("examples/scripts_LS/results", exist_ok=True)
    
    # 处理产品参数和关系矩阵
    products = process_products()

    results = {}
    
    for year in CONFIG["years"]:
        year_results = {}
        
        # 数据处理
        costs = process_cost_data(year)
        ts = process_time_series()
        
        # 三种情景循环进行三次优化，结果保存到results字典中，差异只在负荷的定义上，因此添加负荷的函数应该是根据情景来定义的

        # 对每个情景进行计算
        for scenario in CONFIG["scenario"]:
            # 创建并优化网络
            n = create_network(costs, ts, products, scenario)
            n.optimize(solver_name="gurobi", solver_options={"OutputFlag": 0})

            # 计算负荷切除成本，每个产品发电机的marginal_cost总和,这个放在这里不是很好看，不过它是我们的主要关注点，放在main里也许也可以
            if scenario == 2:
                # 场景2简单切负荷，直接计算负荷切除成本
                # 这里需要注意，因为generators_t.p是出力（单位是功率），而marginal_cost是单位能量的成本，所以需要将出力乘以时间分辨率（CONFIG["resolution"]）来得到能量
                load_shedding_cost = CONFIG["resolution"]*sum(
                            (n.generators_t.p[gen] * n.generators.at[gen, "marginal_cost"]).sum() 
                            for gen in n.generators.index if "industrial load" in gen
                )
            elif scenario == 3:
                # 场景3考虑产业链协同的切负荷，因为marginal_cost除以了resolution, 所以这里需要乘以resolution
                load_shedding_cost = CONFIG["resolution"]*sum(
                            (n.generators_t.p[gen] * n.generators.at[gen, "marginal_cost"]).sum() 
                            for gen in n.generators.index if "industrial load" in gen
                )
            elif scenario == 1:
                load_shedding_cost = 0

            # 保存结果
            year_results[scenario] = {
                'system_cost': n.objective * 1e-9, # Billion
                'generator_capacities': n.generators.p_nom_opt * 1e-3,
                'storage_capacities': n.storage_units.p_nom_opt * 1e-3,               
                'load_shedding_cost': load_shedding_cost * 1e-9
            }
            
            # 根据scenario绘制不同的图像
            plot_results(n, products, scenario)

        results[year] = year_results
    
    # 输出结果表格到excel
    df = pd.DataFrame(results)
    df.to_excel("examples/scripts_LS/results/capacity_expansion_planning.xlsx")
    print_results_table(results, products)

if __name__ == "__main__":
    main()

