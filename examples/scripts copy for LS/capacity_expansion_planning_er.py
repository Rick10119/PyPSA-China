# capacity_expansion_planning.py

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pypsa
from config import CONFIG  # 导入配置
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
    ts.load -= CONFIG["al_demand"] # 负荷需求
    ts['aluminum'] = CONFIG["al_demand"] # 铝需求，10%负荷
    ts.load *= 1e3 # 负荷需求单位GW, 转换为MW
    ts.aluminum *= 1e3 # 铝需求单位转换为MW
    return ts.resample(f"{CONFIG['resolution']}h").first()

def create_network(costs, ts, excess_rate):
    """创建和配置网络"""
    n = pypsa.Network()
    
    # # 计算电解槽容量
    # al_p_nom = CONFIG["al_demand"] * (1 + excess_rate) * 1e3  # 转换为MW
    
    # 添加基础节点
    n.add("Bus", "electricity")
    # n.add("Bus", "aluminum", carrier="AL")
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
    # n.add("Load", "al demand", bus="aluminum", p_set=ts.aluminum)
    
    # 添加发电机组
    add_generators(n, costs, ts)
    
    # 添加存储设施
    add_storage(n, costs)
    
    # # 添加其他组件
    # add_other_components(n, al_p_nom)  # 传入计算得到的电解槽容量
    
    # 添加柔性负荷
    add_FlexibleLoad(n, ts)

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
    
    # =====================
# 添加柔性负荷（产品、仓储和生产过程）
# =====================
def add_FlexibleLoad(n, ts):
    """
    向网络中添加柔性负荷（产品、仓储和生产过程）。
    """
    # 从CSV文件读取产品参数数据
    df = pd.read_csv("examples/data/production_parameters_forLS.csv")
    
    # 提取各列数据
    materials = df['materials'].tolist()
    hierarchy = df['hierarchy'].tolist()
    demand_matrix = df['demand'].tolist()
    net_demand_coefficient = df['net_demand_coefficient'].tolist()
    price_matrix = df['price'].tolist()
    energy_consumption = df['energy_consumption'].tolist()
    initial_inventories = df['initial_inventory'].tolist()
    max_inventory_coefficient = df['max_inventory'].tolist() 

    # 获取关系矩阵（仍从config获取，因为这是一个复杂的矩阵）
    relationship_matrix = CONFIG["relationship_matrix"]

    # 修正价格矩阵单位为1欧元,1欧元=7.55人民币
    price_matrix = np.array(price_matrix) * 10000 / 7.55  

    # 第一部分：添加产品仓储设施
    
    # 添加产品节点（每种物料一个bus，对应文档bus_i）
    for material in materials:
        n.add("Bus", material)

    for i, material in enumerate(materials):
        n.add(
            "Store",
            f"{material}_store",
            bus=material,
            carrier="product",  # 仓储类型
            e_nom=max_inventory_coefficient[i],  # 仓储容量上限（可调）目前假设为7天需求量
            e_initial=initial_inventories[i],  # 初始库存
            e_cyclic=False  # 是否周期性
        )

    # 第二部分：添加生产过程的Link

    # 使用单位用能向量修正关系矩阵，单位向量对应元素除以矩阵的对应列
    relationship_matrix = np.array(relationship_matrix) / np.array(energy_consumption).reshape(-1, 1)

    for j, material in enumerate(materials):

        inputs = []
        efficiencies = []
        for i, input_material in enumerate(materials):
            if relationship_matrix[i][j] > 0 and i != j:
                inputs.append(input_material)
                efficiencies.append(relationship_matrix[i][j])
        # 创建Link元件（每个产品一个生产Link）
        if inputs: 
            link_kwargs = {
                # 电力输入
                "bus0": "electricity",  # 电力输入总线
                "carrier": hierarchy[j],  # 根据配置中的产品层级添加carrier
                "p_nom": energy_consumption[j] * CONFIG["demand_matrix"][j] / 24 * CONFIG["resolution"],  # 每小时的能量需求
                "p_max_pu": 1.0,  # 最大出力比例 # 生产能力上限,假设产能过剩20%
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
            link_kwargs[f"efficiency{len(inputs) + 1}"] = relationship_matrix[j][j]
            link_name = f"{material}_production"
            n.add("Link", link_name, **link_kwargs)

    # 第三部分：添加弹性负荷（负荷切除）
    # 构造需求时间序列：每个产品每天最后一个时段有需求
    AllTime = len(ts.index) # 试验次数，等于时间序列长度
    days = AllTime * CONFIG["resolution"]  // 24  # 天数
    demand_series = np.zeros(AllTime)
    for i in range(days):
        demand_series[i * 24 // CONFIG["resolution"] + 24 // CONFIG["resolution"] - 1] = -1  # 负号表示消耗
    # 为每个产品添加带惩罚成本的负荷
    for i, material in enumerate(materials):
        n.add("Generator",
            f"{material} industrial load",
            bus = material,
            p_max_pu = 0,  # 最大出力为0，只能切负荷
            p_min_pu = demand_series,  # 需求量
            p_nom = demand_matrix[i]*net_demand_coefficient[i],  # 需求规模
            p_nom_extendable = False,  # 允许模型优化此负荷的规模
            marginal_cost = price_matrix[i]  # 切负荷成本（万元）
      )

    # # 添加CO2约束
    # n.add(
    #     "GlobalConstraint",
    #     "CO2Limit",
    #     carrier_attribute="co2_emissions",
    #     sense="<=",
    #     constant=0,
    # )

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
    print(f"{'Excess Rate':^12} | {'Al Capacity':^12} | {'System Cost':^12} | {'Wind':^10} | {'Solar':^10} | {'Battery':^10} | {'H2 Store':^10} | {'Load Shedding Cost':^12}")
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
        
        # 获取切负荷成本
        load_shedding_cost = result['load_shedding_cost']
        # 修正切负荷成本为正数
        load_shedding_cost = (load_shedding_cost + 365*10000/7.55*np.sum(np.array(CONFIG["price_matrix"]) * np.array(CONFIG["demand_matrix"])) ) * 1e-9  # 转换为十亿单位
        
        # 格式化输出
        print(f"{rate * 10:^12.1f} | {result['al_capacity']:^12.2f} | {result['system_cost']:^12.2f} | "
              f"{wind_cap:^10.2f} | {solar_cap:^10.2f} | {battery_cap:^10.2f} | {h2_cap:^10.2f} | {load_shedding_cost:^12.2f}")

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

def plot_loadshedding(n):
    """
    绘制负荷切除相关结果
    按carrier分类，绘制每天的总用能
    仅考虑upstream_product和downstream_product两种carrier
    """
    # 原有的每日电力消耗图
    for carrier in ['upstream_product', 'downstream_product']:
        plt.figure(figsize=(12, 6))
        found_products = False
        
        for material in CONFIG["materials"]:
            link_name = f"{material}_production"
            if (link_name in n.links.index and 
                link_name in n.links_t.p0 and
                "carrier" in n.links.columns and
                n.links.at[link_name, "carrier"] == carrier):
                
                daily_power = n.links_t.p0[link_name].resample('D').sum()
                plt.plot(daily_power.index, daily_power, label=material)
                found_products = True
        
        if found_products:
            plt.rcParams['font.sans-serif'] = ['SimHei']
            carrier_name = "上游产品" if carrier == "upstream_product" else "下游产品"
            plt.title(f'{carrier_name}每日电力消耗')
            plt.xlabel('日期')
            plt.ylabel('每日总电力消耗 (MWh)')
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"examples/results/{carrier}_daily_electricity_usage.png")
        plt.close()
    
    # 新增1: 切负荷情况图 - 每天最后一个时段的标幺化出力（仅第一个月）
    for carrier in ['upstream_product', 'downstream_product']:
        plt.figure(figsize=(14, 8))
        found_products = False
        
        # 获取时间序列信息
        AllTime = len(n.snapshots)
        days_per_month = 30  # 第一个月的天数
        days = min(days_per_month, AllTime * CONFIG["resolution"] // 24)
        
        # 为每个产品收集每天最后时段的标幺化出力
        for material in CONFIG["materials"]:
            gen_name = f"{material} industrial load"
            if (gen_name in n.generators.index and 
                "carrier" not in n.generators.columns or
                any(n.generators.index == gen_name)):
                
                # 检查该产品是否属于当前carrier类型
                material_idx = CONFIG["materials"].index(material)
                if CONFIG["hierarchy"][material_idx] == carrier:
                    
                    # 获取发电机出力和额定功率
                    gen_output = n.generators_t.p[gen_name]
                    p_nom = n.generators.at[gen_name, "p_nom"]
                    
                    # 提取第一个月每天最后一个时段的出力并标幺化
                    daily_last_output = []
                    dates = []
                    
                    for day in range(days):
                        last_hour_idx = day * (24 // CONFIG["resolution"]) + (24 // CONFIG["resolution"]) - 1
                        if last_hour_idx < len(gen_output):
                            # 标幺化：除以p_nom，取绝对值
                            normalized_output = abs(gen_output.iloc[last_hour_idx] / p_nom) if p_nom != 0 else 0
                            daily_last_output.append(normalized_output)
                            dates.append(gen_output.index[last_hour_idx].date())
                    
                    if daily_last_output:
                        plt.plot(dates, daily_last_output, marker='o', label=f'{material}', linewidth=2)
                        found_products = True
        
        if found_products:
            plt.rcParams['font.sans-serif'] = ['SimHei']
            carrier_name = "上游产品" if carrier == "upstream_product" else "下游产品"
            plt.title(f'{carrier_name}切负荷情况（标幺化）- 第一个月每日最后时段')
            plt.xlabel('日期')
            plt.ylabel('标幺化切负荷出力')
            plt.grid(True, alpha=0.3)
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"examples/results/{carrier}_load_shedding_normalized_first_month.png", bbox_inches='tight')
        plt.close()
    
    # 新增2: 产品切负荷成本占比饼图
    # 计算每个产品的切负荷成本
    product_costs = {}
    total_cost = 0
    
    for material in CONFIG["materials"]:
        gen_name = f"{material} industrial load"
        if gen_name in n.generators.index:
            # 计算该产品的切负荷成本（参考433-436行的计算方式）
            cost = (n.generators_t.p[gen_name] * n.generators.at[gen_name, "marginal_cost"]).sum()
            abs_cost = abs(cost)  # 取绝对值
            
            if abs_cost > 0:  # 只包含有切负荷成本的产品
                product_costs[material] = abs_cost
                total_cost += abs_cost
    
    # 绘制饼图
    if product_costs and total_cost > 0:
        plt.figure(figsize=(12, 10))
        
        # 准备饼图数据
        labels = list(product_costs.keys())
        sizes = list(product_costs.values())
        percentages = [cost/total_cost*100 for cost in sizes]
        
        # 设置颜色
        colors = plt.cm.Set3(range(len(labels)))
        
        # 绘制饼图
        wedges, texts, autotexts = plt.pie(sizes, labels=labels, autopct='%1.1f%%', 
                                          colors=colors, startangle=90, 
                                          textprops={'fontsize': 10})
        
        # 设置中文字体和标题
        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.title('产品切负荷成本占比分布', fontsize=16, pad=20)
        
        # 添加图例，显示具体数值
        legend_labels = [f'{label}: {cost:.2e} EUR' for label, cost in product_costs.items()]
        plt.legend(wedges, legend_labels, title="产品及成本", 
                  loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
        
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig("examples/results/load_shedding_cost_distribution.png", bbox_inches='tight')
        plt.close()
        
        # 打印详细的成本信息
        print("\n=== 产品切负荷成本详情 ===")
        print(f"{'产品名称':^15} | {'切负荷成本(EUR)':^15} | {'占比(%)':^10}")
        print("-" * 50)
        for material, cost in sorted(product_costs.items(), key=lambda x: x[1], reverse=True):
            percentage = cost/total_cost*100
            print(f"{material:^15} | {cost:^15.2e} | {percentage:^10.2f}")
        print(f"{'总计':^15} | {total_cost:^15.2e} | {100.0:^10.2f}")

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
                'al_capacity': CONFIG["al_demand"] * (1 + excess_rate),
                # 切负荷成本：每个产品发电机的marginal_cost总和
                'load_shedding_cost': sum(
                        (n.generators_t.p[gen] * n.generators.at[gen, "marginal_cost"]).sum() 
                        for gen in n.generators.index if "industrial load" in gen
                        )   
            }
            
            # 绘制该过剩率下的用电情况
            # plot_results(n, excess_rate)

            # 绘制负荷切除情况
            plot_loadshedding(n)
        
        results[year] = year_results
    
    # 输出结果表格到excel
    df = pd.DataFrame(results)
    df.to_excel("examples/results/capacity_expansion_planning.xlsx")
    print_results_table(results)

if __name__ == "__main__":
    main()

