# =====================
# 结果可视化
# =====================
"""
该脚本用于绘制优化结果，包括总成本、切负荷率和发电量调度堆叠图。
"""

import os
from tkinter import N
import pandas as pd
import matplotlib.pyplot as plt
from config import CONFIG

def plot_results(n):
    """
    绘制优化结果。
    包括：
    - 最优装机容量
    - 发电量调度堆叠图
    """
    results_dir = "planning_results"
    os.makedirs(results_dir, exist_ok=True)
    
    # 从config中获取材料列表
    materials = CONFIG["materials"]
    demand_matrix = CONFIG["demand_matrix"]
    price_matrix = CONFIG["price_matrix"]

    # 1. 获取优化结果并保存为数组，以便后续分析与绘图
    print("\n===== 优化结果获取 =====")
    
    # 1.1 获取各类型发电机组的最终装机容量
    print("\n1. 各类型机组最终装机容量:")
    
    # 风电和光伏装机容量
    renewable_capacity = {}
    for tech in ["solar", "onwind"]:
        if f"{tech}_plant" in n.generators.index:
            renewable_capacity[tech] = n.generators.at[f"{tech}_plant", "p_nom_opt"]
            print(f"  {tech}: {renewable_capacity[tech]:.2f} MW")
    
    # 火电机组总装机容量（按类型分组）
    thermal_capacity = {}
    for coal_type in CONFIG["coal_types"]:
        total_capacity = 0
        unit_count = CONFIG["thermal_units"][coal_type]["count"]
        # 统计每个类型的所有机组容量
        for i in range(unit_count):
            unit_name = f"{coal_type}_unit_{i}"
            if unit_name in n.generators.index:
                total_capacity += n.generators.at[unit_name, "p_nom"]
        thermal_capacity[coal_type] = total_capacity
        print(f"  {coal_type}: {total_capacity:.2f} MW")
    
    # 水电机组总装机容量
    hydro_capacity = {}
    for hydro_type in CONFIG["hydro_types"]:
        if f"{hydro_type}_aggregated" in n.generators.index:
            hydro_capacity[hydro_type] = n.generators.at[f"{hydro_type}_aggregated", "p_nom_opt"]
            print(f"  {hydro_type}: {hydro_capacity[hydro_type]:.2f} MW")
    
    # 储能电站容量（功率和能量）
    storage_capacity = {}
    storage_energy = {}
    for storage_type in CONFIG["storage_types"]:
        # 获取储能容量
        if f"{storage_type}" in n.stores.index:
            # 获取能量容量
            storage_energy[storage_type] = n.stores.at[f"{storage_type}", "e_nom_opt"]
            
            # 获取功率容量（充电和放电）
            charge_power = 0
            discharge_power = 0
            if f"{storage_type}_charger" in n.links.index:
                charge_power = n.links.at[f"{storage_type}_charger", "p_nom"]
            if f"{storage_type}_discharger" in n.links.index:
                discharge_power = n.links.at[f"{storage_type}_discharger", "p_nom"]
            
            storage_capacity[storage_type] = max(charge_power, discharge_power)
            
            print(f"  {storage_type}: 功率 {storage_capacity[storage_type]:.2f} MW, 能量 {storage_energy[storage_type]:.2f} MWh")
    
    # 1.2 获取各机组的实际出力值（时间序列数据）
    print("\n2. 获取各机组的实际出力时间序列")
    
    # 风电和光伏实际出力
    renewable_output = {}
    for tech in ["solar", "onwind"]:
        tech_name = f"{tech}_plant"
        if tech_name in n.generators.index and tech_name in n.generators_t.p:
            renewable_output[tech] = n.generators_t.p[tech_name]
            print(f"  {tech} 平均出力: {renewable_output[tech].mean():.2f} MW")
    
    # 火电机组实际出力（按类型分组）
    thermal_output = {}
    for coal_type in CONFIG["coal_types"]:
        output_series = pd.Series(0, index=n.snapshots)
        for i in range(CONFIG["thermal_units"][coal_type]["count"]):
            unit_name = f"{coal_type}_unit_{i}"
            if unit_name in n.generators.index and unit_name in n.generators_t.p:
                output_series += n.generators_t.p[unit_name]
        thermal_output[coal_type] = output_series
        print(f"  {coal_type} 平均出力: {output_series.mean():.2f} MW")
    
    # 水电机组实际出力
    hydro_output = {}
    for hydro_type in CONFIG["hydro_types"]:
        hydro_name = f"{hydro_type}_aggregated"
        if hydro_name in n.generators.index and hydro_name in n.generators_t.p:
            hydro_output[hydro_type] = n.generators_t.p[hydro_name]
            print(f"  {hydro_type} 平均出力: {hydro_output[hydro_type].mean():.2f} MW")
    
    # 储能实际充放电
    storage_charge = {}
    storage_discharge = {}
    for storage_type in CONFIG["storage_types"]:
        # 充电出力
        charge_name = f"{storage_type}_charger"
        if charge_name in n.links.index and charge_name in n.links_t.p0:
            storage_charge[storage_type] = n.links_t.p0[charge_name]
            print(f"  {storage_type} 平均充电功率: {storage_charge[storage_type].mean():.2f} MW")
        
        # 放电出力
        discharge_name = f"{storage_type}_discharger"
        if discharge_name in n.links.index and discharge_name in n.links_t.p1:
            storage_discharge[storage_type] = n.links_t.p1[discharge_name]
            print(f"  {storage_type} 平均放电功率: {storage_discharge[storage_type].mean():.2f} MW")
    
    # 获取各个产品生产在每一时段的用能，对每一个f"{material}_electricity_input"link，它的输入功率即为产品生产的用能
    print("\n3. 获取各产品生产的用能时间序列")
    product_energy_input = {}
    for material in materials:
        input_link = f"{material}_electricity_input"
        if input_link in n.links.index and input_link in n.links_t.p0:
            product_energy_input[material] = n.links_t.p0[input_link]
            print(f"  {material} 平均用能: {product_energy_input[material].mean():.2f} MW")
     
    
    # 1.3 获取目标函数数值
    print("\n3. 目标函数数值:")
    
    # 原始目标函数值（未修正）
    objective_value = n.objective
    print(f"  未修正的目标函数值: {objective_value:.2f} 万元")
    
    # 计算产品总增加值（所有产品满足全部需求的总成本）
    total_demand_value = 0
    for i, material in enumerate(materials):
        # 每个产品的总需求值 = 需求量 * 价格 * 需求次数
        days = len(n.snapshots) // 24
        material_demand_cost = demand_matrix[i] * price_matrix[i] * days
        total_demand_value += material_demand_cost
    
    # 修正后的目标函数值（总成本）
    corrected_objective = objective_value + total_demand_value
    print(f"  产品总增加值: {total_demand_value:.2f} 万元")
    print(f"  修正后目标函数值: {corrected_objective:.2f} 万元")
    
    # 计算各部分成本，分析成本组成
    print("\n4. 各部分成本分析:")

    # 4.1 计算各类发电机组的投资成本（仅计算新建容量）
    investment_costs = {}

    # 可再生能源（风电和光伏）投资成本
    renewable_investment = {}
    for tech in ["solar", "onwind"]:
        tech_name = f"{tech}_plant"
        if tech_name in n.generators.index:
            # 计算新增容量（当前容量减去初始容量）
            initial_capacity = CONFIG[f"p_nom_{tech}"]
            optimized_capacity = n.generators.at[tech_name, "p_nom_opt"]
            new_capacity = max(0, optimized_capacity - initial_capacity)  # 确保不为负
            
            # 计算投资成本
            unit_cost = CONFIG["costs"][tech]["capital_cost"]  # 单位投资成本 [万元/MW]
            total_cost = new_capacity * unit_cost  # 总投资成本 [万元]
            
            renewable_investment[tech] = total_cost
            print(f"  {tech} 新增容量: {new_capacity:.2f} MW, 投资成本: {total_cost:.2f} 万元")
    
    # 火电机组投资成本
    thermal_investment = {}
    for coal_type in CONFIG["coal_types"]:
        # 火电机组是否有新增容量需要检查每台机组
        total_new_capacity = 0
        unit_count = CONFIG["thermal_units"][coal_type]["count"]
        unit_capacity = CONFIG["thermal_units"][coal_type]["capacity"]
        
        for i in range(unit_count):
            unit_name = f"{coal_type}_unit_{i}"
            if unit_name in n.generators.index and "p_nom_opt" in n.generators.columns:
                # 计算新增容量
                initial_capacity = unit_capacity
                optimized_capacity = n.generators.at[unit_name, "p_nom_opt"]
                new_capacity = max(0, optimized_capacity - initial_capacity)
                total_new_capacity += new_capacity
        
        # 计算投资成本
        unit_cost = CONFIG["costs"]["thermal"][coal_type]["capital_cost"]
        total_cost = total_new_capacity * unit_cost
        
        thermal_investment[coal_type] = total_cost
        if total_new_capacity > 0:
            print(f"  {coal_type} 新增容量: {total_new_capacity:.2f} MW, 投资成本: {total_cost:.2f} 万元")
    
    # 水电机组投资成本
    hydro_investment = {}
    for hydro_type in CONFIG["hydro_types"]:
        hydro_name = f"{hydro_type}_aggregated"
        if hydro_name in n.generators.index:
            # 计算新增容量
            initial_capacity = CONFIG["hydro_units"][hydro_type]["capacity"] * CONFIG["hydro_units"][hydro_type]["count"]
            optimized_capacity = n.generators.at[hydro_name, "p_nom_opt"]
            new_capacity = max(0, optimized_capacity - initial_capacity)
            
            # 计算投资成本
            unit_cost = CONFIG["costs"]["hydro"][hydro_type]["capital_cost"]
            total_cost = new_capacity * unit_cost
            
            hydro_investment[hydro_type] = total_cost
            if new_capacity > 0:
                print(f"  {hydro_type} 新增容量: {new_capacity:.2f} MW, 投资成本: {total_cost:.2f} 万元")
    
    # 储能投资成本
    storage_investment = {}
    for storage_type in CONFIG["storage_types"]:
        if f"{storage_type}" in n.stores.index:
            # 计算新增能量容量
            storage_params = CONFIG["storage_params"]["matrix"]
            storage_index = CONFIG["storage_types"].index(storage_type)
            
            initial_capacity = storage_params["capacity"][storage_index] * storage_params["count"][storage_index]
            optimized_capacity = n.stores.at[f"{storage_type}", "e_nom_opt"]
            new_capacity = max(0, optimized_capacity - initial_capacity)
            
            # 简化计算：储能成本按照能量容量计算，假设每MWh储能成本为50万元
            # 实际项目中应根据具体的储能类型和成本模型进行计算
            energy_unit_cost = 50  # 假设值 [万元/MWh]
            total_cost = new_capacity * energy_unit_cost
            
            storage_investment[storage_type] = total_cost
            if new_capacity > 0:
                print(f"  {storage_type} 新增能量容量: {new_capacity:.2f} MWh, 投资成本: {total_cost:.2f} 万元")
    
    # 4.2 汇总各类投资成本
    total_renewable_investment = sum(renewable_investment.values())
    total_thermal_investment = sum(thermal_investment.values())
    total_hydro_investment = sum(hydro_investment.values())
    total_storage_investment = sum(storage_investment.values())
    
    total_investment = (total_renewable_investment + total_thermal_investment + 
                        total_hydro_investment + total_storage_investment)
    
    # 将投资成本添加到结果字典中
    investment_costs = {
        'renewable_investment': renewable_investment,
        'thermal_investment': thermal_investment,
        'hydro_investment': hydro_investment,
        'storage_investment': storage_investment,
        'total_renewable_investment': total_renewable_investment,
        'total_thermal_investment': total_thermal_investment,
        'total_hydro_investment': total_hydro_investment,
        'total_storage_investment': total_storage_investment,
        'total_investment': total_investment
    }

    # # 绘制饼图，显示各类投资成本占比
    # plt.figure(figsize=(10, 6))
    # labels = list(investment_costs['renewable_investment'].keys()) + \
    #          list(investment_costs['thermal_investment'].keys()) + \
    #          list(investment_costs['hydro_investment'].keys()) + \
    #          list(investment_costs['storage_investment'].keys())
    # sizes = list(investment_costs['renewable_investment'].values()) + \
    #         list(investment_costs['thermal_investment'].values()) + \
    #         list(investment_costs['hydro_investment'].values()) + \
    #         list(investment_costs['storage_investment'].values())
    # plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
    # plt.title('各类发电机组投资成本占比')

    # 返回获取的结果，以便后续分析和绘图
    results = {
        'renewable_capacity': renewable_capacity,
        'thermal_capacity': thermal_capacity,
        'hydro_capacity': hydro_capacity,
        'storage_capacity': storage_capacity,
        'storage_energy': storage_energy,
        'renewable_output': renewable_output,
        'thermal_output': thermal_output,
        'hydro_output': hydro_output,
        'storage_charge': storage_charge,
        'storage_discharge': storage_discharge,
        'objective_value': objective_value,
        'total_demand_value': total_demand_value,
        'corrected_objective': corrected_objective,
        'product_energy_input': product_energy_input
    }
    
    # 下面可以继续使用results进行绘图分析...

    # 绘制新能源机组出力和最大出力折线图
    plt.figure(figsize=(12, 6))
    for tech, output in renewable_output.items():
        # 实际出力
        plt.plot(output.index, output.values, label=f"{tech} 实际出力", alpha=0.7)
        
        # 最大可能出力 = 装机容量 * 容量因子
        if tech == "solar":
            max_output = renewable_capacity[tech] * n.generators_t.p_max_pu[f"{tech}_plant"]
        elif tech == "onwind":
            max_output = renewable_capacity[tech] * n.generators_t.p_max_pu[f"{tech}_plant"]
        plt.plot(output.index, max_output, label=f"{tech} 最大出力", linestyle='--', alpha=0.4)
    plt.rcParams['font.sans-serif'] = ['SimHei']
    plt.title("新能源机组出力与最大出力时间序列")
    plt.xlabel("时间")
    plt.ylabel("出力 (MW)")
    plt.legend()
    plt.grid()
    plt.show()

    # 绘制火电机组出力折线图，不需要绘制最大出力
    plt.figure(figsize=(12, 6))
    for coal_type, output in thermal_output.items():
        plt.plot(output.index, output.values, label=f"{coal_type} 实际出力", alpha=0.7)
    plt.title("火电机组出力时间序列")
    plt.xlabel("时间")
    plt.ylabel("出力 (MW)")
    plt.legend()
    plt.grid()
    plt.show()

    # 绘制水电机组出力折线图
    plt.figure(figsize=(12, 6))
    for hydro_type, output in hydro_output.items():
        plt.plot(output.index, output.values, label=f"{hydro_type} 实际出力", alpha=0.7)
    plt.title("水电机组出力时间序列")
    plt.xlabel("时间")
    plt.ylabel("出力 (MW)")
    plt.legend()
    plt.grid()
    plt.show()

    # 绘制储能充放电功率堆叠图
    plt.figure(figsize=(12, 6))
    plt.stackplot(storage_charge.index, storage_charge.values, storage_discharge.values, labels=['充电功率', '放电功率'], alpha=0.7)
    plt.title("储能充放电功率堆叠图")
    plt.xlabel("时间")
    plt.ylabel("功率 (MW)")
    plt.legend()
    plt.grid()
    plt.show()

    # 绘制产品生产用能折线图（按config中的raw_materials、products等进行绘制）
    plt.figure(figsize=(12, 6))
    for product, energy in product_energy_input.items():
        plt.plot(energy.index, energy.values, label=f"{product} 用能", alpha=0.7)
    plt.title("产品生产用能时间序列")
    plt.xlabel("时间")
    plt.ylabel("用能 (MW)")
    plt.legend()
    plt.grid()
    plt.show()