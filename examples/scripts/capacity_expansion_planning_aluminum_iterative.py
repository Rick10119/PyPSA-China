import matplotlib.pyplot as plt
import pandas as pd
import pypsa
import sys
import os
import argparse
import yaml
import copy
import numpy as np
import time  # 添加时间模块

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
    create_summary_plots,
    plot_aluminum_usage,
    plot_nodal_prices,
    save_and_show_plot,
    plot_iteration_results,
    plot_iterative_results
)

plt.style.use("bmh")

"""
电解铝迭代优化算法 - 重构版本

代码重构改进：
1. 创建了通用的 create_base_network() 函数，统一处理基础网络结构创建
2. 提取了 add_aluminum_components() 函数，模块化处理电解铝相关组件
3. 提取了 add_virtual_generator() 函数，专门处理虚拟发电机
4. 减少了代码重复，提高了可维护性和可读性
5. 支持不同类型的网络创建（完整网络 vs 电解铝优化网络）

模块化设计：
- 网络创建和优化逻辑保留在主文件中
- 所有绘图函数已移动到 plot_capacity_expansion.py 模块中
- 通过导入语句使用绘图函数，保持代码结构清晰

主要函数说明：
- create_base_network(): 创建基础网络结构，支持两种模式
- add_aluminum_components(): 添加电解铝相关组件，支持灵活配置
- create_network(): 创建完整网络（包含所有发电和存储设备）
- solve_aluminum_optimization(): 创建电解铝优化网络（仅包含铝相关组件）
- solve_network_iterative(): 电解铝迭代优化算法主函数

绘图函数（在 plot_capacity_expansion.py 中）：
- plot_aluminum_usage(): 绘制电解铝用能模式
- plot_nodal_prices(): 绘制节点电价时间序列
- plot_iteration_results(): 绘制迭代过程结果
- plot_iterative_results(): 绘制最终迭代结果
"""

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

def create_base_network(ts, costs=None, network_type="full"):
    """
    创建基础网络结构
    network_type: "full" - 完整网络（包含所有发电和存储设备）
                  "aluminum_only" - 仅电解铝优化网络（只包含铝相关组件和虚拟发电机）
    """
    n = pypsa.Network()
    
    # 设置时间序列
    n.set_snapshots(ts.index)
    n.snapshot_weightings.loc[:, :] = CONFIG["resolution"]
    
    # 添加基础节点
    n.add("Bus", "electricity")
    n.add("Bus", "aluminum")
    
    # 添加基础carriers
    n.add("Carrier", "aluminum")
    n.add("Carrier", "aluminum smelter")
    n.add("Carrier", "aluminum storage")
    
    if network_type == "full":
        # 添加完整网络的carriers
        carriers = ["onwind", "offwind", "solar", "OCGT", "hydrogen storage underground", "battery storage"]
        colors = ["dodgerblue", "indianred", "aquamarine", "gold", "indianred", "magenta", "yellowgreen", "gray"]
        
        for i, carrier in enumerate(carriers):
            n.add(
                "Carrier",
                carrier,
                color=colors[i],
                co2_emissions=costs.at[carrier, "CO2 intensity"],
            )
    elif network_type == "aluminum_only":
        # 添加电解铝优化网络的carriers
        n.add("Carrier", "virtual")
    
    return n

def add_aluminum_components(n, al_p_nom, p_min_pu, aluminum_commitment=False, include_storage=True, include_co2_constraint=True, verbose=True):
    """
    添加电解铝相关组件
    n: PyPSA网络对象
    al_p_nom: 电解铝容量 (MW)
    p_min_pu: 最小出力比例
    aluminum_commitment: 是否启用启停约束
    include_storage: 是否包含铝存储
    include_co2_constraint: 是否包含CO2约束
    verbose: 是否输出详细信息
    """
    # 添加铝冶炼设备
    n.add(
        "Link", 
        "smelter", 
        bus0="electricity", 
        bus1="aluminum", 
        carrier="aluminum smelter",
        p_nom=al_p_nom,
        efficiency=1, 
        start_up_cost=CONFIG["al_start_up_cost"] * al_p_nom,
        committable=aluminum_commitment,
        p_min_pu=p_min_pu if aluminum_commitment else 0,
    )
    
    if aluminum_commitment and verbose:
        print(f"电解铝启停约束: {aluminum_commitment}")
        print(f"启动成本: {CONFIG['al_start_up_cost'] * al_p_nom}")
        print(f"最小出力比例: {p_min_pu}")
    
    # 添加铝存储（可选）
    if include_storage:
        n.add(
            "Store", 
            "aluminum storage", 
            bus="aluminum", 
            carrier="aluminum storage",
            e_nom=CONFIG["al_storage_limit"] * al_p_nom,  
            e_cyclic=True, 
            # marginal_cost_storage=CONFIG["al_marginal_cost_storage"]
        )
    
    # 添加CO2约束（可选）
    if include_co2_constraint:
        n.add(
            "GlobalConstraint",
            "CO2Limit",
            carrier_attribute="co2_emissions",
            sense="<=",
            constant=CONFIG["al_co2_limit"],
        )

def add_virtual_generator(n, nodal_prices):
    """添加虚拟发电机"""
    n.add("Generator",
          "virtual_gen_electricity",
          bus="electricity",
          carrier="virtual",
          p_nom=1e6,  # 大容量
          marginal_cost=nodal_prices)  # 使用节点电价作为边际成本

def create_network(costs, ts, p_min_pu, excess_rate, aluminum_commitment=False):
    """创建和配置完整网络"""
    # 使用通用函数创建基础网络
    n = create_base_network(ts, costs, network_type="full")
    
    # 计算电解槽容量
    al_p_nom = CONFIG["al_demand"] * (1 + excess_rate) * 1e3  # MW
    
    # 添加负载
    n.add("Load", "demand", bus="electricity", p_set=ts.load)
    n.add("Load", "al demand", bus="aluminum", p_set=ts.aluminum)
    
    # 添加发电机组
    add_generators(n, costs, ts)
    
    # 添加存储设施
    add_storage(n, costs)
    
    # 添加电解铝相关组件
    add_aluminum_components(n, al_p_nom, p_min_pu, aluminum_commitment, verbose=True)
    
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

def safe_optimize(n, solver_name, solver_options):
    """安全地执行优化，处理版本兼容性问题"""
    # 设置求解器选项来最小化输出
    if solver_name.lower() == 'gurobi':
        # Gurobi求解器选项：关闭输出
        quiet_options = {
            'OutputFlag': 0,  # 关闭输出
            'LogToConsole': 0,  # 关闭控制台日志
            'LogFile': '',  # 关闭日志文件
        }
        # 合并用户选项和静默选项
        final_options = {**quiet_options, **solver_options}
    else:
        # 其他求解器保持原有选项
        final_options = solver_options
    
    # 尝试使用新版本的优化方法
    n.optimize(solver_name=solver_name, **final_options)
    return n.objective

def solve_aluminum_optimization(n, costs, ts, p_min_pu, excess_rate, nodal_prices=None, verbose=True):
    """
    电解铝最优运行问题求解
    基于节点电价，运行以满足铝需求为约束的电解铝最优运行问题
    """
    if verbose:
        print("步骤2: 运行电解铝最优运行问题")
    
    # 使用通用函数创建基础网络
    n_al = create_base_network(ts, network_type="aluminum_only")
    
    # 计算电解槽容量
    al_p_nom = CONFIG["al_demand"] * (1 + excess_rate) * 1e3  # MW
    
    # 添加虚拟发电机来提供电力（基于节点电价）
    add_virtual_generator(n_al, nodal_prices)
    
    # 添加电解铝相关组件（不包含CO2约束）
    add_aluminum_components(n_al, al_p_nom, p_min_pu, aluminum_commitment=True, 
                           include_storage=True, include_co2_constraint=False, verbose=verbose)
    
    # 添加铝负荷
    n_al.add("Load", "al demand", bus="aluminum", p_set=ts.aluminum)
    
    # 求解电解铝优化问题
    if verbose:
        print("求解电解铝优化问题...")
    # 使用安全的优化函数
    safe_optimize(n_al, "gurobi", {})
    
    # 提取电解铝用能模式
    aluminum_usage = n_al.links_t.p0[['smelter']].copy()
    return aluminum_usage

def solve_network_iterative(costs, ts, p_min_pu, excess_rate, max_iterations=10, convergence_tolerance=0.01, verbose=True):
    """
    电解铝迭代优化算法
    1. 使用连续化电解铝模型求解，得到节点电价
    2. 基于节点电价，运行电解铝最优运行问题
    3. 固定电解铝用能，求解剩下的优化问题
    4. 重复2-3，直到收敛
    
    收敛条件：目标函数变化的相对值小于convergence_tolerance（默认1%）
    
    参数:
    verbose: 是否输出详细信息
    """
    if verbose:
        print("开始电解铝迭代优化算法")
    
    # 记录总开始时间
    total_start_time = time.time()
    
    # 初始化电解铝用能模式和目标函数值
    aluminum_usage = None
    previous_objective = None
    iteration = 0
    final_network = None
    final_aluminum_usage = None
    
    # 记录每次迭代的时间
    iteration_times = []
    
    while iteration < max_iterations:
        iteration += 1
        iteration_start_time = time.time()  # 记录每次迭代开始时间
        
        if verbose:
            print(f"开始第 {iteration} 次迭代")
        
        # 步骤1: 使用连续化电解铝模型求解，得到节点电价
        if verbose:
            print("步骤1: 使用连续化电解铝模型求解")
        
        # 创建网络，禁用电解铝启停约束
        n = create_network(costs, ts, p_min_pu, excess_rate, aluminum_commitment=False)
        
        # 如果有固定的电解铝用能，需要添加约束
        if aluminum_usage is not None:
            if verbose:
                print("固定电解铝用能，求解剩余优化问题")
            # 固定电解铝用能 - 使用p_set来固定smelter的link出力
            fixed_aluminum_power = aluminum_usage['smelter'].values
            # 修改al demand设定，和p_set一样
            n.loads_t.p_set['al demand'] = fixed_aluminum_power
            # 修改smelter的link出力，和p_set一样
            n.links_t.p_set['smelter'] = fixed_aluminum_power
        
        # 求解网络
        if verbose:
            print("求解网络...")
        # 使用安全的优化函数，确保计算边际价格
        current_objective = safe_optimize(n, "gurobi", {})
        
        # 提取当前迭代的节点电价（用于电解铝优化）
        current_nodal_prices = None
        if hasattr(n, 'buses_t') and hasattr(n.buses_t, 'marginal_price'):
            electricity_buses = [bus for bus in n.buses.index if bus == "electricity"]
            if len(electricity_buses) > 0:
                price_bus = electricity_buses[0]
                current_nodal_prices = n.buses_t.marginal_price[price_bus]
                if verbose:
                    print(f"提取节点电价: {price_bus}")
        
        # 步骤2: 基于节点电价，运行电解铝最优运行问题
        new_aluminum_usage = solve_aluminum_optimization(n, costs, ts, p_min_pu, excess_rate, current_nodal_prices, verbose)
        
        # 步骤3: 检查收敛性 - 基于目标函数变化的相对值
        if previous_objective is not None and current_objective is not None:
            # 计算目标函数变化的相对值
            objective_change = abs(current_objective - previous_objective)
            relative_change = objective_change / abs(previous_objective) if abs(previous_objective) > 1e-10 else float('inf')
            
            if verbose:
                print(f"目标函数变化统计:")
                print(f"  当前目标函数值: {current_objective:.6e}")
                print(f"  上次目标函数值: {previous_objective:.6e}")
                print(f"  绝对变化: {objective_change:.6e}")
                print(f"  相对变化: {relative_change:.6f} ({relative_change*100:.2f}%)")
                print(f"  收敛阈值: {convergence_tolerance:.6f} ({convergence_tolerance*100:.2f}%)")
            
            if relative_change < convergence_tolerance:
                if verbose:
                    print(f"算法收敛，在第 {iteration} 次迭代后停止")
                    print(f"目标函数相对变化 {relative_change:.6f} ({relative_change*100:.2f}%) < 收敛阈值 {convergence_tolerance:.6f} ({convergence_tolerance*100:.2f}%)")
                final_network = n
                final_aluminum_usage = new_aluminum_usage
                break
        else:
            # 第一次迭代，保存目标函数值用于下次比较
            if current_objective is not None:
                previous_objective = current_objective
        
        # 更新电解铝用能和目标函数值
        aluminum_usage = new_aluminum_usage
        previous_objective = current_objective
        final_network = n
        final_aluminum_usage = new_aluminum_usage
        
        # 记录本次迭代时间
        iteration_time = time.time() - iteration_start_time
        iteration_times.append(iteration_time)
        
        if verbose:
            print(f"第 {iteration} 次迭代完成，耗时: {iteration_time:.2f} 秒")
    
    # 计算总时间
    total_time = time.time() - total_start_time
    
    if iteration >= max_iterations:
        if verbose:
            print(f"达到最大迭代次数 {max_iterations}，算法未完全收敛")
    
    # 输出时间统计
    print(f"\n=== 迭代时间统计 ===")
    print(f"总迭代次数: {iteration}")
    print(f"总耗时: {total_time:.2f} 秒 ({total_time/60:.2f} 分钟)")
    if iteration_times:
        print(f"平均每次迭代耗时: {np.mean(iteration_times):.2f} 秒")
        print(f"最快迭代耗时: {np.min(iteration_times):.2f} 秒")
        print(f"最慢迭代耗时: {np.max(iteration_times):.2f} 秒")
        print(f"迭代时间详情:")
        for i, t in enumerate(iteration_times, 1):
            print(f"  第{i}次迭代: {t:.2f} 秒")
    
    if verbose:
        print(f"电解铝迭代优化算法完成，共进行 {iteration} 次迭代")
    
    # 绘制最后一次迭代的结果
    if final_network is not None and final_aluminum_usage is not None:
        plot_iterative_results(final_network, final_aluminum_usage, ts, p_min_pu, iteration)
    
    # 返回最终的网络结果和时间统计
    return final_network, final_aluminum_usage, {
        'total_time': total_time,
        'iteration_count': iteration,
        'iteration_times': iteration_times,
        'avg_iteration_time': np.mean(iteration_times) if iteration_times else 0
    }

def print_results_table(results):
    """将结果整理成表格形式输出"""
    # 创建表头
    print("\n=== System Results Summary ===")
    print(f"{'Excess Rate':^12} | {'Al Capacity':^12} | {'Saved Cost':^12} | {'CO2 Emissions':^12} | {'Load Shed':^10} | {'Al Load Shed':^12} | {'OCGT':^10} | {'Wind':^10} | {'Solar':^10} | {'Battery':^10} | {'H2 Store':^10} | {'Total Time':^10} | {'Iterations':^10}")
    print("-" * 160)
    
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
        
        # 获取时间统计
        total_time = result.get('total_time', 0)
        iteration_count = result.get('iteration_count', 0)
        
        # 格式化输出
        print(f"{rate * 10:^12.1f} | {result['al_capacity']:^12.2f} | {result['system_cost']:^12.2f} | "
              f"{co2_emissions:^12.2f} | {load_shed:^10.2f} | {al_load_shed:^12.2f} | {ocgt_cap:^10.2f} | {wind_cap:^10.2f} | {solar_cap:^10.2f} | {battery_cap:^10.2f} | {h2_cap:^10.2f} | {total_time:^10.1f} | {iteration_count:^10d}")
    
    # 添加总体时间统计
    print("\n=== 总体时间统计 ===")
    total_times = [result.get('total_time', 0) for result in year_results.values()]
    total_iterations = [result.get('iteration_count', 0) for result in year_results.values()]
    
    if total_times:
        print(f"总计算时间: {sum(total_times):.2f} 秒 ({sum(total_times)/60:.2f} 分钟)")
        print(f"总迭代次数: {sum(total_iterations)}")

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
            print(f"\n=== 测试电解铝迭代优化算法 ===")
            print(f"参数: p_min_pu = {p_min_pu}, excess_rate = {CONFIG['al_excess_rate']}")
            
            # 使用迭代优化算法
            n, aluminum_usage, time_stats = solve_network_iterative(
                costs, ts, p_min_pu, CONFIG["al_excess_rate"],
                max_iterations=CONFIG["max_iterations"],  # 减少迭代次数用于测试
                convergence_tolerance=CONFIG["convergence_tolerance"],  # 电价变化收敛阈值（比电解铝用能变化更严格）
                verbose=CONFIG["verbose"]  # 静默模式，减少输出
            )
            
            if n is None:
                print("迭代优化失败，跳过该参数组合")
                continue
            
            # # 在优化后调用启动分析
            # if 'smelter' in n.links.index:
            #     analyze_startups(n)
            # else:
            #     print("警告：网络中没有冶炼设备，跳过启动分析")
            
            # 分析排放情况
            scenario_name = f"Iterative Min Power {p_min_pu*100:.0f}%"
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
            
            # 获取目标值
            objective_value = n.objective if hasattr(n, 'objective') else 0.0
            
            # 保存结果
            year_results[p_min_pu] = {
                'p_min_pu': p_min_pu,
                'system_cost': CONFIG["original_cost"] - objective_value * 1e-9, # Billion
                'generator_capacities': n.generators.p_nom_opt * 1e-3, # GW
                'storage_capacities': n.storage_units.p_nom_opt * 1e-3, # GW
                'al_capacity': CONFIG["al_demand"] * (1 + CONFIG["al_excess_rate"]),
                'co2_emissions': co2_emissions,  # 添加CO2排放结果
                'load_shedding': load_shedding_data,  # 添加切负荷数据
                'aluminum_usage': aluminum_usage,  # 添加电解铝用能模式
                'total_time': time_stats['total_time'],
                'iteration_count': time_stats['iteration_count'],
                'iteration_times': time_stats['iteration_times'],
                'avg_iteration_time': time_stats['avg_iteration_time']
            }
            
            # 绘制该过剩率下的用电情况
            # 检查网络中是否有冶炼设备
            if 'smelter' in n.links.index:
                plot_results(n, p_min_pu)
            else:
                print("警告：网络中没有冶炼设备，跳过标准绘图")
                # 使用我们自己的绘图函数，传入迭代次数（这里假设是最后一次迭代）
                plot_iterative_results(n, aluminum_usage, ts, p_min_pu, 1)  # 使用1作为默认迭代次数
            
            print(f"迭代优化完成，系统成本: {year_results[p_min_pu]['system_cost']:.2f} Billion")
        
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
    parser = argparse.ArgumentParser(description='Capacity expansion planning with aluminum smelter iterative optimization')
    parser.add_argument('--config', type=str, default="config.yaml", help='Path to config file (default: config.yaml)')
    args = parser.parse_args()
    
    main(config_file=args.config) 