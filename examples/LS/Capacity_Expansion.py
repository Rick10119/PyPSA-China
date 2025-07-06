# simple_planning_problem.py

from calendar import c
import pandas as pd
import numpy as np
import pypsa
import os
from plot_results import plot_results
from config import CONFIG  # 导入配置

# =====================
# 时间序列数据处理
# =====================
def process_time_series():
    """
    处理时间序列数据。
    读取风电、光伏等可再生能源的出力时间序列。
    """
    # 读入安徽风电与光伏数据

    ts = pd.read_csv("examples/data/anhui_renewable_factors_SingleNode_2020.csv", index_col=0, parse_dates=True)
    # 取前三个月（即第一季度）为优化周期
    ts = ts.iloc[:24*91]  # 取前三个月的数据（假设数据是按小时排列的）

    return ts


# =====================
# 网络创建与建模
# =====================
def create_network(ts):
    """
    创建和配置PyPSA网络。
    包括：
    - 添加电力母线
    - 添加电力负荷
    - 添加候选发电机组
    - 添加候选储能单元
    - 设置目标函数（由PyPSA根据成本自动生成）
    """
    n = pypsa.Network()
    n.set_snapshots(ts.index)
    
    # 添加一个全国性的电力母线（"铜板"模型）
    n.add("Bus", "electricity_bus", carrier="AC")
    
    # 添加电力负荷(弹性/刚性负荷)
    # add_load(n, ts)
    add_FlexibleLoad(n,ts)
    
    # 添加候选发电机组（容量是可优化的）
    add_candidate_generators(n, ts)
    
    # 添加候选储能单元（容量是可优化的）
    add_storage_options(n)
    
    return n

# =====================
# 添加电力负荷
# =====================
# 这个地方是关键，工业生产如果定义为刚性负荷，那么就需要在模型中添加一个固定的电力负荷；如果定义为弹性负荷，则替换为我们的模型。
def add_load(n, ts):
    """
    向网络中添加固定的电力负荷。
    """
    n.add("Load",
          "national_load",
          bus="electricity_bus",
          p_set=ts["load"]) # p_set指定了每个时刻的负荷需求

# =====================
# 添加柔性负荷（产品、仓储和生产过程）
# =====================
def add_FlexibleLoad(n, ts):
    """
    向网络中添加柔性负荷（产品、仓储和生产过程）。
    """
    # 获取配置中的产品信息
    materials = CONFIG["materials"]
    initial_inventories = CONFIG["initial_inventories"]
    relationship_matrix = CONFIG["relationship_matrix"]
    energy_consumption = CONFIG["energy_consumption"]
    demand_matrix = CONFIG["demand_matrix"]
    price_matrix = CONFIG["price_matrix"]

    # 第一部分：添加产品仓储设施
    
    # 添加产品节点（每种物料一个bus，对应文档bus_i）
    for material in materials:
        n.add("Bus", material)

    for material in materials:
        n.add(
            "Store",
            f"{material}_store",
            bus=material,
            carrier="product",  # 仓储类型
            e_nom=10e8,  # 仓储容量上限（可调）
            e_initial=initial_inventories[material],  # 初始库存
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
                # 电力输入 - 设置为每小时平均产能
                "bus0": "electricity_bus",  # 输入中间转换节点
                "carrier": "conversion",  # 转换过程
                "p_max_pu": 1.0,  # 最大出力比例
                "p_nom": 1.2*energy_consumption[j] * CONFIG["demand_matrix"][j] / 24,  # 生产能力上限,假设产能过剩20%
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
    # 构造需求时间序列：每个产品每天最后一小时有需求
    days = len(ts.index) // 24  # Calculate actual number of days from ts dimensions
    demand_series = np.zeros(days*24)
    for i in range(days):
        demand_series[i * 24 + 23] = -1  # 负号表示消耗
    # 为每个产品添加带惩罚成本的负荷
    for i, material in enumerate(materials):
        n.add("Generator",
            f"{material} industrial load",
            bus = material,
            p_max_pu = 0,  # 最大出力为0，只能切负荷
            p_min_pu = demand_series,  # 需求量
            p_nom = demand_matrix[i],  # 需求规模
            p_nom_extendable = False,  # 允许模型优化此负荷的规模
            marginal_cost = price_matrix[i]  # 切负荷成本（万元）
      )

# =====================
# 添加候选发电机组
# =====================
def add_candidate_generators(n, ts):
    """
    添加不同技术的候选发电机组。
    这些是模型需要做投资决策的对象。
    - p_nom_extendable=True 是核心，表示容量可被优化器决定。
    - capital_cost 是投资成本。
    - marginal_cost 是运营成本。
    """
    # 添加风电和光伏发电机组 (可再生能源)
    for tech in ["solar", "onwind"]:
        # 从配置中获取成本信息
        costs = CONFIG["costs"][tech]
        
        # 可再生能源需要指定其可用率曲线
        p_max_pu_series = ts[tech]
        
        n.add("Generator",
              f"{tech}_plant",
              bus="electricity_bus",
              carrier=tech,
              p_nom=CONFIG[f"p_nom_{tech}"],  # 初始容量
              p_nom_min=CONFIG[f"p_nom_{tech}"],  # 当前容量即为最小容量
              p_nom_extendable=True,  # 允许模型优化此发电机的容量
              capital_cost=costs["capital_cost"],  # 投资成本 [万元/MW]
              marginal_cost=costs["marginal_cost"],  # 边际成本 [万元/MWh]
              p_max_pu=p_max_pu_series,  # 最大出力比例（时间序列）
              p_min_pu=0.0  # 最小出力比例为0，允许完全关闭
              )
        
    # 添加火电机组（各种不同类型）
    for coal_type in CONFIG["coal_types"]:
        # 获取机组参数
        unit_data = CONFIG["thermal_units"][coal_type]
        costs = CONFIG["costs"]["thermal"][coal_type]
        
        # 为每台机组单独定义
        unit_count = unit_data["count"]
        for i in range(unit_count):
            n.add("Generator",
            f"{coal_type}_unit_{i}",
            bus="electricity_bus",
            carrier="thermal", 
            p_nom=unit_data["capacity"],  # 单机容量
            p_nom_extendable=False,  # 允许容量扩展
            p_nom_min=unit_data["capacity"],  # 当前容量即为最小容量
            capital_cost=costs["capital_cost"],  # 投资成本
            marginal_cost=costs["marginal_cost"],  # 边际成本
            committtable=True,  # 可启停
            start_up_cost=CONFIG["costs"]["thermal"]["startup_cost"],  # 启动成本
            shutdown_cost=CONFIG["costs"]["thermal"]["shutdown_cost"],  # 停机成本
            p_min_pu=0, #unit_data["p_min"],  # 最小出力比例
            ramp_limit_up=unit_data["ramp_rate"] / unit_data["capacity"],  # 单机爬坡率限制
            ramp_limit_down=unit_data["ramp_rate"] / unit_data["capacity"],  # 单机爬坡率限制
            min_up_time=CONFIG["thermal_units"]["min_up_time"],  # 最小运行时间
            min_down_time=CONFIG["thermal_units"]["min_down_time"]  # 最小停机时间
            )
    
    # 添加水电机组（各种不同类型）
    for hydro_type in CONFIG["hydro_types"]:
        # 获取机组参数
        unit_data = CONFIG["hydro_units"][hydro_type]
        costs = CONFIG["costs"]["hydro"][hydro_type]
        
        # 将n台机组聚合为一台大机组
        n.add("Generator",
             f"{hydro_type}_aggregated",
             bus="electricity_bus", 
             carrier="hydro",
             p_nom=unit_data["capacity"]*unit_data["count"],  # 总容量
             p_nom_extendable=True,  # 允许容量扩展
             p_nom_min=unit_data["capacity"]*unit_data["count"],  # 当前容量即为最小容量
            #  capital_cost=costs["capital_cost"],  # 投资成本
             marginal_cost=costs["marginal_cost"],  # 边际成本
             start_up_cost=CONFIG["costs"]["hydro"]["startup_cost"],  # 启动成本
             shutdown_cost=CONFIG["costs"]["hydro"]["shutdown_cost"],  # 停机成本
             p_min_pu=0, #unit_data["p_min"],  # 最小出力比例
             ramp_limit_up=unit_data["ramp_rate"] / unit_data["capacity"] * unit_data["count"],  # 放宽爬坡率限制
             ramp_limit_down=unit_data["ramp_rate"] / unit_data["capacity"] * unit_data["count"]  # 放宽爬坡率限制
             )

# =====================
# 添加候选储能单元
# =====================
def add_storage_options(n):
    """
    添加候选的储能单元。
    根据CONFIG["storage_params"]中的矩阵格式定义，
    为每种储能类型添加相应数量的储能设备。
    所有储能设备的能量容量(e_nom)可优化，功率容量不进行优化。
    """
    # 获取储能参数
    storage_params = CONFIG["storage_params"]
    storage_types = CONFIG["storage_types"]
    
    # 全局参数
    eff_dispatch = storage_params["efficiency_dispatch"]  # 放电效率
    eff_store = storage_params["efficiency_store"]        # 充电效率
    soc_max = storage_params["soc_max"]                   # SOC上限
    soc_min = storage_params["soc_min"]                   # SOC下限
    special_eff = storage_params["special_efficiency"]    # 特殊效率
    
    # 获取矩阵数据
    matrix = storage_params["matrix"]
    
    # 为每种储能类型添加设备
    for i, storage_type in enumerate(storage_types):
        # 获取该类型的参数
        capacity = matrix["capacity"][i]
        p_charge_max = matrix["p_charge_max"][i]
        p_discharge_max = matrix["p_discharge_max"][i]
        degr_cost = matrix["degradation_cost"][i]
        count = matrix["count"][i]
        eff_override = matrix["eff_override"][i]
        
        # 如果需要使用特殊效率
        curr_eff_dispatch = eff_dispatch
        curr_eff_store = eff_store
        if eff_override:
            curr_eff_dispatch = special_eff["efficiency_dispatch"]
            curr_eff_store = special_eff["efficiency_store"]
        
        # 添加一个聚合后的储能单元
        
        # 添加储能bus
        n.add("Bus",
              f"{storage_type}_bus",
              carrier="battery")

        # 添加储能容量
        n.add("Store",
             f"{storage_type}",
             bus=f"{storage_type}_bus",                  # 连接到储能专用bus
             carrier="battery",
             e_nom=capacity*count,                       # 初始能量容量 
             e_nom_extendable=True,                      # 能量容量可扩展
             e_nom_min=capacity*count,                # 当前能量容量即为最小容量
             e_cyclic=True,                              # 周期性能量约束
             e_min_pu=soc_min,                           # SOC下限
             e_max_pu=soc_max                            # SOC上限
            )

        # 添加充电变流器
        n.add("Link",
             f"{storage_type}_charger",
             bus0="electricity_bus",
             bus1=f"{storage_type}_bus",                 # 连接到储能专用bus
             carrier ="battery",
             p_nom=p_charge_max*count,                   # 固定的功率容量
             efficiency=curr_eff_store,                  # 充电效率
             marginal_cost=degr_cost,                    # 退化成本当做边际成本
            )

        # 添加放电变流器 
        n.add("Link", 
             f"{storage_type}_discharger",
             bus0=f"{storage_type}_bus",                 # 连接到储能专用bus
             bus1="electricity_bus",
             carrier ="battery",
             p_nom=p_discharge_max*count,                # 固定的功率容量
             efficiency=curr_eff_dispatch,               # 放电效率
             marginal_cost=degr_cost,                    # 退化成本当做边际成本
            )

# =====================
# 主运行函数
# =====================
def main():
    """
    主运行函数。
    步骤：
    1. 准备时间序列数据
    2. 创建并配置网络模型
    3. 运行优化求解
    4. 打印并绘制结果
    """
    print("Step 1: Preparing time-series data...")
    ts = process_time_series()
    
    print("Step 2: Creating the network model...")
    network = create_network(ts)
    
    print("Step 3: Solving the optimization problem...")
    network.optimize(solver_name="gurobi")
    
    # print("\n" + "="*20 + " RESULTS " + "="*20)
    # print(f"Total System Cost: {network.objective / 1e6:,.2f} M€")
    # print("\nOptimal Capacities (MW):")
    # print(network.generators.p_nom_opt)
    # if not network.storage_units.empty:
    #     print("\nOptimal Storage Capacities:")
    #     print(f"  Power (MW): {network.storage_units.p_nom_opt.iloc[0]:.2f}")
    #     print(f"  Energy (MWh): {network.storage_units.e_nom_opt.iloc[0]:.2f}")
    # print("="*49 + "\n")
    
    print("Step 4: Plotting results...")
    plot_results(network)

    print(f"Done! Results saved to '{os.path.abspath('planning_results')}' directory.")

# =====================
# 程序入口
# =====================
if __name__ == "__main__":
    main()