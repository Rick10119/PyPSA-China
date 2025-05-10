# LoadShedding.py

import pandas as pd
import numpy as np
import pypsa
import os
import matplotlib.pyplot as plt
from config import CONFIG  # 导入配置

# =====================
# 时间序列数据处理
# =====================
def process_time_series():
    """
    处理时间序列数据。
    读取风电、光伏等可再生能源的出力时间序列。
    只取四月份（第4个720小时），与实验数据一致。
    """
    ts = pd.read_csv("examples/data/time-series-lecture-2.csv", index_col=0, parse_dates=True)
    # 取四月份数据，即第四个720小时
    ts = ts.iloc[720*3:720*4]
    return ts

# =====================
# 网络创建与建模
# =====================
def create_network(ts):
    """
    创建和配置PyPSA网络。
    包括：
    - 添加电力总线（electricity）
    - 添加各产品节点（与文档中的bus_i对应）
    - 添加发电机组（风电、光伏等）
    - 添加carrier（能源类型）
    - 添加仓储设施（Store，约束s_{i,t}）
    - 添加生产过程的Link（物料守恒、生产约束）
    - 添加弹性负荷（负荷切除，目标函数相关）
    """
    n = pypsa.Network()
    # 添加基础节点（电力总线）
    n.add("Bus", "electricity")
    n.set_snapshots(ts.index)
    # 添加产品节点（每种物料一个bus，对应文档bus_i）
    for material in CONFIG["materials"]:
        n.add("Bus", material)
    # 添加发电机组（风电、光伏等）
    add_generators(n, ts)
    # 添加carriers（能源类型）
    carriers = ["onwind", "offwind", "solar"]
    n.add(
        "Carrier",
        carriers,
        color=["dodgerblue", "aquamarine", "gold"]
    )
    # 添加产品仓储设施（Store，约束s_{i,t}）
    add_storage(n, CONFIG["materials"], CONFIG["initial_inventories"])
    # 添加生产过程的link（物料守恒、生产约束）
    add_links(n, CONFIG["materials"], CONFIG["relationship_matrix"])
    # 添加带有负荷切除成本的弹性负荷（目标函数相关）
    add_load_shedding(n, CONFIG["materials"], CONFIG["demand_matrix"], CONFIG["price_matrix"])
    return n

# =====================
# 添加发电机组
# =====================
def add_generators(n, ts):
    """
    添加可再生能源发电机组。
    每个发电机组的出力固定为ts[tech]，即时间序列数据。
    对应文档中的g_t。
    """
    for tech in ["onwind", "offwind", "solar"]:
        n.add(
            "Generator",
            tech,
            bus="electricity",
            carrier=tech,
            p_min_pu=ts[tech],  # 最小出力（按时序）
            p_max_pu=ts[tech]   # 最大出力（按时序）
            # 可添加成本、效率等参数
        )

# =====================
# 添加仓储设施（Store）
# =====================
def add_storage(n, materials, initial_inventories):
    """
    为每个产品添加Store元件，用于存储。
    对应文档中的s_{i,t}，并施加仓储上限约束。
    e_nom: 仓储容量上限（s_{i,max}）
    e_initial: 初始库存（s_{i,initial}）
    e_cyclic: 是否循环（周期性）
    """
    for material in materials:
        n.add(
            "Store",
            f"{material}_store",
            bus=material,
            e_nom=1000,  # 仓储容量上限（可调）
            e_initial=initial_inventories[material],  # 初始库存
            e_marginal_cost=-0.1,  # 边际成本（可调）
            e_cyclic=True)  # 是否周期性

# =====================
# 添加生产过程的Link
# =====================
def add_links(n, materials, relationship_matrix):
    """
    添加产品之间的生产Link。
    对应文档中的物料守恒约束、生产约束：
    s_{i,t} = s_{i,t-1} - sum_j(eta_{ij}*e_j*p_{j,t}) + e_i*p_{i,t} - m_{i,t}
    其中relationship_matrix[i][j] = eta_{ij}
    """
    for j, material in enumerate(materials):
        inputs = []
        efficiencies = []
        for i, input_material in enumerate(materials):
            if relationship_matrix[i][j] > 0:
                inputs.append(input_material)
                efficiencies.append(relationship_matrix[i][j])
        # 创建Link元件（每个产品一个生产Link）
        if inputs:
            link_kwargs = {
                "bus0": "electricity",  # 电力输入
                # "efficiency": -1.0,  # 可选
                "p_nom": 20  # 生产能力上限（可调）
            }
            # 设置输入物料的效率（消耗量为负）
            link_kwargs[f"efficiency"] = -efficiencies[0]
            for idx, input_material in enumerate(inputs):
                link_kwargs[f"bus{idx + 1}"] = input_material
                if idx + 1 != 1:
                    link_kwargs[f"efficiency{idx + 1}"] = -efficiencies[idx]
            # 添加产品输出bus
            link_kwargs[f"bus{len(inputs) + 1}"] = material
            # 产品输出效率为1
            link_kwargs[f"efficiency{len(inputs) + 1}"] = 1.0
            link_name = f"{material}_production"
            n.add("Link", link_name, **link_kwargs)
    # 打印Link参数，便于debug
    for link_name in n.links.index:
        print(f"Link: {link_name}")
        for column in n.links.columns:
            print(f"  {column}: {n.links.at[link_name, column]}")
    print('test')

# =====================
# 添加弹性负荷（负荷切除）
# =====================
def add_load_shedding(n, materials, demand_matrix, price_matrix):
    """
    添加弹性负荷（负荷切除），对应目标函数：
    min sum_t sum_i c_i * R_{d,i,t}
    其中R_{d,i,t} = d_{i,t} - m_{i,t}
    这里通过Generator的p_min_pu设置需求，p_max_pu=0表示只能切负荷。
    marginal_cost为切负荷成本c_i。
    """
    # 构造需求时间序列：每个产品每天最后一小时有需求
    demand_series = np.zeros(30 * 24)
    for i in range(30):
        demand_series[i * 24 + 23] = -1  # 负号表示消耗
    # 为每个产品添加带惩罚成本的负荷
    for i, material in enumerate(materials):
        n.add("Generator",
            f"{material} industrial load",
            bus = material,
            p_max_pu = 0,  # 最大出力为0，只能切负荷
            p_min_pu = demand_series,  # 需求量
            p_nom = demand_matrix[i],  # 需求规模
            marginal_cost = price_matrix[i])  # 切负荷成本

# =====================
# 结果绘图
# =====================
def plot_results(n):
    """
    绘制结果图像。
    目前只绘制了links的输入功率。
    可扩展为绘制储能、负荷切除等结果。
    """
    n.links_t.p0.plot(figsize=(9, 7), lw=3)
    plt.tight_layout()
    plt.savefig("examples/LS/results/links_p0.png")
    plt.show()

# =====================
# 主运行函数
# =====================
def main():
    """
    主运行函数。
    步骤：
    1. 创建结果保存目录
    2. 处理时间序列数据
    3. 创建并优化网络
    4. 绘制结果
    """
    os.makedirs("LS/results", exist_ok=True)
    ts = process_time_series()
    n = create_network(ts)
    n.optimize()  # 求解器优化（默认Glpk）
    plot_results(n)

# 运行一次main函数
if __name__ == "__main__":
    main()

