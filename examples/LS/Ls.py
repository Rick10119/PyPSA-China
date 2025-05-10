
# LoadShedding.py

import pandas as pd
import numpy as np
import pypsa
import os
import matplotlib.pyplot as plt
from config import CONFIG  # 导入配置

def process_time_series():
    """处理时间序列数据"""
    ts = pd.read_csv("examples/data/time-series-lecture-2.csv", index_col=0, parse_dates=True)
    
    # 取四月份数据，即第四个720小时
    ts = ts.iloc[720*3:720*4]
    
    return ts

def create_network(ts):
    """创建和配置网络"""
    n = pypsa.Network()
    
    # 添加基础节点
    n.add("Bus", "electricity")
    n.set_snapshots(ts.index)
    
    # 添加产品节点
    for material in CONFIG["materials"]:
        n.add("Bus", material)

    # 添加发电机组
    add_generators(n, ts)
    
    # 添加carriers
    carriers = ["onwind", "offwind", "solar"]
    n.add(
        "Carrier",
        carriers,
        color=["dodgerblue", "aquamarine", "gold"]
    )

    # 添加产品仓储设施
    add_storage(n, CONFIG["materials"], CONFIG["initial_inventories"])
    
    # 添加表示生产过程的link
    add_links(n, CONFIG["materials"], CONFIG["relationship_matrix"])

    # 添加带有负荷切除成本的弹性负荷
    add_load_shedding(n, CONFIG["materials"], CONFIG["demand_matrix"], CONFIG["price_matrix"])
    
    return n

def add_generators(n, ts):
    """添加发电设备"""

    # 添加可再生能源发电机,其出力固定为ts[tech]，即时间序列数据
    for tech in ["onwind", "offwind", "solar"]:
        n.add(
            "Generator",
            tech,
            bus="electricity",
            carrier=tech,
            p_min_pu=ts[tech],
            p_max_pu=ts[tech]
            #capital_cost=costs.at[tech, "capital_cost"],
            #marginal_cost=costs.at[tech, "marginal_cost"],
            #efficiency=costs.at[tech, "efficiency"],
        )

def add_storage(n, materials, initial_inventories):
    """添加产品仓储设施"""
    # for material in materials:
    #     n.add(
    #         "StorageUnit",
    #         f"{material}_storage",  
    #         bus=material,
    #         # carrier="electricity",
    #         p_nom=1000,  # 最大功率容量（MW）
    #         p_nom_extendable=False,  # 不允许扩展功率容量
    #         state_of_charge_initial=initial_inventories[material],  # 初始储能状态（MWh）
    #         cyclic_state_of_charge=False,  # 允许循环储能状态
    #     )

    # 为每个产品添加store元件, 用于存储
    for material in materials:
        n.add(
            "Store",
            f"{material}_store",
            bus=material,
            e_nom=1000,
            e_initial=initial_inventories[material],
            e_marginal_cost=-0.1,
            e_cyclic=False)
    
         
def add_links(n, materials, relationship_matrix):
    """添加产品链接"""

   # 添加Link元件
    for j, material in enumerate(materials):
        inputs = []
        efficiencies = []
        for i, input_material in enumerate(materials):
            if relationship_matrix[i][j] > 0:
                inputs.append(input_material)
                efficiencies.append(relationship_matrix[i][j])

        # 创建Link元件
        if inputs:
            link_kwargs = {
                "bus0": "electricity",
                # "efficiency": -1.0,
                "p_nom": 20
            }
            
            link_kwargs[f"efficiency"] = -efficiencies[0]
            # 动态添加输入节点和效率
            for idx, input_material in enumerate(inputs):
                link_kwargs[f"bus{idx + 1}"] = input_material
                # 如果idx+1不是1则执行
                if idx + 1 != 1:
                    link_kwargs[f"efficiency{idx + 1}"] = -efficiencies[idx]
            
            # 添加产品的bus
            link_kwargs[f"bus{len(inputs) + 1}"] = material
            # 添加产品的效率
            link_kwargs[f"efficiency{len(inputs) + 1}"] = 1.0

            link_name = f"{material}_production"
            n.add("Link", link_name, **link_kwargs)


    # 打印links的每一行与每一列的参数
    for link_name in n.links.index:
        print(f"Link: {link_name}")
        for column in n.links.columns:
            print(f"  {column}: {n.links.at[link_name, column]}")
    print('test')

def add_load_shedding(n, materials, demand_matrix, price_matrix):
    """添加弹性负荷"""
    
    # 定义每个产品的需求量，30天中，在每天的最后一小时有需求，其他时段无需求
    # 定义一个长度为30*24的序列，代表每个产品在每个时间点的需求量，在每个24的倍数时刻有需求，为1，其他时刻为0
    demand_series = np.zeros(30 * 24)
    for i in range(30):
        demand_series[i * 24 + 23] = -1  # 每天的最后一个小时有需求
    # 这个数组作为负荷的p_min_pu, 代表每个产品在每个时间点的需求量

    # 为每个产品添加带惩罚成本的负荷
    for i, material in enumerate(materials):
        # 添加负荷元件
        n.add("Generator",
            f"{material} industrial load",
            bus = material,
            p_max_pu = 0,
            p_min_pu = demand_series,
            p_nom = demand_matrix[i],
            marginal_cost = price_matrix[i])

# 绘制一些图像
def plot_results(n):
    """绘制结果图像"""

    # 这个函数还没有写好
    # 绘制links的输入功率
    n.links_t.p0.plot(figsize=(9, 7), lw=3)
    plt.tight_layout()
    plt.savefig("examples/LS/results/links_p0.png")
    print('test')

def main():
    """主运行函数"""
    # 创建结果保存目录
    os.makedirs("LS/results", exist_ok=True)
     
    # 数据处理
    # costs = process_cost_data(year)
    ts = process_time_series()
        
    # 创建并优化网络
    n = create_network(ts)

    # 设置求解器,并求解网络
    n.optimize(solver_name="gurobi", solver_options={
        "OutputFlag": 0,  # 不显示求解过程
        "IntFeasTol": 1e-2  # 整数可行性容差
    })


# 运行一次main函数
if __name__ == "__main__":
    main()

