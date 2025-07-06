# LoadShedding.py
"""
 使用安徽省2020年真实装机数据，创建一个简单的单节点模型
 认为安徽省的电力系统是一个铜板模型，只有一个电力总线。
 安徽省产品需求按照2020年安徽省GDP占比从全国需求按比例分配。

"""
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
    """
    # 读入安徽省2020年风电与光伏数据
    ts = pd.read_csv("examples/data/anhui_renewable_factors_SingleNode_2020.csv", index_col=0, parse_dates=True)
    # 取新能源出力最低的连续7天，11-23至11-29
    ts = ts.loc["2020-11-23 00:00":"2020-11-29 23:00"]

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
    carriers = ["onwind", "solar"]
    n.add(
        "Carrier",
        carriers,
        color=["dodgerblue", "gold"]
    )
    # 添加产品仓储设施（Store，约束s_{i,t}）
    add_storage(n, CONFIG["materials"], CONFIG["initial_inventories"])
    # 添加生产过程的link（物料守恒、生产约束）
    add_links(n, CONFIG["materials"], CONFIG["relationship_matrix"], CONFIG["energy_consumption"])
    # 添加带有负荷切除成本的弹性负荷（目标函数相关）
    # add_load_shedding(n, CONFIG["materials"], CONFIG["demand_matrix"], CONFIG["price_matrix"])
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
    for tech in ["onwind", "solar"]:
        n.add(
            "Generator",
            tech,
            bus="electricity",
            carrier=tech,
            p_min_pu=0,  # 最小出力（按时序）
            p_max_pu=ts[tech],   # 最大出力（按时序）
            p_nom=CONFIG[f"p_nom_{tech}"]*10, # 发电机组容量（可调）
            # p_nom_extendable=True  # 容量可扩展
            # 可添加成本、效率等参数
            # marginal_cost=CONFIG["costs"][tech]["marginal_cost"]  # 边际成本（可调）
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
            e_nom=10e8,  # 仓储容量上限（可调）
            e_initial=initial_inventories[material],  # 初始库存
            # e_marginal_cost=-0.1,  # 边际成本（可调）
            e_cyclic=False)  # 是否周期性

# =====================
# 添加生产过程的Link
# =====================
def add_links(n, materials, relationship_matrix, energy_consumption):
    """
    添加产品之间的生产Link。
    对应文档中的物料守恒约束、生产约束：
    s_{i,t} = s_{i,t-1} - sum_j(eta_{ij}*e_j*p_{j,t}) + e_i*p_{i,t} - m_{i,t}
    其中relationship_matrix[i][j] = eta_{ij}
    """
    # 使用单位用能向量修正关系矩阵，单位向量对应元素除以矩阵的对应列
    relationship_matrix = np.array(relationship_matrix) / np.array(energy_consumption).reshape(-1, 1)
    # 打印关系矩阵，便于debug
    # print("Relationship Matrix:")
    # print(relationship_matrix)

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
                "bus0": "electricity",  # 输入电能
                "p_nom": 100, # 生产能力上限（可调）
                "p_max_pu": 1.0,  # 最大出力比例
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

    # # 打印Link参数，便于debug
    # for link_name in n.links.index:
    #     print(f"Link: {link_name}")
    #     for column in n.links.columns:
    #         print(f"  {column}: {n.links.at[link_name, column]}")
    # print('test')

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
    demand_series = np.zeros(7 * 24)
    for i in range(7):
        demand_series[i * 24 + 23] = -1  # 负号表示消耗
    # 为每个产品添加带惩罚成本的负荷
    for i, material in enumerate(materials):
        n.add("Generator",
            f"{material} industrial load",
            bus = material,
            p_max_pu = 0,  # 最大出力为0，只能切负荷
            p_min_pu = demand_series,  # 需求量
            p_nom = demand_matrix[i],  # 需求规模
            marginal_cost = price_matrix[i])  # 切负荷成本（万元）

# =====================
# 产品分类定义
# =====================
def categorize_materials(materials):
    """
    将产品分为四类：
    1. Upstream raw materials - 上游原材料
    2. Primary products - 初级产品
    3. Downstream raw materials - 下游原材料
    4. Downstream products - 下游产品
    
    返回一个字典，包含四个分类的产品列表
    """
    categories = {
        "Upstream raw materials": [],
        "Primary products": [],
        "Downstream raw materials": [],
        "Downstream products": []
    }
    
    # 根据config.py中的分类来分组产品
    for i, material in enumerate(materials):
        if i < 8:  # 前8个是上游原材料
            categories["Upstream raw materials"].append(material)
        elif i < 13:  # 接下来5个是初级产品
            categories["Primary products"].append(material)
        elif i < 14:  # 接下来1个是下游原材料
            categories["Downstream raw materials"].append(material)
        else:  # 剩余的是下游产品
            categories["Downstream products"].append(material)
    
    return categories

# 在绘图函数中使用分类信息
def plot_results(n):
    """
    绘制结果图像，按产品分类进行绘制。
    """
    # 获取产品分类
    categories = categorize_materials(CONFIG["materials"])
    
    # 创建保存结果的目录
    os.makedirs("examples/LS/results", exist_ok=True)

    # 绘制发电机输出
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 设置中文
    plt.figure(figsize=(12, 6))
    renewable_gens = ['onwind', 'solar']
    n.generators_t.p[renewable_gens].plot(lw=2)
    plt.title("可再生能源发电量")
    plt.ylabel("发电量 (MW)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("examples/LS/results/renewable_generators.png")
    plt.close()
    
    # 打印目标函数值
    print(f"目标函数值: {n.objective}")
    # 计算总增加值（商品需求量* 商品价格*天数）
    total_value = sum(CONFIG["demand_matrix"][i] * CONFIG["price_matrix"][i] * 7 for i, material in enumerate(CONFIG["materials"]))
    print(f"总增加值: {total_value} 万元")

    # # 按分类绘制产品库存
    # for category, materials_list in categories.items():
    #     if not materials_list:
    #         continue
            
    #     plt.figure(figsize=(12, 6))
    #     store_names = [f"{material}_store" for material in materials_list]
    #     data = n.stores_t.e[store_names]
        
        # # 对大数值进行缩放以便更好显示
        # for col in data.columns:
        #     if data[col].max() > 1e6:
        #         scale = 1e6
        #         data[col] = data[col] / scale
        #         col_parts = col.split('_')
        #         new_col = f"{col_parts[0]} (百万吨)"
        #         data.rename(columns={col: new_col}, inplace=True)
        
        # data.plot(lw=2)
        # plt.title(f"{category} - 库存变化")
        # plt.ylabel("库存量")
        # plt.grid(True)
        # plt.tight_layout()
        # plt.savefig(f"examples/LS/results/{category.lower().replace(' ', '_')}_inventory.png")
        # plt.show()
    
    # # 按分类绘制产品生产用能
    # for category, materials_list in categories.items():
    #     if category in ["Upstream raw materials", "Downstream raw materials"]:
    #         continue  # 上下游原材料都没有生产过程
            
    #     plt.figure(figsize=(12, 6))
    #     electricity_links = [f"{material}_electricity_input" for material in materials_list]
    #     electricity_links = [name for name in electricity_links if name in n.links_t.p0.columns]
        
    #     if not electricity_links:
    #         continue
            
    #     n.links_t.p0[electricity_links].plot(lw=2)
    #     plt.title(f"{category} - Production Energy Consumption")
    #     plt.ylabel("Power (MW)") 
    #     plt.grid(True)
    #     plt.tight_layout()
    #     plt.savefig(f"examples/LS/results/{category.lower().replace(' ', '_')}_production_energy.png")
    #     plt.show()
    
    # # 绘制负荷切除情况
    # plt.figure(figsize=(12, 6))
    # # 获取所有电力转换Link的输入功率，用于绘制工业生产用能
    # electricity_links = [link for link in n.links.index if 'electricity_input' in link]
    # electricity_consumption = n.links_t.p0[electricity_links]
    # # 按分类绘制工业生产用能
    # for category, materials_list in categories.items():
    #     category_links = [f"{material}_electricity_input" for material in materials_list]
    #     category_links = [link for link in category_links if link in electricity_links]
        
    #     if not category_links:
    #         continue
            
    #     plt.figure(figsize=(12, 6))
    #     electricity_consumption[category_links].plot(lw=2)
    #     plt.title(f"{category} - Industrial Energy Consumption")
    #     plt.ylabel("Power (MW)")
    #     plt.grid(True)
    #     plt.tight_layout()
    #     plt.savefig(f"examples/LS/results/{category.lower().replace(' ', '_')}_energy_consumption.png")

        

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
    os.makedirs("examples/LS/results", exist_ok=True)
    ts = process_time_series()
    n = create_network(ts)
    n.optimize()  # 求解器优化（默认Glpk）
    plot_results(n)

# 运行一次main函数
if __name__ == "__main__":
    main()

