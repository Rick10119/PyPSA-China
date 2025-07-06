
import matplotlib.pyplot as plt
import numpy as np
from config import CONFIG

def plot_results(n, products, scenario=None):
    """
    绘制结果图表
    对于不同的场景，绘制不同的图表。
    """

    materials = products['materials']
    hierarchy = products['hierarchy']
    demand = products['demand']
    price = products['price']
    net_demand_rate = products['net_demand_rate']

    # 共有图像
    # 1.绘制发电机功率堆叠图（全周期）
    plt.figure(figsize=(16, 10))
    found_generators = False
    
    # 获取所有发电机类型（除了工业负荷相关的）
    generator_types = ['onwind', 'offwind', 'solar', 'OCGT']
    generation_data = {}
    
    # 设置中文名称映射
    gen_name_map = {
        'onwind': '陆上风电',
        'offwind': '海上风电', 
        'solar': '太阳能',
        'OCGT': '燃气轮机'
    }
    
    for gen_type in generator_types:
        if gen_type in n.generators.index:
            # 检查装机容量是否不为0
            p_nom_opt = n.generators.at[gen_type, 'p_nom_opt']
            if p_nom_opt > 1e-6:  # 避免数值误差，设置小阈值
                # 获取发电机出力时间序列（全周期）
                gen_output = n.generators_t.p[gen_type]
                
                gen_name_cn = gen_name_map.get(gen_type, gen_type)
                generation_data[gen_name_cn] = gen_output.values  # 实际功率值（MW）
                found_generators = True

    # 添加储能放电到发电堆叠中
    storage_types = ['battery storage', 'hydrogen storage underground']
    storage_name_map = {
        'battery storage': '电池储能放电',
        'hydrogen storage underground': '氢储能放电'
    }

    for storage_type in storage_types:
        if storage_type in n.storage_units.index:
            p_nom_opt = n.storage_units.at[storage_type, 'p_nom_opt']
            if p_nom_opt > 1e-6:
                # 获取储能功率（全周期）
                storage_power = n.storage_units_t.p[storage_type]
                # 只取放电部分（正值）
                discharge_power = storage_power.where(storage_power > 0, 0).values
                
                # 如果有显著的放电，则添加到堆叠中
                if discharge_power.sum() > 1e-3:
                    storage_name_cn = storage_name_map.get(storage_type, storage_type)
                    generation_data[storage_name_cn] = discharge_power
                    found_generators = True

    if found_generators and generation_data:
        # 准备堆叠数据
        time_index = n.generators_t.p.index
        labels = list(generation_data.keys())
        data_matrix = np.array([generation_data[label] for label in labels])
        
        # 设置颜色
        colors = ['dodgerblue', 'aquamarine', 'gold', 'indianred', 'purple', 'orange', 'green', 'pink']
        
        # 绘制堆叠面积图
        plt.stackplot(time_index, *data_matrix, labels=labels, 
                     colors=colors[:len(labels)], alpha=0.8)
        
        # 添加总负荷线（如果是scenario 1）
        if scenario == 1:
            # 计算总负荷
            total_load = sum(products["energy_consumption"][i] * products["demand"][i] 
                            for i in range(len(products["materials"]))) * 366 / (365.25 * 24)  # 平均功率
            total_load_line = np.full(len(time_index), total_load)
            plt.plot(time_index, total_load_line, 'r--', linewidth=3, label='总负荷需求')
        
        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.title('发电机功率堆叠图（全周期）')
        plt.xlabel('时间')
        plt.ylabel('功率 (MW)')
        plt.grid(True, alpha=0.3)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        if scenario is not None:
            plt.savefig(f"examples/scripts_LS/results/power_stack_full_scenario{scenario}.png", bbox_inches='tight')
        else:
            plt.savefig("examples/scripts_LS/results/power_stack_full.png", bbox_inches='tight')

    plt.close()

    # 2.绘制已建储能系统的充放电情况
    plt.figure(figsize=(14, 8))
    found_storage = False
    
    # 获取所有储能类型
    storage_types = ['battery storage', 'hydrogen storage underground']
    
    for storage_type in storage_types:
        if storage_type in n.storage_units.index:
            # 检查装机容量是否不为0
            p_nom_opt = n.storage_units.at[storage_type, 'p_nom_opt']
            if p_nom_opt > 1e-6:  # 避免数值误差，设置小阈值
                # 获取储能充放电功率时间序列
                storage_power = n.storage_units_t.p[storage_type]
                
                # 设置中文名称映射
                storage_name_map = {
                    'battery storage': '电池储能',
                    'hydrogen storage underground': '地下氢储能'
                }
                
                storage_name_cn = storage_name_map.get(storage_type, storage_type)
                plt.plot(storage_power.index, storage_power, 
                        label=f'{storage_name_cn} ({p_nom_opt:.1f} GW)', linewidth=1.5)
                found_storage = True
    
    if found_storage:
        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.title('储能系统充放电情况')
        plt.xlabel('时间')
        plt.ylabel('充放电功率 (MW)')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.xticks(rotation=45)
        plt.axhline(y=0, color='black', linestyle='--', alpha=0.5)  # 添加零线
        plt.tight_layout()
        if scenario is not None:
            plt.savefig(f"examples/scripts_LS/results/storage_charging_discharging_scenario{scenario}.png")
        else:
            plt.savefig("examples/scripts_LS/results/storage_charging_discharging.png")
    plt.close()

    # 3.绘制储能系统荷电状态（SOC）
    plt.figure(figsize=(14, 8))
    found_soc = False
    
    for storage_type in storage_types:
        if storage_type in n.storage_units.index:
            # 检查装机容量是否不为0
            p_nom_opt = n.storage_units.at[storage_type, 'p_nom_opt']
            if p_nom_opt > 1e-6:
                # 获取储能荷电状态
                if storage_type in n.storage_units_t.state_of_charge.columns:
                    soc = n.storage_units_t.state_of_charge[storage_type]
                    
                    # 获取储能容量
                    max_hours = n.storage_units.at[storage_type, 'max_hours']
                    energy_capacity = p_nom_opt * max_hours  # MWh
                    
                    # 标幺化SOC
                    normalized_soc = soc / energy_capacity if energy_capacity > 0 else soc
                    
                    storage_name_cn = storage_name_map.get(storage_type, storage_type)
                    plt.plot(normalized_soc.index, normalized_soc, 
                            label=f'{storage_name_cn} ({energy_capacity:.1f} GWh)', linewidth=1.5)
                    found_soc = True
    
    if found_soc:
        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.title('储能系统荷电状态（标幺值）')
        plt.xlabel('时间')
        plt.ylabel('荷电状态 (p.u.)')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.xticks(rotation=45)
        plt.ylim(0, 1.1)  # 设置y轴范围
        plt.tight_layout()
        if scenario is not None:
            plt.savefig(f"examples/scripts_LS/results/storage_state_of_charge_scenario{scenario}.png")
        else:
            plt.savefig("examples/scripts_LS/results/storage_state_of_charge.png")
    plt.close()

    # 各Scenario独有图像
    if scenario is None:
        pass
    elif scenario == 1:
        pass
    elif scenario == 2:
        pass
    elif scenario == 3:
        # 按产品类型（type）绘制每日电力消耗图
        product_types = CONFIG["type"]  # 从config获取所有产品类型
    
        for product_type in product_types:
            plt.figure(figsize=(12, 6))
            found_products = False
            
            for i, material in enumerate(materials):
                # 检查该产品是否属于当前类型
                if i < len(products["type"]) and products["type"][i] == product_type:
                    link_name = f"{material}_production"
                    if (link_name in n.links.index and 
                        link_name in n.links_t.p0.columns):
                        
                        # 获取Link的功率时间序列（MW）
                        power_series = n.links_t.p0[link_name]
                        
                        # 转换为每日能耗（MWh）
                        # 方法：将功率乘以时间分辨率，然后按天求和
                        daily_energy = (power_series * CONFIG["resolution"]).resample('D').sum()
                        
                        # 标幺化
                        daily_energy = daily_energy / (products["energy_consumption"][i] * products["demand"][i] )

                        # 绘制该产品的每日能耗
                        plt.plot(daily_energy.index, daily_energy, 
                                label=material, linewidth=1.5)
                        found_products = True
            
            # 如果该类型有产品，则保存图表
            if found_products:
                plt.rcParams['font.sans-serif'] = ['SimHei']
                
                # 设置中文类型名称映射
                type_name_map = {
                    "mineral": "矿物原料",
                    "chemical": "化工产品", 
                    "building": "建筑材料",
                    "vehicle": "交通工具",
                    "component": "零部件",
                    "fertilizer": "化肥产品",
                    "electrical_equipment": "电气设备",
                    "industrial_equipment": "工业设备"
                }
                
                type_name_cn = type_name_map.get(product_type, product_type)
                plt.title(f'{type_name_cn}每日电力消耗')
                plt.xlabel('日期')
                plt.ylabel('每日总电力消耗 (MWh)')
                plt.grid(True, alpha=0.3)
                plt.legend()
                plt.xticks(rotation=45)
                plt.tight_layout()
                plt.savefig(f"examples/scripts_LS/results/{product_type}_daily_electricity_usage.png")
            plt.close()
        
        # 新增1: 按产品类型绘制切负荷情况图 - 每天最后一个时段的标幺化出力
        for product_type in product_types:
            plt.figure(figsize=(14, 8))
            found_products = False
            
            # 获取时间序列信息
            AllTime = len(n.snapshots)
            total_days = AllTime * CONFIG["resolution"] // 24
            
            # 为每个产品收集每天最后时段的标幺化供给
            for i, material in enumerate(materials):
                if i < len(products["type"]) and products["type"][i] == product_type:
                    gen_name = f"{material} industrial load"
                    if gen_name in n.generators.index:
                        
                        # 获取发电机出力和额定功率
                        gen_output = n.generators_t.p[gen_name]
                        p_nom = n.generators.at[gen_name, "p_nom"]

                        # 提取每天最后一个时段的供给并标幺化
                        daily_last_output = []
                        dates = []
                        
                        for day in range(total_days):
                            last_hour_idx = day * (24 // CONFIG["resolution"]) + (24 // CONFIG["resolution"]) - 1
                            if last_hour_idx < len(gen_output):
                                # 标幺化：除以p_nom，取绝对值
                                normalized_output = abs(gen_output.iloc[last_hour_idx] / p_nom) if p_nom != 0 else 0
                                daily_last_output.append(normalized_output)
                                dates.append(gen_output.index[last_hour_idx].date())
                        
                        if daily_last_output:
                            plt.plot(dates, daily_last_output, 
                                    label=f'{material}', linewidth=1.5)
                            found_products = True
            
            if found_products:
                plt.rcParams['font.sans-serif'] = ['SimHei']
                type_name_cn = type_name_map.get(product_type, product_type)
                plt.title(f'{type_name_cn}产品供给情况（标幺化）- 每日最后时段')
                plt.xlabel('日期')
                plt.ylabel('标幺化产品供给 (p.u.)')
                plt.grid(True, alpha=0.3)
                plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.xticks(rotation=45)
                plt.tight_layout()
                plt.savefig(f"examples/scripts_LS/results/{product_type}_load_shedding_normalized.png", 
                        bbox_inches='tight')
            plt.close()
        
        # 新增2: 产品切负荷成本占比饼图
        # 计算每个产品的切负荷成本
        product_costs = {}
        total_cost = 0

        for i, material in enumerate(materials):
            gen_name = f"{material} industrial load"
            if gen_name in n.generators.index:
                # 计算该产品的切负荷成本
                cost = CONFIG['resolution']*(n.generators_t.p[gen_name] * n.generators.at[gen_name, "marginal_cost"]).sum()
                cost = products['product_total_profit'][material] + cost  # 产品的切负荷成本（欧元）
                abs_cost = abs(cost)  # 取绝对值
                
                if abs_cost > 0:  # 只包含有切负荷成本的产品
                    product_costs[material] = abs_cost
                    total_cost += abs_cost
            
        # 打印详细的成本信息
        print("\n=== 产品切负荷成本详情 ===")
        print(f"{'产品名称':^15} | {'切负荷成本(EUR)':^15} | {'占比(%)':^10}")
        print("-" * 50)
        for material, cost in sorted(product_costs.items(), key=lambda x: x[1], reverse=True):
            percentage = cost/total_cost*100
            print(f"{material:^15} | {cost:^15.2e} | {percentage:^10.2f}")
        print(f"{'总计':^15} | {total_cost:^15.2e} | {100.0:^10.2f}")

        # 按hierarchy分别绘制各个产品的store的存储量折线图
        for carrier in ['upstream_product', 'downstream_product']:
            plt.figure(figsize=(12, 6))
            found_products = False
            
            for i, material in enumerate(materials):
                if hierarchy[i] == carrier:  # 只绘制当前carrier类型的产品
                    store_name = f"{material}_store"
                    if store_name in n.stores.index:
                        store_data = n.stores_t.e[store_name] / n.stores.at[store_name, "e_nom"] *14 # 标幺化存储量
                        plt.plot(store_data.index, store_data, label=material)
                        found_products = True
            
            if found_products:
                plt.rcParams['font.sans-serif'] = ['SimHei']
                carrier_name = "上游产品" if carrier == "upstream_product" else "下游产品"
                plt.title(f'{carrier_name}仓储存储量变化')
                plt.xlabel('时间')
                plt.ylabel('存储量 (标幺化)')
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.xticks(rotation=45)
                plt.tight_layout()
                plt.savefig(f"examples/scripts_LS/results/{carrier}_storage_levels.png")
            plt.close()
        
        # 单独绘制铜材生产的用能（全周期）
            if "copper" in materials:
                plt.figure(figsize=(12, 6))
                gen_name = "copper_production"
                if gen_name in n.links.index:
                # 获取铜材生产Link的出力数据（全周期）
                    production_data = n.links_t.p0[gen_name]
                # 标幺化：除以铜材生产的额定功率
                    production_data = production_data / n.links.at[gen_name, "p_nom"]  # 标幺化生产用电量
                plt.plot(production_data.index, production_data, label='copper production', linewidth=2, color='orange')
                plt.rcParams['font.sans-serif'] = ['SimHei']
                plt.title('铜材生产用电量变化（全周期）')
                plt.xlabel('时间')
                plt.ylabel('用电量 (MW)')  # 用电量单位为MW
                plt.legend()
                plt.grid(True, alpha=0.3)
                plt.xticks(rotation=45)
                plt.tight_layout()
                plt.savefig("examples/scripts_LS/results/copper_production_usage_full_period.png")
                plt.close()

        # # 单独绘制氧化铝的仓储变化图（前三天）
        # if "alumina" in materials:
        #     plt.figure(figsize=(12, 6))
        #     store_name = "alumina_store"
        #     if store_name in n.stores.index:
        #         store_data = abs(n.stores_t.e[store_name]-2e8)   # 获取存储量数据，除以1e8显示较小数值
        #         # 计算前三天的数据点数：3天 * 24小时 / resolution
        #         three_days_points = 14 * 24 // CONFIG["resolution"]
        #         store_data_14days = store_data.iloc[:three_days_points]  # 取前三天数据

        #         plt.plot(store_data_14days.index, store_data_14days, label='alumina', linewidth=2, color='red')

        #         plt.rcParams['font.sans-serif'] = ['SimHei']
        #         plt.title('氧化铝仓储存储量变化（前三天）')
        #         plt.xlabel('时间')
        #         plt.ylabel('存储量')  # 更改单位说明已缩放
                
        #         # 设置y轴格式，确保数字正确显示
        #         plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.2f}'))
                
        #         plt.legend()
        #         plt.grid(True, alpha=0.3)
        #         plt.xticks(rotation=45)
        #         plt.tight_layout()
        #         plt.savefig("examples/scripts_LS/results/alumina_storage_levels_3days.png")
        #     plt.close()
        
        # 单独绘制铜材的仓储变化图（全周期）
        if "copper" in materials:
            plt.figure(figsize=(12, 6))
            store_name = "copper_store"
            if store_name in n.stores.index:
                store_data = n.stores_t.e[store_name]  # 获取存储量数据（全周期）
            # 标幺化
                store_data = store_data / n.stores.at[store_name, "e_nom"] * 14 / products["net_demand_rate"][15]  # 标幺化存储量，净需求为额定值

            plt.plot(store_data.index, store_data, label='copper', linewidth=2, color='blue')

            plt.rcParams['font.sans-serif'] = ['SimHei']
            plt.title('铜材仓储存储量变化（全周期）')
            plt.xlabel('时间')
            plt.ylabel('存储量 (标幺化)')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig("examples/scripts_LS/results/copper_storage_levels_full_period.png")
            plt.close()
