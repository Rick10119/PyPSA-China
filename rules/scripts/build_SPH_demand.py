import pandas as pd

# 读取人口数据来计算各省的供暖需求
with pd.HDFStore("data/population/population.h5", mode='r') as store:
    population = store['population']
    
# 打印可用的省份名称，以便检查
print("可用的省份名称:", population.index.tolist())

# 计算每个省份的供暖需求
def calculate_heating_demand(population_data):
    # 总供暖需求分配
    total_sph = 777  # TWh (分散式供暖)
    
    # 秦淮线以北省份（权重1.0）
    northern_provinces = [
        'Heilongjiang', 'Jilin', 'Liaoning',  # 东北三省
        'Beijing', 'Tianjin', 'Hebei', 'Shanxi', 'Shandong',  # 华北地区
        'InnerMongolia',  # 内蒙古
        'Gansu', 'Qinghai', 'Ningxia', 'Xinjiang',  # 西北地区
        'Shaanxi', 'Henan', 'Tibet'  # 其他北方省份和西藏
    ]
    
    # 部分供暖省份（权重0.2）
    partial_heating_provinces = [
        'Sichuan', 'Hubei', 'Hunan', 'Jiangxi',  # 中部省份
        'Jiangsu', 'Shanghai', 'Zhejiang'  # 东部沿海省份
    ]
    
    # 打印匹配检查
    print("\n检查省份匹配情况:")
    print("北方省份:")
    for province in northern_provinces:
        exists = province in population_data.index
        print(f"{province}: {'存在' if exists else '不存在'}")
    print("\n部分供暖省份:")
    for province in partial_heating_provinces:
        exists = province in population_data.index
        print(f"{province}: {'存在' if exists else '不存在'}")
    
    # 计算加权人口总和
    weighted_population = (
        population_data.loc[northern_provinces].sum() +  # 权重1.0
        population_data.loc[partial_heating_provinces].sum() * 0.2  # 权重0.2
    )
    
    sph_demand = pd.Series(index=population_data.index, dtype=float)
    
    for province in population_data.index:
        if province in northern_provinces:
            # 北方省份获得全额供暖需求
            sph_demand[province] = (population_data[province] / weighted_population) * total_sph
        elif province in partial_heating_provinces:
            # 部分供暖省份获得0.2权重的供暖需求
            sph_demand[province] = (population_data[province] * 0.2 / weighted_population) * total_sph
        else:
            # 其他省份设为0
            sph_demand[province] = 0
    
    return sph_demand

# 计算供暖需求
sph_demand = calculate_heating_demand(population)

# 保存为CSV
sph_demand.to_csv('data/heating/SPH_2020.csv')

print("\n已生成 SPH_2020.csv 文件")
print("总供暖需求:", sph_demand.sum(), "TWh")
print("\n各省份供暖需求:")
for province in sph_demand.index:
    if sph_demand[province] > 0:
        print(f"{province}: {sph_demand[province]:.2f} TWh")