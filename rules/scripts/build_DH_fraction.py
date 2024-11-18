import pandas as pd

# 读取SPH需求数据
sph_demand = pd.read_csv('data/heating/SPH_2020.csv', index_col=0).squeeze()

# 读取人口数据
with pd.HDFStore("data/population/population.h5", mode='r') as store:
    population = store['population']

def calculate_dh_fraction():
    """根据SPH需求计算集中供热比例"""
    # 集中供暖省份
    central_heating_provinces = [
        'Heilongjiang', 'Jilin', 'Liaoning',  # 东北三省
        'Beijing', 'Hebei', 'Shanxi', 'Shandong',  # 华北地区
        'InnerMongolia',  # 内蒙古
        'Gansu', 'Qinghai', 'Ningxia', 'Xinjiang'  # 西北地区
    ]
    
    # 目标集中供暖需求
    target_central_demand = 1402  # TWh
    
    # 计算集中供暖省份的总SPH需求
    total_central_sph = sph_demand[central_heating_provinces].sum()
    
    dh_fraction = pd.Series(index=population.index, dtype=float)
    
    for province in population.index:
        if province in central_heating_provinces:
            # 根据该省份占总SPH需求的比例分配集中供暖
            province_sph = sph_demand.loc[province]
            dh_fraction[province] = (province_sph / total_central_sph) * (target_central_demand / province_sph)
        else:
            # 非集中供暖省份设为0
            dh_fraction[province] = 0.0
    
    return dh_fraction

# 计算集中供暖比例
dh_fraction = calculate_dh_fraction()

# 保存为h5文件
with pd.HDFStore('data/heating/DH_percent2020.h5', mode='w') as store:
    store['central_fraction'] = dh_fraction

# 打印结果和验证
print("已生成 DH_percent2020.h5 文件")
print("\n各省份集中供暖比例:")
for province in dh_fraction.index:
    if dh_fraction[province] > 0:
        print(f"{province}: {dh_fraction[province]:.2f}")

# 验证
print("\n验证:")
print(f"目标集中供暖需求: 1402 TWh")
central_provinces_sph = sph_demand[dh_fraction > 0].sum()
print(f"集中供暖省份总SPH需求: {central_provinces_sph:.2f} TWh")
estimated_central = sum(dh_fraction * sph_demand)
print(f"估算的集中供暖需求: {estimated_central:.2f} TWh")