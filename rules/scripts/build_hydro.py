import pandas as pd
import numpy as np
from pathlib import Path
import yaml

def create_hydro_profiles(config_path="config.yaml"):
    """
    生成2020年水电出力特性时间序列
    并保存为 hydro_p_max_pu.h5 文件
    """
    # 读取配置文件
    with open(config_path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    # 读取必要的数据
    print("Reading input data...")
    
    # 读取大坝数据
    df = pd.read_csv('data/hydro/dams_large.csv', index_col=0)
    
    # 读取入流数据
    print("Reading inflow data...")
    # 使用2016年的数据，不设置时区
    date_range = pd.date_range('2016-01-01', '2016-12-31 23:00', freq='h')
    inflow = pd.read_pickle('data/hydro/daily_hydro_inflow_per_dam_1979_2016_m3.pickle')
    inflow = inflow.reindex(date_range, method='ffill')  # 重采样到小时数据
    inflow.columns = df.index
    
    # 读取水耗系数
    water_consumption_factor = df.loc[:, 'Water_consumption_factor_avg'] * 1e3  # m^3/KWh -> m^3/MWh
    
    # 确定入流电站
    bus0s = [0, 21, 11, 19, 22, 29, 8, 40, 25, 1, 7, 4, 10, 15, 12, 20, 26, 6, 3, 39]
    bus1s = [5, 11, 19, 22, 32, 8, 40, 25, 35, 2, 4, 10, 9, 12, 20, 23, 6, 17, 14, 16]
    inflow_stations = [dam for dam in range(len(df.index)) if not dam in bus1s]
    
    print(f"Processing {len(inflow_stations)} inflow stations...")
    
    # 计算所有入流电站的p_pu
    all_p_pu = []
    for inflow_station in inflow_stations:
        p_nom = (inflow/water_consumption_factor).iloc[:,inflow_station].max()
        p_pu = (inflow/water_consumption_factor).iloc[:,inflow_station] / p_nom
        all_p_pu.append(p_pu)
    
    # 计算平均出力特性
    print("Calculating average profiles...")
    avg_p_pu = pd.concat(all_p_pu, axis=1).mean(axis=1)
    
    # 创建2020年的出力特性，添加时区信息
    date_range_2020 = pd.date_range('2020-01-01', '2020-12-31 23:00', freq='h', tz='Asia/Shanghai')
    
    # 为所有省份使用相同的出力特性
    year_profile = pd.DataFrame(
        {province: avg_p_pu.values for province in df.index},
        index=date_range_2020
    )
    
    # 确保输出目录存在
    Path('data/p_nom').mkdir(parents=True, exist_ok=True)
    
    # 保存为h5文件
    print("Saving to hydro_p_max_pu.h5...")
    year_profile.to_hdf(
        'data/p_nom/hydro_p_max_pu.h5',
        key='hydro_p_max_pu',
        mode='w'
    )
    
    print("Done! File saved to data/p_nom/hydro_p_max_pu.h5")
    print(f"Shape of the final profile: {year_profile.shape}")
    print(f"Time range: {year_profile.index[0]} to {year_profile.index[-1]}")

if __name__ == "__main__":
    create_hydro_profiles()