import pandas as pd
import numpy as np
import os

def interpolate_load_data(year_before=2040, year_after=2055, target_years=[2045, 2050]):
    """
    通过线性插值生成目标年份的负荷数据
    """
    # 读取参考文件，检查其属性
    file_before = f'data/load/load_{year_before}_weatheryears_2020_2060_TWh.h5'
    
    with pd.HDFStore(file_before, mode='r') as store:
        load_before = store['load']
        print("Reference file info:")
        print(store.info())
        print("\nDataFrame info:")
        print(load_before.info())
        
    file_after = f'data/load/load_{year_after}_weatheryears_2020_2060_TWh.h5'
    with pd.HDFStore(file_after, mode='r') as store:
        load_after = store['load']
    
    # 对每个目标年份进行插值
    for target_year in target_years:
        weight_after = (target_year - year_before) / (year_after - year_before)
        weight_before = 1 - weight_after
        
        load_target = weight_before * load_before + weight_after * load_after
        
        # 确保数据类型和精度与原文件一致
        load_target = load_target.astype(load_before.dtypes.iloc[0])
        
        output_file = f'data/load/load_{target_year}_weatheryears_2020_2060_TWh.h5'
        
        # 使用fixed格式但不使用压缩
        with pd.HDFStore(output_file, mode='w') as store:
            store.put('load', load_target, format='fixed')
            
        print(f"\nCreated {output_file}")
        print(f"File size: {os.path.getsize(output_file)/1024:.2f} KB")

    # 打印统计信息
    print("\nFile sizes (KB):")
    print(f"{year_before}: {os.path.getsize(file_before)/1024:.2f}")
    for target_year in target_years:
        filename = f'data/load/load_{target_year}_weatheryears_2020_2060_TWh.h5'
        print(f"{target_year}: {os.path.getsize(filename)/1024:.2f}")
    print(f"{year_after}: {os.path.getsize(file_after)/1024:.2f}")

if __name__ == "__main__":
    interpolate_load_data()