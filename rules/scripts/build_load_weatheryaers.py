import pandas as pd
import numpy as np
import os

def interpolate_load_data(year_before=2040, year_after=2055, target_years=[2045, 2050]):
    """
    通过线性插值生成目标年份的负荷数据
    Parameters:
    -----------
    year_before : int
        起始年份
    year_after : int
        结束年份
    target_years : list
        需要插值的目标年份列表
    """
    # 读取2040和2055年的数据
    file_before = f'data/load/load_{year_before}_weatheryears_2020_2060_TWh.h5'
    file_after = f'data/load/load_{year_after}_weatheryears_2020_2060_TWh.h5'
    
    with pd.HDFStore(file_before, mode='r') as store:
        load_before = store['load']
        print(f"Loaded {file_before}")
    
    with pd.HDFStore(file_after, mode='r') as store:
        load_after = store['load']
        print(f"Loaded {file_after}")
    
    # 对每个目标年份进行插值
    for target_year in target_years:
        # 计算权重
        weight_after = (target_year - year_before) / (year_after - year_before)
        weight_before = 1 - weight_after
        
        # 线性插值
        load_target = weight_before * load_before + weight_after * load_after
        
        # 保存插值结果
        os.makedirs('data/load', exist_ok=True)
        output_file = f'data/load/load_{target_year}_weatheryears_2020_2060_TWh.h5'
        
        with pd.HDFStore(output_file, mode='w') as store:
            store['load'] = load_target
            
        print(f"\nCreated {output_file}")
        print(f"Year {target_year} Statistics:")
        print(f"Weight before ({year_before}): {weight_before:.2f}")
        print(f"Weight after ({year_after}): {weight_after:.2f}")
        print(f"Mean load: {load_target.mean().mean():.2f} TWh")

    # 打印总体统计信息
    print("\nOverall Load Statistics (TWh):")
    print(f"{year_before} mean: {load_before.mean().mean():.2f}")
    for target_year in target_years:
        with pd.HDFStore(f'data/load/load_{target_year}_weatheryears_2020_2060_TWh.h5', mode='r') as store:
            load = store['load']
            print(f"{target_year} mean: {load.mean().mean():.2f}")
    print(f"{year_after} mean: {load_after.mean().mean():.2f}")

if __name__ == "__main__":
    interpolate_load_data(year_before=2040, 
                         year_after=2055, 
                         target_years=[2045, 2050])