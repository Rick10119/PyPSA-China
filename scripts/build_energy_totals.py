import pandas as pd
import numpy as np
from _helpers import configure_logging

def build_energy_totals():
    """
    Create energy demand totals for heating sector based on planning horizons.
    """
    planning_horizons = int(snakemake.wildcards.planning_horizons)
    
    # 读取省级人口数据
    pop = pd.read_csv("data/population/population_from_National_Data_2020.csv", index_col=0)
    
    # 读取省级供暖需求数据
    sph = pd.read_hdf(snakemake.input.daily_heat_demand)
    print("Available columns:", sph.columns.tolist())
    
    try:
        # 尝试读取数据，如果列名不同则使用替代方案
        space_heating_per_hdd = sph['space_heating_demand']
    except KeyError:
        # 可能的替代列名
        alternative_names = ['heating_demand', 'heat_demand', 'space_heating']
        for name in alternative_names:
            if name in sph.columns:
                space_heating_per_hdd = sph[name]
                print(f"Using alternative column name: {name}")
                break
        else:
            raise KeyError(f"Could not find space heating demand column. Available columns: {sph.columns.tolist()}")
    
    # 计算每天的生活热水需求 (MWh/day)
    # 假设每人每天热水需求为X kWh
    daily_hw_per_person = 1.5  # kWh/person/day
    hot_water_per_day = pop['population'] * daily_hw_per_person / 1000  # 转换为MWh
    
    # 保存到HDF5文件
    with pd.HDFStore(snakemake.output.energy_totals, mode='w') as store:
        store['space_heating_per_hdd'] = space_heating_per_hdd
        store['hot_water_per_day'] = hot_water_per_day
        
        # 保存元数据
        store['metadata'] = pd.Series({
            'year': planning_horizons,
            'unit_space_heating': 'MWh/HDD',
            'unit_hot_water': 'MWh/day'
        })

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_energy_totals', planning_horizons='2020')
    configure_logging(snakemake)
    build_energy_totals() 