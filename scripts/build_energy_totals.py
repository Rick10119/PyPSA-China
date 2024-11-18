import pandas as pd
import numpy as np
from _helpers import configure_logging

def build_energy_totals():
    """
    Create energy demand totals for heating sector based on planning horizons.
    """
    planning_horizons = int(snakemake.wildcards.planning_horizons)
    
    # 读取省级人口数据
    pop = pd.read_csv(snakemake.input.population, index_col=0)
    
    # 读取省级供暖需求数据
    with pd.HDFStore(snakemake.input.heat_demand_profile, mode='r') as store:
        # 从heat_demand_profiles读取数据
        heat_demand = store['heat_demand_profiles']
        print("数据形状:", heat_demand.shape)
        print("可用列:", heat_demand.columns.tolist())
    
    # 计算每天的生活热水需求 (MWh/day)
    daily_hw_per_person = 1.5  # kWh/person/day
    hot_water_per_day = pop['population'] * daily_hw_per_person / 1000  # 转换为MWh
    
    # 保存到HDF5文件
    with pd.HDFStore(snakemake.output.energy_totals, mode='w') as store:
        store['space_heating_per_hdd'] = heat_demand
        store['hot_water_per_day'] = hot_water_per_day
        
        # 保存元数据
        store['metadata'] = pd.Series({
            'year': planning_horizons,
            'unit_space_heating': 'normalized heat demand per province',
            'unit_hot_water': 'MWh/day'
        })

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_energy_totals', planning_horizons='2020')
    configure_logging(snakemake)
    build_energy_totals() 