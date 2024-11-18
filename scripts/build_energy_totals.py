import pandas as pd
import numpy as np
from scripts._helpers import configure_logging

def build_energy_totals():
    """
    Create energy demand totals for heating sector based on planning horizons.
    Historical data from 2020:
    - District heating (central): 1402 TWh
    - Space heating (decentral): 777 TWh
    - Domestic hot water: 516 TWh
    Total: 2695 TWh
    """
    planning_horizons = int(snakemake.wildcards.planning_horizons)
    heating_demand = snakemake.wildcards.heating_demand
    
    # 读取省级人口数据
    pop = pd.read_csv(snakemake.input.population, index_col=0)
    population_values = pop['2020']  # 基准年份人口数据
    
    # 计算每天的生活热水需求 (MWh/day)
    if heating_demand == 'positive':
        # 增长情景：到2060年达到1000 kWh/person/year
        if planning_horizons == 2020:
            annual_hw_per_person = 516 / population_values.sum()  # TWh/person/year in 2020
        elif planning_horizons == 2060:
            annual_hw_per_person = 1.0  # MWh/person/year in 2060
        else:
            # 线性插值
            annual_hw_per_person = np.interp(
                planning_horizons,
                [2020, 2060],
                [516 / population_values.sum(), 1.0]
            )
    else:
        # 保持不变情景：维持2020年水平
        annual_hw_per_person = 516 / population_values.sum()  # TWh/person/year
    
    # 转换为每日需求 (MWh/day)
    daily_hw_per_person = annual_hw_per_person * 1000 / 365  # 转换为kWh/day
    hot_water_per_day = population_values * daily_hw_per_person / 1000  # 转换为MWh/day
    
    # 读取省级供暖需求数据
    with pd.HDFStore(snakemake.input.heat_demand_profile, mode='r') as store:
        heat_demand = store['heat_demand_profiles']
        print("供暖需求数据形状:", heat_demand.shape)
        print("可用列:", heat_demand.columns.tolist())
    
    # 保存到HDF5文件
    with pd.HDFStore(snakemake.output.energy_totals, mode='w') as store:
        store['space_heating_per_hdd'] = heat_demand
        store['hot_water_per_day'] = hot_water_per_day
        
        # 保存元数据
        store['metadata'] = pd.Series({
            'year': planning_horizons,
            'scenario': heating_demand,
            'unit_space_heating': 'normalized heat demand per province',
            'unit_hot_water': 'MWh/day',
            'annual_hot_water_per_person': annual_hw_per_person,
            'total_population': population_values.sum()
        })

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_energy_totals', 
                                  planning_horizons='2020',
                                  heating_demand='positive')
    configure_logging(snakemake)
    build_energy_totals() 