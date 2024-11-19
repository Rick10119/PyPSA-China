# SPDX-FileCopyrightText: : 2024 The PyPSA-China Authors
#
# SPDX-License-Identifier: MIT

"""
此文件主要提供三个核心功能：
1. calculate_annuity: 计算年金系数
2. load_costs: 加载和处理成本数据（被多个其他脚本调用）
3. update_transmission_costs: 更新输电成本
"""

import logging
import pandas as pd

idx = pd.IndexSlice

logger = logging.getLogger(__name__)

def calculate_annuity(n, r):
    """Calculate the annuity factor for an asset with lifetime n years and
    discount rate of r, e.g. annuity(20, 0.05) * 20 = 1.6"""
    """计算年金系数，用于将投资成本转换为年化成本
    n: 资产寿命（年）
    r: 贴现率
    """

    if isinstance(r, pd.Series):
        return pd.Series(1/n, index=r.index).where(r == 0, r/(1. - 1./(1.+r)**n))
    elif r > 0:
        return r / (1. - 1./(1.+r)**n)
    else:
        return 1 / n

def load_costs(tech_costs, config, elec_config, cost_year, Nyears):
    """加载和处理技术成本数据
    1. 读取成本数据文件
    2. 单位转换（kW到MW，USD到EUR）
    3. 计算年化成本
    4. 设置默认参数
    """

    # set all asset costs and other parameters
    costs = pd.read_csv(tech_costs, index_col=list(range(3))).sort_index()

    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"),"value"] *= 1e3
    costs.loc[costs.unit.str.contains("USD"),"value"] *= config['USD2013_to_EUR2013']

    cost_year = float(cost_year)
    costs = (costs.loc[idx[:,cost_year,:], "value"]
             .unstack(level=2).groupby("technology").sum(min_count=1))


    costs = costs.fillna({"CO2 intensity" : 0,
                          "FOM" : 0,
                          "VOM" : 0,
                          "discount rate" : config['discountrate'],
                          "efficiency" : 1,
                          "fuel" : 0,
                          "investment" : 0,
                          "lifetime" : 25})

    costs["capital_cost"] = ((calculate_annuity(costs["lifetime"], costs["discount rate"]) +
                             costs["FOM"]/100.) *
                             costs["investment"] * Nyears)

    costs.at['OCGT', 'fuel'] = costs.at['gas', 'fuel']
    costs.at['CCGT', 'fuel'] = costs.at['gas', 'fuel']

    costs['marginal_cost'] = costs['VOM'] + costs['fuel'] / costs['efficiency']

    costs = costs.rename(columns={"CO2 intensity": "co2_emissions"})

    costs.at['OCGT', 'co2_emissions'] = costs.at['gas', 'co2_emissions']
    costs.at['CCGT', 'co2_emissions'] = costs.at['gas', 'co2_emissions']

    costs.at['solar', 'capital_cost'] = 0.5*(costs.at['solar-rooftop', 'capital_cost'] +
                                             costs.at['solar-utility', 'capital_cost'])

    def costs_for_storage(store, link1, link2=None, max_hours=1.):
        capital_cost = link1['capital_cost'] + max_hours * store['capital_cost']
        if link2 is not None:
            capital_cost += link2['capital_cost']
        return pd.Series(dict(capital_cost=capital_cost,
                              marginal_cost=0.,
                              co2_emissions=0.))

    max_hours = elec_config['max_hours']
    costs.loc["battery"] = \
        costs_for_storage(costs.loc["battery storage"], costs.loc["battery inverter"],
                          max_hours=max_hours['battery'])
    costs.loc["H2"] = \
        costs_for_storage(costs.loc["hydrogen storage tank type 1"], costs.loc["fuel cell"],
                          costs.loc["electrolysis"], max_hours=max_hours['H2'])

    for attr in ('marginal_cost', 'capital_cost'):
        overwrites = config.get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites

    return costs

def update_transmission_costs(n, costs, length_factor=1.0):
    """更新输电线路成本
    1. 计算交流线路成本
    2. 计算直流线路成本（考虑海底电缆）
    """
    # TODO: line length factor of lines is applied to lines and links.
    # Separate the function to distinguish.

    n.lines['capital_cost'] = (n.lines['length'] * length_factor *
                               costs.at['HVAC overhead', 'capital_cost'])

    if n.links.empty: return

    dc_b = n.links.carrier == 'DC'

    # If there are no dc links, then the 'underwater_fraction' column
    # may be missing. Therefore we have to return here.
    if n.links.loc[dc_b].empty: return

    costs = (n.links.loc[dc_b, 'length'] * length_factor *
            ((1. - n.links.loc[dc_b, 'underwater_fraction']) *
            costs.at['HVDC overhead', 'capital_cost'] +
            n.links.loc[dc_b, 'underwater_fraction'] *
            costs.at['HVDC submarine', 'capital_cost']) +
            costs.at['HVDC inverter pair', 'capital_cost'])
    n.links.loc[dc_b, 'capital_cost'] = costs

def add_build_year_to_new_assets(n, year):
    """
    为网络中的新增资产添加建设年份属性
    
    参数:
        n: pypsa.Network - 需要处理的网络对象
        year: int - 建设年份
        
    功能:
        1. 遍历网络中的Link、Generator、Store组件
        2. 检查是否已有build_year属性
        3. 如果没有，则添加并设置为指定年份
    """
    # 遍历需要处理的组件类型
    for c in n.iterate_components(["Link", "Generator", "Store"]):
        # 如果组件没有build_year属性
        if "build_year" not in c.df:
            # 添加build_year列并设置为指定年份
            c.df["build_year"] = year