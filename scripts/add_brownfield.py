# SPDX-FileCopyrightText: : 2024 The PyPSA-China Authors
#
# SPDX-License-Identifier: MIT

"""
处理多时期规划中的"棕地"情况，考虑已有设施的演变。

输入文件:
    network: 当前时间点的基础网络 (.nc)
    network_p: 上一时期的已求解网络 (.nc)
    costs: 成本参数 (costs_{planning_horizons}.csv)
    profile_{tech}: 可再生能源出力曲线 (.nc)
    overrides: 网络属性覆盖

输出文件:
    network_name: 更新后的网络文件 (.nc)

工作流程:
1. 读取上一时期的网络结果
2. 保留未到期的发电和储能设施
3. 更新输电线路容量
4. 更新可再生能源开发潜力
5. 输出更新后的网络
"""

import logging
logger = logging.getLogger(__name__)

import pandas as pd
idx = pd.IndexSlice

import pypsa
import numpy as np
import xarray as xr

from add_existing_baseyear import add_build_year_to_new_assets
from _helpers import override_component_attrs
from functions import pro_names, offwind_nodes


def basename(x):
    return x.split("-2")[0]

def add_brownfield(n, n_p, year):
    """
    将上一时期的设施添加到当前网络中。
    
    参数:
        n: pypsa.Network, 当前时间点的网络
        n_p: pypsa.Network, 上一时期的已求解网络
        year: int, 规划年份
    
    主要步骤:
    1. 保留上一时期的输电容量作为下限
    2. 移除寿命到期的设施
    3. 固定仍在运行的设施容量
    4. 更新可再生能源潜力
    """

    print("adding brownfield")

    # 1. 电网输电容量处理：设置上一期优化结果为最小值
    n.lines.s_nom_min = n_p.lines.s_nom_opt
    n.links.loc[(n.links.length>0) & (n.links.lifetime==np.inf),"p_nom"] = n_p.links.loc[(n_p.links.carrier=='AC') & (n_p.links.build_year==0),"p_nom_opt"]
    n.links.loc[(n.links.length>0) & (n.links.lifetime==np.inf),"p_nom_min"] = n_p.links.loc[(n_p.links.carrier=='AC') & (n_p.links.build_year==0),"p_nom_opt"]

    # 2025年特殊处理
    if year == 2025:
        add_build_year_to_new_assets(n_p, 2023)

    # 2. 处理发电、储能和连接设施
    for c in n_p.iterate_components(["Link", "Generator", "Store"]):
        attr = "e" if c.name == "Store" else "p"

        # first, remove generators, links and stores that track
        # CO2 or global EU values since these are already in n
        n_p.mremove(
            c.name,
            c.df.index[c.df.lifetime==np.inf]
        )
        # 移除寿命到期的设施
        n_p.mremove(
            c.name,
            c.df.index[c.df.build_year + c.df.lifetime < year]
        )
        # remove assets if their optimized nominal capacity is lower than a threshold
        # since CHP heat Link is proportional to CHP electric Link, make sure threshold is compatible
        chp_heat = c.df.index[(
            c.df[attr + "_nom_extendable"]
            & c.df.index.str.contains("CHP")
        )]

        threshold = snakemake.config['existing_capacities']['threshold_capacity']

        if not chp_heat.empty:
            threshold_chp_heat = threshold*c.df.loc[chp_heat].efficiency2/c.df.loc[chp_heat].efficiency
            n_p.mremove(
                c.name,
                chp_heat[c.df.loc[chp_heat, attr + "_nom_opt"] < threshold_chp_heat]
            )

        # 移除低于阈值的设施
        n_p.mremove(
            c.name,
            c.df.index[c.df[attr + "_nom_extendable"] & ~c.df.index.isin(chp_heat) & (c.df[attr + "_nom_opt"] < threshold)]
        )

        # 3. 固定保留设施的容量
        c.df[attr + "_nom"] = c.df[attr + "_nom_opt"]
        c.df[attr + "_nom_extendable"] = False
        c.df[attr + "_nom_max"] = np.inf

        # 导入组件到新网络
        n.import_components_from_dataframe(c.df, c.name)

        # 导入时间序列数据
        selection = (
            n.component_attrs[c.name].type.str.contains("series")
            & n.component_attrs[c.name].status.str.contains("Input")
        )

        for tattr in n.component_attrs[c.name].index[selection]:
            n.import_series_from_dataframe(c.pnl[tattr].set_index(n.snapshots), c.name, tattr)

    # 4. 更新可再生能源潜力
    for tech in ['onwind', 'offwind', 'solar']:
        ds_tech = xr.open_dataset(snakemake.input['profile_' + tech])
        p_nom_max_initial = ds_tech['p_nom_max'].to_pandas()

        if tech == 'offwind':
            for node in offwind_nodes:
                # 更新海上风电潜力
                n.generators.loc[(n.generators.bus == node) & (n.generators.carrier == tech) & (n.generators.build_year == year),"p_nom_max"] = \
                p_nom_max_initial[node] - n_p.generators[(n_p.generators.bus == node) & (n_p.generators.carrier == tech)].p_nom_opt.sum()
        else:
            for node in pro_names:
                # 更新陆上风电和光伏潜力
                n.generators.loc[(n.generators.bus == node) & (n.generators.carrier == tech) & (n.generators.build_year == year),"p_nom_max"] = \
                p_nom_max_initial[node] - n_p.generators[(n_p.generators.bus == node) & (n_p.generators.carrier == tech)].p_nom_opt.sum()

    # 确保潜力不为负
    n.generators.loc[(n.generators.p_nom_max < 0), "p_nom_max"] = 0

    # 5. 煤电厂碳捕集改造设置
    n.generators.loc[n.generators.carrier == 'coal power plant', 'p_nom_extendable'] = True
    n.generators.loc[n.generators.index.str.contains("retrofit") & ~n.generators.index.str.contains(
        str(year)), "p_nom_extendable"] = False


if __name__ == '__main__':
    # 测试代码部分
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('add_brownfield',
                                   opts='ll',
                                   topology='current+FCG',
                                   pathway='linear-275',
                                   co2_reduction='1.0',
                                   planning_horizons=2025)

    print(snakemake.input.network_p)
    logging.basicConfig(level=snakemake.config['logging_level'])

    year = int(snakemake.wildcards.planning_horizons)

    # 读取网络
    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

    add_build_year_to_new_assets(n, year)

    n_p = pypsa.Network(snakemake.input.network_p, override_component_attrs=overrides)

    # 执行brownfield添加
    add_brownfield(n, n_p, year)

    # 保存结果
    n.export_to_netcdf(snakemake.output.network_name)
