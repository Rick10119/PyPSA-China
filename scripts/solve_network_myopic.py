# SPDX-FileCopyrightText: : 2022 The PyPSA-China Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8

import logging
import re
import copy
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
import os
import time
from _helpers import (
    configure_logging,
    override_component_attrs,
)

# 添加并行计算相关的导入
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import functools

# 允许传递给PyPSA optimize的参数
ALLOWED_OPTIMIZE_KWARGS = [
    "solver_name", "solver_options", "extra_functionality",
    "track_iterations", "min_iterations", "max_iterations"
]

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)

# 设置gurobipy和linopy的日志级别，避免输出info信息
import logging
gurobipy_logger = logging.getLogger('gurobipy')
gurobipy_logger.setLevel(logging.WARNING)

linopy_logger = logging.getLogger('linopy')
linopy_logger.setLevel(logging.WARNING)

# 设置linopy.model的日志级别
linopy_model_logger = logging.getLogger('linopy.model')
linopy_model_logger.setLevel(logging.WARNING)

from pypsa.descriptors import get_switchable_as_dense as get_as_dense

def prepare_network(
        n,
        solve_opts=None,
        using_single_node=False,
        single_node_province="Shandong"
):
    # 检查是否启用单节点模式
    if using_single_node:
        # Filter to keep only specified province components
        province_buses = n.buses[n.buses.index.str.contains(single_node_province)].index
        
        # Remove non-province buses and their components
        non_province_buses = n.buses[~n.buses.index.isin(province_buses)].index
        
        # Remove generators not in specified province
        non_province_generators = n.generators[~n.generators.bus.isin(province_buses)].index
        n.mremove("Generator", non_province_generators)
        
        # Remove loads not in specified province
        non_province_loads = n.loads[~n.loads.bus.isin(province_buses)].index
        n.mremove("Load", non_province_loads)
        
        # Remove storage units not in specified province
        non_province_storage = n.storage_units[~n.storage_units.bus.isin(province_buses)].index
        n.mremove("StorageUnit", non_province_storage)
        
        # Remove stores not in specified province
        non_province_stores = n.stores[~n.stores.bus.isin(province_buses)].index
        n.mremove("Store", non_province_stores)
        
        # Remove links not connected to specified province
        non_province_links = n.links[~(n.links.bus0.isin(province_buses) | n.links.bus1.isin(province_buses))].index
        n.mremove("Link", non_province_links)
        
        # Remove lines not connected to specified province
        non_province_lines = n.lines[~(n.lines.bus0.isin(province_buses) | n.lines.bus1.isin(province_buses))].index
        n.mremove("Line", non_province_lines)
        
        # Finally remove non-province buses
        n.mremove("Bus", non_province_buses)
    
    # Fix any remaining links that might have undefined buses
    for link in n.links.index:
        if not (n.links.at[link, "bus0"] in n.buses.index and n.links.at[link, "bus1"] in n.buses.index):
            n.mremove("Link", [link])

    if "clip_p_max_pu" in solve_opts:
        for df in (
            n.generators_t.p_max_pu,
            n.generators_t.p_min_pu,
            n.storage_units_t.inflow,
        ):
            df.where(df > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    load_shedding = solve_opts.get("load_shedding")
    if load_shedding:
        # intersect between macroeconomic and surveybased willingness to pay
        # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
        # TODO: retrieve color and nice name from config
        n.add("Carrier", "load", color="#dd2e23", nice_name="Load shedding")
        buses_i = n.buses.query("carrier == 'AC'").index
        if not np.isscalar(load_shedding):
            # TODO: do not scale via sign attribute (use Eur/MWh instead of Eur/kWh)
            load_shedding = 1e2  # Eur/kWh

        n.madd(
            "Generator",
            buses_i,
            " load",
            bus=buses_i,
            carrier="load",
            sign=1e-3,  # Adjust sign to measure p and p_nom in kW instead of MW
            marginal_cost=load_shedding,  # Eur/kWh
            p_nom=1e9,  # kW
        )

    if solve_opts.get("noisy_costs"):
        for t in n.iterate_components():
            # if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if "marginal_cost" in t.df:
                t.df["marginal_cost"] += 1e-2 + 2e-3 * (
                    np.random.random(len(t.df)) - 0.5
                )

        for t in n.iterate_components(["Line", "Link"]):
            t.df["capital_cost"] += (
                1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)
            ) * t.df["length"]

    if solve_opts.get('nhours'):
        nhours = solve_opts['nhours']
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760. / nhours

    return n

def add_chp_constraints(n):
    electric = (
        n.links.index.str.contains("CHP")
        & n.links.index.str.contains("generator")
    )
    heat = (
        n.links.index.str.contains("CHP")
        & n.links.index.str.contains("boiler")
    )

    electric_ext = n.links[electric].query("p_nom_extendable").index
    heat_ext = n.links[heat].query("p_nom_extendable").index

    electric_fix = n.links[electric].query("~p_nom_extendable").index
    heat_fix = n.links[heat].query("~p_nom_extendable").index

    p = n.model["Link-p"]  # dimension: [time, link]

    # output ratio between heat and electricity and top_iso_fuel_line for extendable
    if not electric_ext.empty:
        p_nom = n.model["Link-p_nom"]
        
        # Scale factors to improve numerical stability
        scale_factor = 1e-3  # Scale down by 1000
        
        # Get efficiency ratios with scaling
        elec_eff = (n.links.p_nom_ratio * n.links.efficiency)[electric_ext].values * scale_factor
        heat_eff = n.links.efficiency[heat_ext].values * scale_factor
        
        # Add constraint with scaled values
        lhs = (
            p_nom.loc[electric_ext] * elec_eff
            - p_nom.loc[heat_ext] * heat_eff
        )
        n.model.add_constraints(lhs == 0, name="chplink-fix_p_nom_ratio")

        # Scale the top_iso_fuel_line constraint
        rename = {"Link-ext": "Link"}
        lhs = (
            p.loc[:, electric_ext] * scale_factor
            + p.loc[:, heat_ext] * scale_factor
            - p_nom.rename(rename).loc[electric_ext] * scale_factor
        )
        n.model.add_constraints(lhs <= 0, name="chplink-top_iso_fuel_line_ext")

    # top_iso_fuel_line for fixed
    if not electric_fix.empty:
        # Scale the fixed capacity constraint
        scale_factor = 1e-3
        lhs = p.loc[:, electric_fix] * scale_factor + p.loc[:, heat_fix] * scale_factor
        rhs = n.links.p_nom[electric_fix] * scale_factor
        n.model.add_constraints(lhs <= rhs, name="chplink-top_iso_fuel_line_fix")

    # back-pressure
    if not n.links[electric].index.empty:
        # Scale the back-pressure constraint
        scale_factor = 1e-3
        lhs = (
            p.loc[:, heat] * (n.links.efficiency[heat] * n.links.c_b[electric].values) * scale_factor
            - p.loc[:, electric] * n.links.efficiency[electric] * scale_factor
        )
        n.model.add_constraints(lhs <= 0, name="chplink-backpressure")

def add_transimission_constraints(n):
    """
    Add constraints for transmission lines that allow for asymmetric capacities
    while maintaining reasonable limits on the difference between directions.
    """
    if not n.links.p_nom_extendable.any():
        return

    positive_bool = n.links.index.str.contains("positive")
    negative_bool = n.links.index.str.contains("reversed")

    positive_ext = n.links[positive_bool].query("p_nom_extendable").index
    negative_ext = n.links[negative_bool].query("p_nom_extendable").index

    # Allow asymmetric capacities but limit the difference
    for pos in positive_ext:
        neg = pos.replace("positive", "reversed")
        if neg not in negative_ext:
            continue
            
        # Get the current capacities
        pos_cap = n.links.at[pos, "p_nom"]
        neg_cap = n.links.at[neg, "p_nom"]
        
        # Allow up to 20% difference between directions
        max_diff = 0.2 * max(pos_cap, neg_cap)
        
        # Add constraint: |pos_cap - neg_cap| <= max_diff
        lhs = n.model["Link-p_nom"].loc[pos] - n.model["Link-p_nom"].loc[neg]
        n.model.add_constraints(lhs <= max_diff, name=f"Link-transmission-{pos}-max")
        n.model.add_constraints(lhs >= -max_diff, name=f"Link-transmission-{pos}-min")

def add_retrofit_constraints(n):
    p_nom_max = pd.read_csv("data/p_nom/p_nom_max_cc.csv",index_col=0)
    planning_horizon = snakemake.wildcards.planning_horizons
    for year in range(int(planning_horizon) - 40, 2021, 5):
        coal = n.generators[(n.generators.carrier=="coal power plant") & (n.generators.build_year==year)].query("p_nom_extendable").index
        Bus = n.generators[(n.generators.carrier == "coal power plant") & (n.generators.build_year == year)].query(
            "p_nom_extendable").bus.values
        coal_retrofit = n.generators[n.generators.index.str.contains("retrofit")& (n.generators.build_year==year) & n.generators.bus.isin(Bus)].query("p_nom_extendable").index
        coal_retrofitted = n.generators[n.generators.index.str.contains("retrofit") & (n.generators.build_year==year) & n.generators.bus.isin(Bus)].query("~p_nom_extendable").groupby("bus").sum().p_nom_opt

        # Create a Series with proper index for the available capacity
        available_capacity = pd.Series(
            (p_nom_max[str(year)].loc[Bus] - coal_retrofitted.reindex(p_nom_max[str(year)].loc[Bus].index,fill_value=0)).values,
            index=Bus
        )

        # Create constraint with proper dimension handling
        for bus in Bus:
            if bus in coal_retrofitted.index:
                retrofit_cap = coal_retrofitted[bus]
            else:
                retrofit_cap = 0
                
            max_cap = p_nom_max[str(year)].loc[bus]
            coal_bus = coal[n.generators.loc[coal].bus == bus]
            coal_retrofit_bus = coal_retrofit[n.generators.loc[coal_retrofit].bus == bus]
            
            if len(coal_bus) > 0 or len(coal_retrofit_bus) > 0:
                lhs = n.model["Generator-p_nom"].loc[coal_bus].sum() + n.model["Generator-p_nom"].loc[coal_retrofit_bus].sum()
                rhs = max_cap - retrofit_cap
                n.model.add_constraints(lhs == rhs, name=f"Generator-coal-retrofit-{year}-{bus}")

def extra_functionality(n, snapshots, fixed_aluminum_usage=None):
    """
    Collects supplementary constraints which will be passed to ``pypsa.linopf.network_lopf``.
    If you want to enforce additional custom constraints, this is a good location to add them.
    The arguments ``opts`` and ``snakemake.config`` are expected to be attached to the network.
    """
    opts = n.opts
    config = n.config
    add_chp_constraints(n)
    add_transimission_constraints(n)
    if snakemake.wildcards.planning_horizons != "2020":
        add_retrofit_constraints(n)

def solve_aluminum_optimization(n, config, solving, opts="", nodal_prices=None, target_province=None, national_smelter_production=None, **kwargs):
    """
    电解铝最优运行问题求解 - 单省份版本
    基于节点电价，运行以满足铝需求为约束的电解铝最优运行问题
    只求解指定省份的电解铝优化问题
    
    Parameters:
    -----------
    target_province : str
        目标省份名称，如"Shandong"、"Henan"等
    """
    set_of_options = solving["solver"]["options"]
    solver_options = solving["solver_options"][set_of_options] if set_of_options else {}
    solver_name = solving["solver"]["name"]
    
    # 为MILP问题使用专门的求解器设置
    if "MILP" in solving["solver_options"]:
        milp_solver_options = solving["solver_options"]["MILP"]
    else:
        milp_solver_options = solver_options
    
    # 读取铝冶炼厂年产量数据并过滤
    # 从snakemake.input中获取al_smelter_p_max文件路径
    if 'snakemake' in globals():
        al_smelter_p_nom_path = snakemake.input.al_smelter_p_max
    else:
        # 如果没有snakemake，使用默认路径
        al_smelter_p_nom_path = "data/p_nom/al_smelter_p_max.csv"
    
    al_smelter_annual_production = pd.read_csv(al_smelter_p_nom_path)
    al_smelter_annual_production = al_smelter_annual_production.set_index('Province')['p_nom']
    
    # 过滤出年产量大于0.01 10kt/year的省份（与prepare_base_network保持一致）
    al_smelter_annual_production = al_smelter_annual_production[al_smelter_annual_production > 0.01]
    
    # Convert annual production (10kt/year) to power capacity (MW)
    # 1 ton of aluminum requires ~13.3 MWh of electricity
    # Convert 10kt/year to MW: (10kt/year * 10000 * 13.3 MWh/ton) / (8760 hours/year) = MW
    al_smelter_p_nom = al_smelter_annual_production * 10000 * 13.3 / 8760  # Convert to MW
    
    # 如果没有指定目标省份，返回None
    if target_province is None:
        return None
    
    # 检查目标省份是否在铝冶炼厂年产量列表中
    if target_province not in al_smelter_p_nom.index:
        return None
    
    # 找到指定省份的电解铝相关组件
    aluminum_buses = n.buses[n.buses.carrier == "aluminum"].index
    aluminum_smelters = n.links[n.links.carrier == "aluminum"].index
    aluminum_stores = n.stores[n.stores.carrier == "aluminum"].index
    aluminum_loads = n.loads[n.loads.bus.isin(aluminum_buses)].index
    
    # 过滤出指定省份的电解铝组件
    target_aluminum_buses = [bus for bus in aluminum_buses if target_province in bus]
    # 过滤出指定省份的电解槽
    target_aluminum_smelters = [smelter for smelter in aluminum_smelters if target_province in smelter]
    target_aluminum_stores = [store for store in aluminum_stores if target_province in store]
    target_aluminum_loads = [load for load in aluminum_loads if target_province in load]
    
    # 如果没有找到该省份的电解铝组件，返回None
    if not target_aluminum_smelters:
        return None
    
    # 获取电解铝厂运行参数
    from scripts.scenario_utils import get_aluminum_smelter_operational_params
    
    # 重新设置电解铝冶炼设备的参数
    for smelter in target_aluminum_smelters:
        n.links.at[smelter, 'committable'] = config['aluminum_commitment']
        # 获取运行参数
        operational_params = get_aluminum_smelter_operational_params(config, al_smelter_p_nom=al_smelter_p_nom[target_province])
        n.links.at[smelter, 'p_min_pu'] = operational_params['p_min_pu'] if config['aluminum_commitment'] else 0
    
    # 移除所有非电解铝相关的组件
    for component_type in ["Generator", "StorageUnit", "Store", "Link", "Load"]:
        if component_type == "Store":
            # 保留指定省份的电解铝存储
            other_stores = n.stores[~n.stores.index.isin(target_aluminum_stores)].index
            n.mremove(component_type, other_stores)
        elif component_type == "Link":
            # 保留指定省份的电解铝冶炼设备
            other_links = n.links[~n.links.index.isin(target_aluminum_smelters)].index
            n.mremove(component_type, other_links)
        elif component_type == "Load":
            # 保留指定省份的电解铝负荷
            other_loads = n.loads[~n.loads.index.isin(target_aluminum_loads)].index
            n.mremove(component_type, other_loads)
        else:
            # 移除所有其他组件
            n.mremove(component_type, n.df(component_type).index)
    
    # 移除aluminum和AC以外的节点
    non_aluminum_buses = n.buses[(n.buses.carrier != "aluminum") & (n.buses.carrier != "AC")].index
    n.mremove("Bus", non_aluminum_buses)
    
    # 确保虚拟carrier存在
    if "virtual" not in n.carriers.index:
        n.add("Carrier", "virtual")
    
    # 添加虚拟发电机来提供电力（基于节点电价）
    for bus in n.buses.index:
        if bus.endswith(" aluminum"):
            # 这是电解铝节点，不需要虚拟发电机
            continue
            
        # 根据节点类型确定边际成本
        if bus in nodal_prices.columns:
            # 如果该节点有对应的边际电价，使用该节点的电价
            marginal_cost = nodal_prices[bus]
            
        n.add("Generator",
            f"virtual_gen_{bus}",
            bus=bus,
            carrier="virtual",
            p_nom=1e6,  # 大容量
            marginal_cost=marginal_cost)  # 使用节点边际电价作为边际成本
    
    # 使用传入的全国优化结果中该省份的铝冶炼厂产量平均值作为负荷
    target_aluminum_load = target_aluminum_loads[0] if target_aluminum_loads else None
    
    if target_aluminum_load and national_smelter_production and target_province in national_smelter_production:
        # 更新负荷值为平均产量
        n.loads.at[target_aluminum_load, 'p_set'] = national_smelter_production[target_province]
    
    # 定义铝优化专用的extra_functionality函数
    def aluminum_extra_functionality(n_al, snapshots_al):
        # 找到所有铝冶炼厂
        aluminum_smelters = n_al.links[n_al.links.carrier == "aluminum"].index
        
        if len(aluminum_smelters) > 0 and len(snapshots_al) > 0:
            first_snapshot = snapshots_al[0]
            last_snapshot = snapshots_al[-1]
            
            # 为每个铝冶炼厂添加第一个和最后一个时段功率相等约束
            for smelter in aluminum_smelters:
                # 获取该冶炼厂的功率变量 - 注意索引顺序：先时间，再组件
                first_smelter_p = n_al.model.variables["Link-p"].loc[first_snapshot, smelter]
                last_smelter_p = n_al.model.variables["Link-p"].loc[last_snapshot, smelter]
                
                # 添加约束：第一个时段的功率等于最后一个时段的功率
                n_al.model.add_constraints(
                    first_smelter_p == last_smelter_p, 
                    name=f"aluminum-first-last-equal-{smelter}"
                )    
    # 求解电解铝优化问题
    # 只保留PyPSA支持的参数
    optimize_kwargs = {k: v for k, v in kwargs.items() if k in ALLOWED_OPTIMIZE_KWARGS}
    status, condition = n.optimize(
        solver_name=solver_name,
        extra_functionality=aluminum_extra_functionality,
        **milp_solver_options,  # 使用MILP专用参数
        **optimize_kwargs,
    )
    
    # 提取电解铝用能模式
    aluminum_usage = n.links_t.p0[target_aluminum_smelters].copy()
    return aluminum_usage

def solve_aluminum_optimization_parallel_wrapper(args):
    """
    并行化电解铝优化的包装函数
    这个函数需要在独立的进程中运行，因此需要接收所有必要的参数
    
    Parameters:
    -----------
    args : tuple
        包含所有必要参数的元组:
        (original_network_path, original_overrides, config, solving, opts, 
         nodal_prices, target_province, solve_opts, using_single_node, 
         single_node_province, overrides_path, national_smelter_production, **kwargs)
    
    Returns:
    --------
    tuple : (province, aluminum_usage) 或 (province, None)
    """
    # 解包参数
    (original_network_path, original_overrides, config, solving, opts, 
     nodal_prices, target_province, solve_opts, using_single_node, 
     single_node_province, overrides_path, national_smelter_production, kwargs_dict) = args
    
    # 重新创建网络
    if original_overrides:
        n_province = pypsa.Network(original_network_path, override_component_attrs=original_overrides)
    else:
        n_province = pypsa.Network(original_network_path)
    
    # 设置网络文件路径
    n_province._network_path = original_network_path
    if original_overrides:
        n_province._overrides_path = overrides_path
    
    # 重新应用网络准备
    n_province = prepare_network(
        n_province,
        solve_opts,
        using_single_node=using_single_node,
        single_node_province=single_node_province
    )
    
    # 设置配置
    n_province.config = config
    n_province.opts = opts
    
    # 求解单个省份的电解铝优化问题
    province_aluminum_usage = solve_aluminum_optimization(
        n_province, 
        config, 
        solving, 
        opts, 
        nodal_prices=nodal_prices, 
        target_province=target_province,
        national_smelter_production=national_smelter_production,
        **kwargs_dict
    )
    
    return (target_province, province_aluminum_usage)

def solve_aluminum_optimization_parallel(target_provinces, original_network_path, original_overrides, 
                                       config, solving, opts, current_nodal_prices, 
                                       solve_opts, using_single_node, single_node_province, 
                                       overrides_path, national_smelter_production=None, max_workers=None, **kwargs):
    """
    并行求解多个省份的电解铝优化问题
    
    Parameters:
    -----------
    target_provinces : list
        需要优化的省份列表
    original_network_path : str
        原始网络文件路径
    original_overrides : dict or None
        原始覆盖配置
    config : dict
        配置字典
    solving : dict
        求解配置
    opts : list
        选项列表
    current_nodal_prices : pd.DataFrame
        当前节点电价
    solve_opts : dict
        求解选项
    using_single_node : bool
        是否使用单节点模式
    single_node_province : str
        单节点省份名称
    overrides_path : str
        覆盖文件路径
    max_workers : int, optional
        最大并行进程数，默认为CPU核心数
    **kwargs : dict
        其他参数
    
    Returns:
    --------
    dict : {province: aluminum_usage} 或 None
    """
    if not target_provinces:
        return {}
    
    # 设置最大并行进程数
    if max_workers is None:
        max_workers = min(len(target_provinces), mp.cpu_count())
    
    # 准备参数
    args_list = []
    for province in target_provinces:
        args = (original_network_path, original_overrides, config, solving, opts,
                current_nodal_prices, province, solve_opts, using_single_node,
                single_node_province, overrides_path, national_smelter_production, kwargs)
        args_list.append(args)
    
    # 使用进程池并行执行
    all_aluminum_usage = {}
    
    start_time = time.time()
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # 提交所有任务
        future_to_province = {
            executor.submit(solve_aluminum_optimization_parallel_wrapper, args): args[6]  # args[6] 是 target_province
            for args in args_list
        }
        
        # 收集结果
        for future in as_completed(future_to_province):
            result_province, province_aluminum_usage = future.result()
            
            if province_aluminum_usage is not None:
                all_aluminum_usage[result_province] = province_aluminum_usage
    
    total_time = time.time() - start_time
    
    return all_aluminum_usage

def solve_network_iterative(n, config, solving, opts="", max_iterations=10, convergence_tolerance=0.01, **kwargs):
    """
    电解铝迭代优化算法
    1. 使用连续化电解铝模型求解，得到节点电价和目标函数值
    2. 检查收敛性：比较步骤一的目标函数值变化
    3. 基于节点电价，运行电解铝最优运行问题，得到新的电解铝用能模式
    4. 固定电解铝用能，求解剩下的优化问题
    5. 重复1-4，直到收敛
    
    收敛条件：步骤一目标函数变化的相对值小于convergence_tolerance（默认1%）
    """
    import time
    
    set_of_options = solving["solver"]["options"]
    solver_options = solving["solver_options"][set_of_options] if set_of_options else {}
    solver_name = solving["solver"]["name"]
    cf_solving = solving["options"]
    track_iterations = cf_solving.get("track_iterations", False)
    min_iterations = cf_solving.get("min_iterations", 4)
    max_transmission_iterations = cf_solving.get("max_iterations", 6)

    skip_iterations = cf_solving.get("skip_iterations", False)
    if not n.lines.s_nom_extendable.any():
        skip_iterations = True
    

    
    # 读取铝冶炼厂年产量数据，获取需要优化的省份列表
    if 'snakemake' in globals():
        al_smelter_p_nom_path = snakemake.input.al_smelter_p_max
    else:
        al_smelter_p_nom_path = "data/p_nom/al_smelter_p_max.csv"
    
    al_smelter_annual_production = pd.read_csv(al_smelter_p_nom_path)
    al_smelter_annual_production = al_smelter_annual_production.set_index('Province')['p_nom']
    al_smelter_annual_production = al_smelter_annual_production[al_smelter_annual_production > 0.01]
    
    # Convert annual production (10kt/year) to power capacity (MW)
    al_smelter_p_nom = al_smelter_annual_production * 10000 * 13.3 / 8760  # Convert to MW
    
    # 检查哪些省份在网络中实际存在电解铝组件
    available_provinces = []
    for province in al_smelter_p_nom.index:
        # 检查该省份是否有电解铝冶炼设备
        aluminum_smelters = n.links[n.links.carrier == "aluminum"].index
        # 过滤出该省份的电解槽
        province_smelters = [smelter for smelter in aluminum_smelters if province in smelter]
        
        if province_smelters:
            available_provinces.append(province)
    
    target_provinces = available_provinces
    
    if not target_provinces:
        return None
    
    # 记录总开始时间
    total_start_time = time.time()
    
    # 初始化电解铝用能模式和目标函数值
    aluminum_usage = None
    previous_objective = None
    iteration = 0
    final_network = None
    
    # 记录每次迭代的时间
    iteration_times = []
    
    # 记录求解器性能统计
    solver_performance = {
        'iteration_times': [],
        'objective_values': [],
        'convergence_info': []
    }
    
    # 保存原始网络文件路径和overrides，用于重新加载
    original_network_path = None
    original_overrides = None
    if hasattr(n, '_network_path'):
        original_network_path = n._network_path
    if hasattr(n, '_overrides_path'):
        original_overrides = override_component_attrs(n._overrides_path)
    
    # 获取单节点参数
    using_single_node = kwargs.get("using_single_node", False)
    single_node_province = kwargs.get("single_node_province", "Shandong")
    

    
    # 初始化当前网络 - 只在第一次迭代时创建
    n_current = None
    
    while iteration < max_iterations:
        iteration += 1
        iteration_start_time = time.time()  # 记录每次迭代开始时间
        

        
        # 网络初始化策略：使用上一次的结果作为初值
        if iteration == 1:
            # 第一次迭代：重新加载网络
            if original_network_path:
                if original_overrides:
                    n_current = pypsa.Network(original_network_path, override_component_attrs=original_overrides)
                else:
                    n_current = pypsa.Network(original_network_path)
                
                # 设置网络文件路径，用于后续可能的重新加载
                n_current._network_path = original_network_path
                if original_overrides:
                    n_current._overrides_path = kwargs.get("overrides_path", "data/override_component_attrs")
            
            # 重新应用网络准备
            n_current = prepare_network(
                n_current,
                kwargs.get("solve_opts", {}),
                using_single_node=using_single_node,
                single_node_province=single_node_province
            )
            
            # 设置配置
            n_current.config = config
            n_current.opts = opts
        else:
            # 后续迭代：使用上一次的网络结果作为初值
            # n_current已经是上一次迭代的结果，直接使用
            
            # 清除上一次的求解结果，但保留网络结构
            if hasattr(n_current, 'model'):
                del n_current.model
            if hasattr(n_current, 'objective'):
                del n_current.objective
        
        # 重新创建电解铝冶炼设备，确保在步骤1中禁用电解铝启停约束
        aluminum_smelters = n_current.links[n_current.links.carrier == "aluminum"].index
        
        # 保存原始参数
        smelter_params = {}
        for smelter in aluminum_smelters:
            smelter_params[smelter] = {
                'bus0': n_current.links.at[smelter, 'bus0'],
                'bus1': n_current.links.at[smelter, 'bus1'],
                'carrier': n_current.links.at[smelter, 'carrier'],
                'p_nom': n_current.links.at[smelter, 'p_nom'],
                'p_nom_extendable': n_current.links.at[smelter, 'p_nom_extendable'],
                'efficiency': n_current.links.at[smelter, 'efficiency'],
                'marginal_cost': n_current.links.at[smelter, 'marginal_cost'],
                'capital_cost': n_current.links.at[smelter, 'capital_cost'],
                'stand_by_cost': n_current.links.at[smelter, 'stand_by_cost'],
                'shut_down_cost': n_current.links.at[smelter, 'shut_down_cost'],
                'start_up_cost': n_current.links.at[smelter, 'start_up_cost'],
                'committable': n_current.links.at[smelter, 'committable'],
                'p_min_pu': n_current.links.at[smelter, 'p_min_pu']
            }
        
        # 移除现有的smelter
        n_current.mremove("Link", aluminum_smelters)
        
        # 重新添加smelter，在步骤1中禁用电解铝启停约束
        for smelter_name, params in smelter_params.items():
            n_current.add("Link",
                smelter_name,
                bus0=params['bus0'],
                bus1=params['bus1'],
                carrier=params['carrier'],
                p_nom=params['p_nom'],
                p_nom_extendable=params['p_nom_extendable'],
                efficiency=params['efficiency'],
                marginal_cost=params['marginal_cost'],
                capital_cost=params['capital_cost'],
                start_up_cost=params['start_up_cost'],
                shut_down_cost=params['shut_down_cost'],
                stand_by_cost=params['stand_by_cost'],
                committable=False,  # 在步骤1中禁用电解铝启停约束
                p_min_pu=0  # 在步骤1中设置最小出力为0
            )
        
        # 临时修改配置，禁用电解铝启停约束
        original_commitment = config.get("aluminum_commitment", False)
        config["aluminum_commitment"] = False
        
        # 如果有固定的电解铝用能，需要添加约束
        if aluminum_usage is not None:
            # 获取电解铝厂运行参数
            from scripts.scenario_utils import get_aluminum_smelter_operational_params
            
            # 根据电解铝用能模式动态设置约束
            for smelter in aluminum_usage.columns:
                if smelter in n_current.links.index:
                    aluminum_power = aluminum_usage[smelter].values
                    
                    # 从冶炼设备名称中提取省份信息
                    # 冶炼设备名称格式: "Province aluminum smelter"
                    province = smelter.replace(" aluminum smelter", "")
                    
                    # 获取该省份的运行参数
                    if province in al_smelter_p_nom.index:
                        operational_params = get_aluminum_smelter_operational_params(config, al_smelter_p_nom=al_smelter_p_nom[province])
                        min_pu = operational_params['p_min_pu']
                    else:
                        # 如果找不到省份，使用默认值
                        min_pu = config['aluminum'].get('al_p_min_pu', 0.7)
                                     
                    # 根据用能模式设置约束
                    for i, power in enumerate(aluminum_power):
                        if power < 1:# threshold set by me
                            # 当用能为0时，固定p_max_pu为0，p_min_pu为0
                            n_current.links_t.p_max_pu.at[n_current.snapshots[i], smelter] = 0
                            n_current.links_t.p_min_pu.at[n_current.snapshots[i], smelter] = 0
                        else:
                            # 当用能不为0时，不固定p_max_pu，但设置p_min_pu为最小出力比例
                            n_current.links_t.p_max_pu.at[n_current.snapshots[i], smelter] = 1
                            n_current.links_t.p_min_pu.at[n_current.snapshots[i], smelter] = min_pu
        
        # 求解网络
        # 只保留PyPSA支持的参数
        optimize_kwargs = {k: v for k, v in kwargs.items() if k in ALLOWED_OPTIMIZE_KWARGS}
        
        current_solver_options = solver_options
        
        if skip_iterations:
            status, condition = n_current.optimize(
                solver_name=solver_name,
                solver_options=current_solver_options,  # 使用solver_options参数传递求解器特定选项
                extra_functionality=extra_functionality,
                **optimize_kwargs,
            )
        else:
            status, condition = n_current.optimize.optimize_transmission_expansion_iteratively(
                solver_name=solver_name,
                solver_options=current_solver_options,  # 使用solver_options参数传递求解器特定选项
                track_iterations=track_iterations,
                min_iterations=min_iterations,
                max_iterations=max_transmission_iterations,
                extra_functionality=extra_functionality,
                **optimize_kwargs,
            )
        
        # 记录求解性能统计
        solver_performance['iteration_times'].append(time.time() - iteration_start_time)
        
        # 获取当前目标函数值
        current_objective = None
        if hasattr(n_current, 'objective'):
            current_objective = n_current.objective
        elif hasattr(n_current, 'model') and hasattr(n_current.model, 'objective_value'):
            current_objective = n_current.model.objective_value
        
        # 提取当前迭代的节点电价（用于电解铝优化）
        current_nodal_prices = None
        if hasattr(n_current, 'buses_t') and hasattr(n_current.buses_t, 'marginal_price'):
            # 获取所有carrier为"AC"的电力节点的边际电价
            electricity_buses = n_current.buses[n_current.buses.carrier == "AC"].index
            if len(electricity_buses) > 0:
                # 使用所有AC节点的边际电价
                current_nodal_prices = n_current.buses_t.marginal_price[electricity_buses]
        
        # 步骤2: 检查收敛性 - 基于步骤一的目标函数变化的相对值
        if previous_objective is not None and current_objective is not None:
            # 计算目标函数变化的相对值
            objective_change = abs(current_objective - previous_objective)
            relative_change = objective_change / abs(previous_objective) if abs(previous_objective) > 1e-10 else float('inf')
            
            if relative_change < convergence_tolerance:
                final_network = n_current
                break
        else:
            # 第一次迭代，保存目标函数值用于下次比较
            if current_objective is not None:
                previous_objective = current_objective
        
        # 恢复电解铝启停约束
        config["aluminum_commitment"] = True
        
        # 读取铝冶炼厂年产量数据，获取需要优化的省份列表
        if 'snakemake' in globals():
            al_smelter_p_nom_path = snakemake.input.al_smelter_p_max
        else:
            al_smelter_p_nom_path = "data/p_nom/al_smelter_p_max.csv"
        
        al_smelter_annual_production = pd.read_csv(al_smelter_p_nom_path)
        al_smelter_annual_production = al_smelter_annual_production.set_index('Province')['p_nom']
        al_smelter_annual_production = al_smelter_annual_production[al_smelter_annual_production > 0.01]
        
        # Convert annual production (10kt/year) to power capacity (MW)
        al_smelter_p_nom = al_smelter_annual_production * 10000 * 13.3 / 8760  # Convert to MW
        
        # 检查哪些省份在网络中实际存在电解铝组件
        available_provinces = []
        for province in al_smelter_p_nom.index:
            # 检查该省份是否有电解铝冶炼设备
            aluminum_smelters = n_current.links[n_current.links.carrier == "aluminum"].index
            # 过滤出该省份的电解槽
            province_smelters = [smelter for smelter in aluminum_smelters if province in smelter]
            
            if province_smelters:
                available_provinces.append(province)
        
        target_provinces = available_provinces
        
        if not target_provinces:
            break
        
        # 获取全国优化结果中各省份的铝冶炼厂产量平均值
        national_smelter_production = {}
        # 获取所有铝冶炼厂的产量时间序列
        all_aluminum_smelters = n_current.links[n_current.links.carrier == "aluminum"].index
        
        for province in target_provinces:
            # 找到该省份的电解槽
            province_smelters = [smelter for smelter in all_aluminum_smelters if province in smelter]
            if province_smelters:
                # 获取该省份所有电解槽的产量时间序列
                province_smelter_production = n_current.links_t.p1[province_smelters]
                # 计算平均产量（取绝对值，因为产量通常为负值）
                average_production = province_smelter_production.abs().mean().sum()
                national_smelter_production[province] = average_production if average_production >= 1 else 0
        
        # 求解多个省份的电解铝优化问题（支持并行和串行）
        # 获取并行计算参数
        max_workers = kwargs.get("max_workers", None)  # 可以从kwargs中获取最大进程数
        overrides_path = kwargs.get("overrides_path", "data/override_component_attrs")
        use_parallel = kwargs.get("max_workers") is not None  # 如果指定了max_workers，则使用并行
        
        if use_parallel:
            # 使用并行函数求解
            all_aluminum_usage = solve_aluminum_optimization_parallel(
                target_provinces=target_provinces,
                original_network_path=original_network_path,
                original_overrides=original_overrides,
                config=config,
                solving=solving,
                opts=opts,
                current_nodal_prices=current_nodal_prices,
                solve_opts=kwargs.get("solve_opts", {}),
                using_single_node=using_single_node,
                single_node_province=single_node_province,
                overrides_path=overrides_path,
                max_workers=max_workers,
                national_smelter_production=national_smelter_production,
                **kwargs
            )
        else:
            # 使用原来的串行方法
            all_aluminum_usage = {}
            
            for province in target_provinces:
                # 为每个省份创建网络副本
                if original_network_path:
                    if original_overrides:
                        n_province = pypsa.Network(original_network_path, override_component_attrs=original_overrides)
                    else:
                        n_province = pypsa.Network(original_network_path)
                    
                    # 设置网络文件路径
                    n_province._network_path = original_network_path
                    if original_overrides:
                        n_province._overrides_path = overrides_path
                
                # 重新应用网络准备
                n_province = prepare_network(
                    n_province,
                    kwargs.get("solve_opts", {}),
                    using_single_node=using_single_node,
                    single_node_province=single_node_province
                )
                
                # 设置配置
                n_province.config = config
                n_province.opts = opts
                
                # 求解单个省份的电解铝优化问题
                province_aluminum_usage = solve_aluminum_optimization(
                    n_province, 
                    config, 
                    solving, 
                    opts, 
                    nodal_prices=current_nodal_prices, 
                    target_province=province,
                    national_smelter_production=national_smelter_production,
                    **kwargs
                )
                
                if province_aluminum_usage is not None:
                    all_aluminum_usage[province] = province_aluminum_usage
        
        if not all_aluminum_usage:
            break
        
        # 合并所有省份的电解铝用能结果
        # 获取所有冶炼设备名称
        all_smelters = []
        for province_usage in all_aluminum_usage.values():
            all_smelters.extend(province_usage.columns.tolist())
        
        # 创建合并后的DataFrame
        merged_aluminum_usage = pd.DataFrame(index=current_nodal_prices.index, columns=all_smelters)
        merged_aluminum_usage = merged_aluminum_usage.fillna(0).infer_objects(copy=False)
        
        # 填充各省份的数据
        for province, province_usage in all_aluminum_usage.items():
            for smelter in province_usage.columns:
                if smelter in merged_aluminum_usage.columns:
                    merged_aluminum_usage[smelter] = province_usage[smelter]
        
        # 更新电解铝用能和目标函数值
        aluminum_usage = merged_aluminum_usage
        previous_objective = current_objective
        final_network = n_current
        
        # 记录本次迭代时间
        iteration_time = time.time() - iteration_start_time
        iteration_times.append(iteration_time)
    
    # 恢复原始配置
    config["aluminum_commitment"] = original_commitment
    

    
    # 返回最终的网络结果
    return final_network

def solve_network_standard(n, config, solving, opts="", **kwargs):
    """
    标准网络求解方法（非迭代）
    """
    set_of_options = solving["solver"]["options"]
    solver_options = solving["solver_options"][set_of_options] if set_of_options else {}
    solver_name = solving["solver"]["name"]
    cf_solving = solving["options"]
    track_iterations = cf_solving.get("track_iterations", False)
    min_iterations = cf_solving.get("min_iterations", 4)
    max_iterations = cf_solving.get("max_iterations", 6)

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    skip_iterations = cf_solving.get("skip_iterations", False)
    if not n.lines.s_nom_extendable.any():
        skip_iterations = True
    
    # 使用标准参数求解
    # 只保留PyPSA支持的参数
    optimize_kwargs = {k: v for k, v in kwargs.items() if k in ALLOWED_OPTIMIZE_KWARGS}
    if skip_iterations:
        status, condition = n.optimize(
            solver_name=solver_name,
            solver_options=solver_options,  # 使用solver_options参数传递求解器特定选项
            extra_functionality=extra_functionality,
            **optimize_kwargs,
        )
    else:
        status, condition = n.optimize.optimize_transmission_expansion_iteratively(
            solver_name=solver_name,
            solver_options=solver_options,  # 使用solver_options参数传递求解器特定选项
            track_iterations=track_iterations,
            min_iterations=min_iterations,
            max_iterations=max_iterations,
            extra_functionality=extra_functionality,
            **optimize_kwargs,
        )

    # Store the objective value from the model (兼容性处理)
    if hasattr(n.model, 'objective_value'):
        n.objective = n.model.objective_value
    elif hasattr(n.model, 'objective'):
        n.objective = n.model.objective.value

    return n

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('solve_network_myopic-1',
                                   opts='ll',
                                   topology='current+Neighbor',
                                   pathway='linear2050',
                                   co2_reduction='0.0',
                                   planning_horizons="2050")

    configure_logging(snakemake)

    opts = snakemake.wildcards.opts
    if "sector_opts" in snakemake.wildcards.keys():
        opts += "-" + snakemake.wildcards.sector_opts
    opts = [o for o in opts.split("-") if o != ""]
    solve_opts = snakemake.params.solving["options"]

    np.random.seed(solve_opts.get("seed", 123))

    if "overrides" in snakemake.input.keys():
        overrides = override_component_attrs(snakemake.input.overrides)
        n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)
    else:
        n = pypsa.Network(snakemake.input.network)
    
    # 设置网络文件路径，用于迭代优化中的网络重新加载
    n._network_path = snakemake.input.network
    if "overrides" in snakemake.input.keys():
        n._overrides_path = snakemake.input.overrides

    n = prepare_network(
        n,
        solve_opts,
        using_single_node=snakemake.params.using_single_node,
        single_node_province=snakemake.params.single_node_province
    )

    # 检查是否启用电解铝迭代优化
    # 条件1: 检查snakemake.params.iterative_optimization
    # 条件2: 检查config["aluminum"]["grid_interaction"][planning_horizons]
    # 条件3: 检查config.get("add_aluminum", False)
    
    planning_horizons = snakemake.wildcards.planning_horizons
    
    # 检查电解铝过剩率条件
    aluminum_grid_interaction_condition = False
    if (snakemake.config.get("aluminum", {}).get("grid_interaction", {}).get(planning_horizons, False)):
        aluminum_grid_interaction_condition = True
    
    # 检查电解铝功能启用条件
    aluminum_enabled_condition = snakemake.config.get("add_aluminum", False)
    
    # 综合判断是否启用电解铝迭代优化
    if (snakemake.params.iterative_optimization and 
        aluminum_grid_interaction_condition and 
        aluminum_enabled_condition):
        # 获取迭代优化参数
        max_iterations = snakemake.config.get("aluminum_max_iterations", 10)
        convergence_tolerance = snakemake.config.get("aluminum_convergence_tolerance", 0.01)
        
        # 准备传递给迭代函数的参数
        iteration_kwargs = {
            "log_fn": snakemake.log.solver,
            "solve_opts": solve_opts,
            "using_single_node": snakemake.params.using_single_node,
            "single_node_province": snakemake.params.single_node_province
        }
        
        # 如果有overrides，也传递
        if "overrides" in snakemake.input.keys():
            iteration_kwargs["overrides"] = overrides
        
        # 添加并行计算配置
        if snakemake.config.get("aluminum_parallel", True):  # 默认启用并行计算
            max_workers = snakemake.config.get("aluminum_max_workers", None)  # 默认使用CPU核心数
            iteration_kwargs["max_workers"] = max_workers
        
        # 使用迭代优化算法
        n = solve_network_iterative(
            n,
            config=snakemake.config,
            solving=snakemake.params.solving,
            opts=opts,
            max_iterations=max_iterations,
            convergence_tolerance=convergence_tolerance,
            **iteration_kwargs,
        )
    else:
        # 使用标准求解方法
        n = solve_network_standard(
            n,
            config=snakemake.config,
            solving=snakemake.params.solving,
            opts=opts,
            log_fn=snakemake.log.solver,
        )

    # 确保输出目录存在
    output_dir = os.path.dirname(snakemake.output.network_name)
    os.makedirs(output_dir, exist_ok=True)

    #n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.links_t.p2 = n.links_t.p2.astype(float)
    # Clean up links_t_p3 data before export to avoid dtype issues
    if hasattr(n, 'links_t') and hasattr(n.links_t, 'p3'):
        # Convert DataFrame to numeric, handling any non-numeric values
        n.links_t.p3 = n.links_t.p3.apply(pd.to_numeric, errors='coerce').fillna(0.0).infer_objects(copy=False)
    
    # 导出结果
    n.export_to_netcdf(snakemake.output.network_name) 