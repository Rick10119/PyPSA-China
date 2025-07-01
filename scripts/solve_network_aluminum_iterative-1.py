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
import yaml
import os
from _helpers import (
    configure_logging,
    override_component_attrs,
)

# 允许传递给PyPSA optimize的参数
ALLOWED_OPTIMIZE_KWARGS = [
    "solver_name", "solver_options", "extra_functionality",
    "track_iterations", "min_iterations", "max_iterations"
]

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)
from pypsa.descriptors import get_switchable_as_dense as get_as_dense

def prepare_network(
        n,
        solve_opts=None,
        single_node=False,
        single_node_province="Shandong"
):
    # 检查是否启用单节点模式
    if single_node:
        logger.info(f"启用单节点模式，只保留 {single_node_province} 省份的组件")
        
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
        
        logger.info(f"单节点过滤完成，保留了 {len(province_buses)} 个节点")
    
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
    # add_chp_constraints(n)
    add_transimission_constraints(n)
    if snakemake.wildcards.planning_horizons != "2020":
        add_retrofit_constraints(n)

def solve_aluminum_optimization(n, config, solving, opts="", nodal_prices=None, **kwargs):
    """
    电解铝最优运行问题求解
    基于节点电价，运行以满足铝需求为约束的电解铝最优运行问题
    """
    set_of_options = solving["solver"]["options"]
    solver_options = solving["solver_options"][set_of_options] if set_of_options else {}
    solver_name = solving["solver"]["name"]
    
    # 不复制网络，而是直接使用传入的网络
    # 找到所有电解铝相关的组件
    aluminum_buses = n.buses[n.buses.carrier == "aluminum"].index
    aluminum_smelters = n.links[n.links.carrier == "aluminum smelter"].index
    aluminum_stores = n.stores[n.stores.carrier == "aluminum storage"].index
    aluminum_loads = n.loads[n.loads.bus.isin(aluminum_buses)].index
    
    # 重新设置电解铝冶炼设备的参数
    for smelter in aluminum_smelters:
        n.links.at[smelter, 'committable'] = config['aluminum_commitment']
        n.links.at[smelter, 'p_min_pu'] = config['aluminum']['al_p_min_pu'] if config['aluminum_commitment'] else 0
    
    # 移除所有非电解铝相关的组件
    for component_type in ["Generator", "StorageUnit", "Store", "Link", "Load"]:
        if component_type == "Store":
            # 保留电解铝存储
            other_stores = n.stores[~n.stores.index.isin(aluminum_stores)].index
            n.mremove(component_type, other_stores)
        elif component_type == "Link":
            # 保留电解铝冶炼设备
            other_links = n.links[~n.links.index.isin(aluminum_smelters)].index
            n.mremove(component_type, other_links)
        elif component_type == "Load":
            # 保留电解铝负荷
            other_loads = n.loads[~n.loads.index.isin(aluminum_loads)].index
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
            logger.info(f"为节点 {bus} 使用其对应的边际电价")
            
            n.add("Generator",
                f"virtual_gen_{bus}",
                bus=bus,
                carrier="virtual",
                p_nom=1e6,  # 大容量
                marginal_cost=marginal_cost)  # 使用节点边际电价作为边际成本
            
    # 打印所有bus:
    logger.info(f"所有bus: {n.buses.index}")
    
    # 求解电解铝优化问题
    try:
        # 只保留PyPSA支持的参数
        optimize_kwargs = {k: v for k, v in kwargs.items() if k in ALLOWED_OPTIMIZE_KWARGS}
        n.optimize(
            solver_name=solver_name,
            extra_functionality=extra_functionality,
            **solver_options,
            **optimize_kwargs,
        )
        
        logger.info("电解铝优化问题求解成功")
        # 提取电解铝用能模式
        aluminum_usage = n.links_t.p0[aluminum_smelters].copy()
        return aluminum_usage
            
    except Exception as e:
        logger.error(f"电解铝优化问题求解出错: {e}")
        return None

def solve_network_iterative(n, config, solving, opts="", max_iterations=10, convergence_tolerance=0.01, **kwargs):
    """
    电解铝迭代优化算法
    1. 使用连续化电解铝模型求解，得到节点电价
    2. 基于节点电价，运行电解铝最优运行问题
    3. 固定电解铝用能，求解剩下的优化问题
    4. 重复2-3，直到收敛
    
    收敛条件：目标函数变化的相对值小于convergence_tolerance（默认1%）
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
        logger.info("No expandable lines found. Skipping iterative solving.")
    
    # 检查是否启用电解铝
    if not config.get("add_aluminum", False):
        logger.info("电解铝功能未启用，使用标准求解方法")
        return solve_network_standard(n, config, solving, opts, **kwargs)
    
    logger.info("开始电解铝迭代优化算法")
    
    # 记录总开始时间
    total_start_time = time.time()
    
    # 初始化电解铝用能模式和目标函数值
    aluminum_usage = None
    previous_objective = None
    iteration = 0
    final_network = None
    
    # 记录每次迭代的时间
    iteration_times = []
    
    # 保存原始网络文件路径，用于重新加载
    original_network_path = None
    if hasattr(n, '_network_path'):
        original_network_path = n._network_path
    else:
        # 如果没有网络路径，保存网络对象用于复制
        original_network = n
    
    while iteration < max_iterations:
        iteration += 1
        iteration_start_time = time.time()  # 记录每次迭代开始时间
        
        logger.info(f"开始第 {iteration} 次迭代")
        
        # 步骤1: 使用连续化电解铝模型求解，得到节点电价
        logger.info("步骤1: 使用连续化电解铝模型求解")
        
        # 重新加载网络，禁用电解铝启停约束
        if original_network_path:
            # 从文件重新加载网络
            if "overrides" in kwargs:
                overrides = kwargs["overrides"]
                n_current = pypsa.Network(original_network_path, override_component_attrs=overrides)
            else:
                n_current = pypsa.Network(original_network_path)
        else:
            # 如果没有网络路径，直接使用原始网络（不复制）
            n_current = n
        
        # 重新应用网络准备
        n_current = prepare_network(
            n_current,
            kwargs.get("solve_opts", {}),
            kwargs.get("single_node", False),
            kwargs.get("single_node_province", "Shandong")
        )
        
        # 设置配置
        n_current.config = config
        n_current.opts = opts
        
        # 重新设置电解铝冶炼设备的参数
        aluminum_smelters = n_current.links[n_current.links.carrier == "aluminum smelter"].index
        for smelter in aluminum_smelters:
            n_current.links.at[smelter, 'committable'] = config['aluminum_commitment']
            n_current.links.at[smelter, 'p_min_pu'] = config['aluminum']['al_p_min_pu'] if config['aluminum_commitment'] else 0
        
        # 临时修改配置，禁用电解铝启停约束
        original_commitment = config.get("aluminum_commitment", False)
        config["aluminum_commitment"] = False
        
        # 如果有固定的电解铝用能，需要添加约束
        if aluminum_usage is not None:
            logger.info("固定电解铝用能，求解剩余优化问题")
            # 固定电解铝用能 - 使用p_set来固定smelter的link出力
            for smelter in aluminum_usage.columns:
                if smelter in n_current.links.index:
                    fixed_aluminum_power = aluminum_usage[smelter].values
                    # 确保p_set存在
                    if not hasattr(n_current.links_t, 'p_set'):
                        n_current.links_t.p_set = pd.DataFrame(index=n_current.snapshots, columns=n_current.links.index)
                    # 修改smelter的link出力
                    n_current.links_t.p_set[smelter] = fixed_aluminum_power
                    logger.info(f"固定电解铝冶炼设备 {smelter} 的用能模式")
                    
                    # 同时修改对应的负荷设定
                    aluminum_loads = n_current.loads[n_current.loads.bus.isin(n_current.buses[n_current.buses.carrier == "aluminum"].index)].index
                    for load in aluminum_loads:
                        if hasattr(n_current.loads_t, 'p_set'):
                            n_current.loads_t.p_set[load] = fixed_aluminum_power
                            logger.info(f"固定电解铝负荷 {load} 的用能模式")
        
        # 求解网络
        try:
            # 只保留PyPSA支持的参数
            optimize_kwargs = {k: v for k, v in kwargs.items() if k in ALLOWED_OPTIMIZE_KWARGS}
            if skip_iterations:
                n_current.optimize(
                    solver_name=solver_name,
                    extra_functionality=extra_functionality,
                    **solver_options,
                    **optimize_kwargs,
                )
            else:
                n_current.optimize.optimize_transmission_expansion_iteratively(
                    solver_name=solver_name,
                    track_iterations=track_iterations,
                    min_iterations=min_iterations,
                    max_iterations=max_transmission_iterations,
                    extra_functionality=extra_functionality,
                    **solver_options,
                    **optimize_kwargs,
                )
            
            logger.info("网络求解成功，获得节点电价")
            
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
                    logger.info(f"提取节点电价: 共 {len(electricity_buses)} 个AC节点")
                    logger.info(f"节点列表: {list(electricity_buses)}")
                else:
                    logger.warning("未找到carrier为AC的节点，无法提取节点电价")
            
        except Exception as e:
            logger.error(f"网络求解出错: {e}")
            break
        
        # 步骤2: 基于节点电价，运行电解铝最优运行问题
        logger.info("步骤2: 运行电解铝最优运行问题")
        
        # 恢复电解铝启停约束
        config["aluminum_commitment"] = True
        
        # 重新设置电解铝冶炼设备的参数（启用启停约束）
        aluminum_smelters = n_current.links[n_current.links.carrier == "aluminum smelter"].index
        for smelter in aluminum_smelters:
            n_current.links.at[smelter, 'committable'] = config['aluminum_commitment']
            n_current.links.at[smelter, 'p_min_pu'] = config['aluminum']['al_p_min_pu'] if config['aluminum_commitment'] else 0
        
        # 求解电解铝优化问题，传递节点电价
        new_aluminum_usage = solve_aluminum_optimization(n_current, config, solving, opts, nodal_prices=current_nodal_prices, **kwargs)
        
        if new_aluminum_usage is None:
            logger.error("电解铝优化问题求解失败")
            break
        
        # 步骤3: 检查收敛性 - 基于目标函数变化的相对值
        if previous_objective is not None and current_objective is not None:
            # 计算目标函数变化的相对值
            objective_change = abs(current_objective - previous_objective)
            relative_change = objective_change / abs(previous_objective) if abs(previous_objective) > 1e-10 else float('inf')
            
            logger.info(f"目标函数变化统计:")
            logger.info(f"  当前目标函数值: {current_objective:.6e}")
            logger.info(f"  上次目标函数值: {previous_objective:.6e}")
            logger.info(f"  绝对变化: {objective_change:.6e}")
            logger.info(f"  相对变化: {relative_change:.6f} ({relative_change*100:.2f}%)")
            logger.info(f"  收敛阈值: {convergence_tolerance:.6f} ({convergence_tolerance*100:.2f}%)")
            
            if relative_change < convergence_tolerance:
                logger.info(f"算法收敛，在第 {iteration} 次迭代后停止")
                logger.info(f"目标函数相对变化 {relative_change:.6f} ({relative_change*100:.2f}%) < 收敛阈值 {convergence_tolerance:.6f} ({convergence_tolerance*100:.2f}%)")
                final_network = n_current
                break
        else:
            # 第一次迭代，保存目标函数值用于下次比较
            if current_objective is not None:
                previous_objective = current_objective
        
        # 更新电解铝用能和目标函数值
        aluminum_usage = new_aluminum_usage
        previous_objective = current_objective
        final_network = n_current
        
        # 记录本次迭代时间
        iteration_time = time.time() - iteration_start_time
        iteration_times.append(iteration_time)
        
        logger.info(f"第 {iteration} 次迭代完成，耗时: {iteration_time:.2f} 秒")
    
    # 恢复原始配置
    config["aluminum_commitment"] = original_commitment
    
    # 计算总时间
    total_time = time.time() - total_start_time
    
    if iteration >= max_iterations:
        logger.warning(f"达到最大迭代次数 {max_iterations}，算法未完全收敛")
    
    # 输出时间统计
    logger.info(f"\n=== 迭代时间统计 ===")
    logger.info(f"总迭代次数: {iteration}")
    logger.info(f"总耗时: {total_time:.2f} 秒 ({total_time/60:.2f} 分钟)")
    if iteration_times:
        logger.info(f"平均每次迭代耗时: {np.mean(iteration_times):.2f} 秒")
        logger.info(f"最快迭代耗时: {np.min(iteration_times):.2f} 秒")
        logger.info(f"最慢迭代耗时: {np.max(iteration_times):.2f} 秒")
        logger.info(f"迭代时间详情:")
        for i, t in enumerate(iteration_times, 1):
            logger.info(f"  第{i}次迭代: {t:.2f} 秒")
    
    logger.info(f"电解铝迭代优化算法完成，共进行 {iteration} 次迭代")
    
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
        logger.info("No expandable lines found. Skipping iterative solving.")
    
    # 使用标准参数求解
    try:
        # 只保留PyPSA支持的参数
        optimize_kwargs = {k: v for k, v in kwargs.items() if k in ALLOWED_OPTIMIZE_KWARGS}
        if skip_iterations:
            n.optimize(
                solver_name=solver_name,
                extra_functionality=extra_functionality,
                **solver_options,
                **optimize_kwargs,
            )
        else:
            n.optimize.optimize_transmission_expansion_iteratively(
                solver_name=solver_name,
                track_iterations=track_iterations,
                min_iterations=min_iterations,
                max_iterations=max_iterations,
                extra_functionality=extra_functionality,
                **solver_options,
                **optimize_kwargs,
            )

        logger.info("标准求解成功")

    except Exception as e:
        logger.error(f"标准求解失败: {e}")
        raise e

    # Store the objective value from the model (兼容性处理)
    try:
        if hasattr(n.model, 'objective_value'):
            n.objective = n.model.objective_value
        elif hasattr(n.model, 'objective'):
            n.objective = n.model.objective.value
        else:
            logger.warning("Could not find objective value in model")
    except Exception as e:
        logger.warning(f"Error accessing objective value: {e}")

    return n

if __name__ == '__main__':
    # 加载配置文件
    config_path = "config.yaml"
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"配置文件 {config_path} 不存在")
    
    with open(config_path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    # 设置默认参数
    default_params = {
        'opts': 'll',
        'topology': 'current+Neighbor',
        'pathway': 'linear2050',
        'co2_reduction': '0.0',
        'planning_horizons': '2050',
        'heating_demand': 'positive'
    }
    
    # 创建模拟的snakemake对象
    class MockSnakemake:
        def __init__(self, config, default_params):
            self.config = config
            self.params = {
                'solving': config['solving']
            }
            self.wildcards = type('Wildcards', (), default_params)()
            
            # 设置输入文件路径
            version = config['version']
            results_dir = config['results_dir']
            
            # 构建网络文件路径
            network_path = f"{results_dir}version-{version}/prenetworks-brownfield/{default_params['heating_demand']}/prenetwork-{default_params['opts']}-{default_params['topology']}-{default_params['pathway']}-{default_params['planning_horizons']}.nc"
            
            self.input = type('Input', (), {
                'network': network_path,
                'overrides': 'data/override_component_attrs'
            })()
            
            # 设置输出文件路径
            output_path = f"{results_dir}version-{version}/postnetworks/{default_params['heating_demand']}/postnetwork-{default_params['opts']}-{default_params['topology']}-{default_params['pathway']}-{default_params['planning_horizons']}.nc"
            
            self.output = type('Output', (), {
                'network_name': output_path
            })()
            
            # 设置日志文件路径
            log_dir = f"logs/solve_operations_network/{default_params['heating_demand']}"
            os.makedirs(log_dir, exist_ok=True)
            log_path = f"{log_dir}/postnetwork-{default_params['opts']}-{default_params['topology']}-{default_params['pathway']}-{default_params['planning_horizons']}.log"
            
            self.log = type('Log', (), {
                'solver': log_path
            })()
    
    # 创建模拟的snakemake对象
    snakemake = MockSnakemake(config, default_params)
    
    # 配置日志 - 简化版本，不依赖snakemake.rule
    import logging
    logging.basicConfig(
        level=config.get('logging', {}).get('level', 'INFO'),
        format='%(levelname)s:%(name)s:%(message)s',
        handlers=[
            logging.FileHandler(snakemake.log.solver),
            logging.StreamHandler()
        ]
    )

    opts = snakemake.wildcards.opts
    if "sector_opts" in snakemake.wildcards.__dict__.keys():
        opts += "-" + snakemake.wildcards.sector_opts
    opts = [o for o in opts.split("-") if o != ""]
    solve_opts = snakemake.params['solving']["options"]

    np.random.seed(solve_opts.get("seed", 123))

    # 检查输入文件是否存在
    if not os.path.exists(snakemake.input.network):
        raise FileNotFoundError(f"网络文件不存在: {snakemake.input.network}")
    
    if "overrides" in snakemake.input.__dict__.keys():
        overrides = override_component_attrs(snakemake.input.overrides)
        n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)
    else:
        n = pypsa.Network(snakemake.input.network)

    n = prepare_network(
        n,
        solve_opts,
        single_node=snakemake.config.get("single_node", False),
        single_node_province=snakemake.config.get("single_node_province", "Shandong")
    )

    # 检查是否启用电解铝迭代优化
    if snakemake.config.get("add_aluminum", False):
        # 获取迭代优化参数
        max_iterations = snakemake.config.get("aluminum_max_iterations", 10)
        convergence_tolerance = snakemake.config.get("aluminum_convergence_tolerance", 0.01)
        
        logger.info(f"启用电解铝迭代优化算法")
        logger.info(f"最大迭代次数: {max_iterations}")
        logger.info(f"收敛阈值: {convergence_tolerance}")
        
        # 准备传递给迭代函数的参数
        iteration_kwargs = {
            "log_fn": snakemake.log.solver,
            "solve_opts": solve_opts,
            "single_node": snakemake.config.get("single_node", False),
            "single_node_province": snakemake.config.get("single_node_province", "Shandong")
        }
        
        # 如果有overrides，也传递
        if "overrides" in snakemake.input.__dict__.keys():
            iteration_kwargs["overrides"] = overrides
        
        # 使用迭代优化算法
        n = solve_network_iterative(
            n,
            config=snakemake.config,
            solving=snakemake.params['solving'],
            opts=opts,
            max_iterations=max_iterations,
            convergence_tolerance=convergence_tolerance,
            **iteration_kwargs,
        )
    else:
        logger.info("使用标准求解方法")
        # 使用标准求解方法
        n = solve_network_standard(
            n,
            config=snakemake.config,
            solving=snakemake.params['solving'],
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
        n.links_t.p3 = n.links_t.p3.apply(pd.to_numeric, errors='coerce').fillna(0.0)
    
    # 导出结果
    n.export_to_netcdf(snakemake.output.network_name)
    logger.info(f"结果已保存到: {snakemake.output.network_name}") 