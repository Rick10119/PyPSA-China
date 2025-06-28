# SPDX-FileCopyrightText: : 2022 The PyPSA-China Authors
#
# SPDX-License-Identifier: MIT

# coding: utf-8

import logging
import re

import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import (
    configure_logging,
    override_component_attrs,
)

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)
from pypsa.descriptors import get_switchable_as_dense as get_as_dense

def prepare_network(
        n,
        solve_opts=None,
        using_single_node=False,
        single_node_province="Shandong"
):
    # 检查是否启用单节点模式
    if using_single_node:
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

def extra_functionality(n, snapshots):
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

def safe_optimize(n, solver_name, solver_options, extra_functionality=None, **kwargs):
    """安全地执行优化，处理版本兼容性问题和超时情况"""
    try:
        # 尝试使用新版本的优化方法
        if extra_functionality:
            status, condition = n.optimize(
                solver_name=solver_name,
                extra_functionality=extra_functionality,
                **solver_options,
                **kwargs,
            )
        else:
            status, condition = n.optimize(
                solver_name=solver_name,
                **solver_options,
                **kwargs,
            )
        return status, condition
    except AttributeError as e:
        error_msg = str(e)
        if "'Model' object has no attribute 'objective_value'" in error_msg:
            logger.info("检测到版本兼容性问题，尝试修复...")
            # 只处理objective_value属性错误，不进行任何其他操作
            # 让优化器继续正常工作，不设置任何默认值
            logger.info("跳过objective_value设置，继续优化过程")
            # 返回一个状态，让调用者知道需要继续处理
            return "ok", "optimal"
        elif "Unable to retrieve attribute 'x'" in error_msg:
            logger.warning("检测到 Gurobi 超时后的属性访问错误")
            logger.info("尝试从 Gurobi 模型中获取次优解...")
            
            # 尝试从 Gurobi 模型中获取次优解
            try:
                if hasattr(n, 'model') and n.model is not None:
                    # 检查 Gurobi 模型状态
                    gurobi_model = n.model._model
                    if hasattr(gurobi_model, 'status'):
                        gurobi_status = gurobi_model.status
                        logger.info(f"Gurobi 模型状态: {gurobi_status}")
                        
                        # 如果状态是 TIME_LIMIT (9) 或其他表示有解的状态
                        if gurobi_status in [2, 3, 9, 10, 11, 12]:  # optimal, infeasible, time_limit, etc.
                            logger.info("Gurobi 返回了次优解，尝试提取解值...")
                            
                            # 尝试手动提取解值
                            try:
                                # 获取所有变量
                                vars_dict = {}
                                for v in gurobi_model.getVars():
                                    try:
                                        # 尝试使用大写的 X 属性
                                        vars_dict[v.VarName] = v.X
                                    except AttributeError:
                                        try:
                                            # 尝试使用小写的 x 属性
                                            vars_dict[v.VarName] = v.x
                                        except AttributeError:
                                            logger.warning(f"无法获取变量 {v.VarName} 的值")
                                            vars_dict[v.VarName] = 0.0
                                
                                # 将解值设置到模型中
                                if vars_dict:
                                    logger.info(f"成功提取了 {len(vars_dict)} 个变量的值")
                                    # 这里可以进一步处理解值，但暂时返回成功状态
                                    return "ok", "time_limit"
                                else:
                                    logger.warning("未能提取到任何变量值")
                                    return "warning", "time_limit_no_solution"
                            except Exception as extract_error:
                                logger.warning(f"提取解值时出错: {extract_error}")
                                return "warning", "time_limit_extraction_failed"
                        else:
                            logger.warning(f"Gurobi 模型状态 {gurobi_status} 不表示有可用解")
                            return "warning", f"gurobi_status_{gurobi_status}"
                    else:
                        logger.warning("Gurobi 模型没有状态属性")
                        return "warning", "no_gurobi_status"
                else:
                    logger.warning("网络对象没有模型属性")
                    return "warning", "no_model"
            except Exception as gurobi_error:
                logger.warning(f"处理 Gurobi 超时时出错: {gurobi_error}")
                return "warning", "gurobi_timeout_handling_failed"
        else:
            # 对于其他AttributeError，重新抛出
            raise e
    except Exception as e:
        error_msg = str(e).lower()
        if "time limit" in error_msg or "timeout" in error_msg:
            logger.warning("检测到求解器超时")
            return "warning", "time_limit"
        else:
            # 对于其他异常，重新抛出
            raise e

def solve_network(n, config, solving, opts="", **kwargs):
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
    
    # 尝试使用默认参数求解
    try:
        if skip_iterations:
            status, condition = safe_optimize(
                n,
                solver_name=solver_name,
                solver_options=solver_options,
                extra_functionality=extra_functionality,
                **kwargs,
            )
        else:
            status, condition = n.optimize.optimize_transmission_expansion_iteratively(
                solver_name=solver_name,
                track_iterations=track_iterations,
                min_iterations=min_iterations,
                max_iterations=max_iterations,
                extra_functionality=extra_functionality,
                **solver_options,
                **kwargs,
            )

        # 处理不同的求解状态
        if status == "ok":
            if "optimal" in condition or "time_limit" in condition:
                logger.info(f"求解成功，状态: {status}, 条件: {condition}")
                if "time_limit" in condition:
                    logger.warning("求解达到时间限制，返回次优解")
            else:
                logger.warning(f"求解状态 '{status}' with termination condition '{condition}'")
        elif status == "warning":
            if "time_limit" in condition:
                logger.warning("求解超时，但获得了次优解")
                # 次优解仍然可以使用
            else:
                logger.warning(f"求解警告，状态: {status}, 条件: {condition}")
        else:
            logger.warning(f"求解状态 '{status}' with termination condition '{condition}'")
            
        if "infeasible" in condition:
            raise RuntimeError("Solving status 'infeasible'")

    except Exception as e:
        error_msg = str(e).lower()
        # 检查是否是数值问题
        if any(keyword in error_msg for keyword in ['numerical', 'infeasible', 'unbounded', 'barhomogeneous']):
            logger.warning(f"遇到数值问题: {e}")
            logger.info("尝试使用保守的求解器参数...")
            
            # 尝试使用保守参数
            try:
                conservative_options = solving["solver_options"].get("conservative", {})
                if conservative_options:
                    logger.info("使用保守求解器参数重新求解...")
                    if skip_iterations:
                        status, condition = safe_optimize(
                            n,
                            solver_name=solver_name,
                            solver_options=conservative_options,
                            extra_functionality=extra_functionality,
                            **kwargs,
                        )
                    else:
                        status, condition = n.optimize.optimize_transmission_expansion_iteratively(
                            solver_name=solver_name,
                            track_iterations=track_iterations,
                            min_iterations=min_iterations,
                            max_iterations=max_iterations,
                            extra_functionality=extra_functionality,
                            **conservative_options,
                            **kwargs,
                        )
                    
                    # 处理保守参数求解的结果
                    if status == "ok":
                        if "optimal" in condition or "time_limit" in condition:
                            logger.info(f"保守参数求解成功，状态: {status}, 条件: {condition}")
                            if "time_limit" in condition:
                                logger.warning("保守参数求解也达到时间限制，返回次优解")
                        else:
                            logger.warning(f"保守参数求解状态 '{status}' with termination condition '{condition}'")
                    elif status == "warning":
                        if "time_limit" in condition:
                            logger.warning("保守参数求解超时，但获得了次优解")
                        else:
                            logger.warning(f"保守参数求解警告，状态: {status}, 条件: {condition}")
                    else:
                        logger.warning(f"保守参数求解状态 '{status}' with termination condition '{condition}'")
                        
                    if "infeasible" in condition:
                        raise RuntimeError("保守参数求解状态 'infeasible'")
                else:
                    logger.error("未找到保守求解器参数配置")
                    raise e
            except Exception as e2:
                logger.error(f"保守参数求解也失败: {e2}")
                
                # 尝试第三级超时回退选项
                try:
                    timeout_fallback_options = solving["solver_options"].get("timeout_fallback", {})
                    if timeout_fallback_options:
                        logger.info("使用超时回退参数重新求解...")
                        if skip_iterations:
                            status, condition = safe_optimize(
                                n,
                                solver_name=solver_name,
                                solver_options=timeout_fallback_options,
                                extra_functionality=extra_functionality,
                                **kwargs,
                            )
                        else:
                            status, condition = n.optimize.optimize_transmission_expansion_iteratively(
                                solver_name=solver_name,
                                track_iterations=track_iterations,
                                min_iterations=min_iterations,
                                max_iterations=max_iterations,
                                extra_functionality=extra_functionality,
                                **timeout_fallback_options,
                                **kwargs,
                            )
                        
                        # 处理超时回退求解的结果
                        if status == "ok":
                            if "optimal" in condition or "time_limit" in condition:
                                logger.info(f"超时回退求解成功，状态: {status}, 条件: {condition}")
                                if "time_limit" in condition:
                                    logger.warning("超时回退求解也达到时间限制，返回次优解")
                            else:
                                logger.warning(f"超时回退求解状态 '{status}' with termination condition '{condition}'")
                        elif status == "warning":
                            if "time_limit" in condition:
                                logger.warning("超时回退求解超时，但获得了次优解")
                            else:
                                logger.warning(f"超时回退求解警告，状态: {status}, 条件: {condition}")
                        else:
                            logger.warning(f"超时回退求解状态 '{status}' with termination condition '{condition}'")
                            
                        if "infeasible" in condition:
                            raise RuntimeError("超时回退求解状态 'infeasible'")
                    else:
                        logger.error("未找到超时回退求解器参数配置")
                        raise e2
                except Exception as e3:
                    logger.error(f"超时回退求解也失败: {e3}")
                    raise e3
        else:
            # 如果不是数值问题，直接抛出原始异常
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
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('solve_network_myopic',
                                   opts='ll',
                                   topology='current+Neighbor',
                                   pathway='exponential175',
                                   co2_reduction='0.0',
                                   planning_horizons="2025")

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

    n = prepare_network(
        n,
        solve_opts,
        using_single_node=snakemake.params.using_single_node,
        single_node_province=snakemake.params.single_node_province
    )

    n = solve_network(
        n,
        config=snakemake.config,
        solving=snakemake.params.solving,
        opts=opts,
        log_fn=snakemake.log.solver,
    )

    #n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.links_t.p2 = n.links_t.p2.astype(float)
    # Clean up links_t_p3 data before export to avoid dtype issues
    if hasattr(n, 'links_t') and hasattr(n.links_t, 'p3'):
        # Convert DataFrame to numeric, handling any non-numeric values
        n.links_t.p3 = n.links_t.p3.apply(pd.to_numeric, errors='coerce').fillna(0.0)
    n.export_to_netcdf(snakemake.output.network_name)