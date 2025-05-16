# SPDX-FileCopyrightText: : 2024 The PyPSA-China Authors
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
):

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

def add_battery_constraints(n):
    """
    Add constraint ensuring that charger = discharger, i.e.
    1 * charger_size - efficiency * discharger_size = 0
    """
    if not n.links.p_nom_extendable.any():
        return

    discharger_bool = n.links.index.str.contains("battery discharger")
    charger_bool = n.links.index.str.contains("battery charger")

    dischargers_ext = n.links[discharger_bool].query("p_nom_extendable").index
    chargers_ext = n.links[charger_bool].query("p_nom_extendable").index

    eff = n.links.efficiency[dischargers_ext].values
    lhs = (
        n.model["Link-p_nom"].loc[chargers_ext]
        - n.model["Link-p_nom"].loc[dischargers_ext] * eff
    )

    n.model.add_constraints(lhs == 0, name="Link-charger_ratio")

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

        lhs = (
            p_nom.loc[electric_ext]
            * (n.links.p_nom_ratio * n.links.efficiency)[electric_ext].values
            - p_nom.loc[heat_ext] * n.links.efficiency[heat_ext].values
        )
        n.model.add_constraints(lhs == 0, name="chplink-fix_p_nom_ratio")

        rename = {"Link-ext": "Link"}
        lhs = (
            p.loc[:, electric_ext]
            + p.loc[:, heat_ext]
            - p_nom.rename(rename).loc[electric_ext]
        )
        n.model.add_constraints(lhs <= 0, name="chplink-top_iso_fuel_line_ext")

    # top_iso_fuel_line for fixed
    if not electric_fix.empty:
        lhs = p.loc[:, electric_fix] + p.loc[:, heat_fix]
        rhs = n.links.p_nom[electric_fix]
        n.model.add_constraints(lhs <= rhs, name="chplink-top_iso_fuel_line_fix")

    # back-pressure
    if not n.links[electric].index.empty:
        lhs = (
            p.loc[:, heat] * (n.links.efficiency[heat] * n.links.c_b[electric].values)
            - p.loc[:, electric] * n.links.efficiency[electric]
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

        lhs  = (
            n.model["Generator-p_nom"].loc[coal]
            + n.model["Generator-p_nom"].loc[coal_retrofit]
            - (p_nom_max[str(year)].loc[Bus] - coal_retrofitted.reindex(p_nom_max[str(year)].loc[Bus].index,fill_value=0)).values
        )

        n.model.add_constraints(lhs == 0, name="Generator-coal-retrofit-" + str(year))

def extra_functionality(n, snapshots):
    """
    Collects supplementary constraints which will be passed to ``pypsa.linopf.network_lopf``.
    If you want to enforce additional custom constraints, this is a good location to add them.
    The arguments ``opts`` and ``snakemake.config`` are expected to be attached to the network.
    """
    opts = n.opts
    config = n.config
    add_chp_constraints(n)
    add_battery_constraints(n)
    add_transimission_constraints(n)
    if snakemake.wildcards.planning_horizons != "2020":
        add_retrofit_constraints(n)


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
    
    if skip_iterations:
        status, condition = n.optimize(
            solver_name=solver_name,
            extra_functionality=extra_functionality,
            **solver_options,
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

    if status != "ok":
        logger.warning(
            f"Solving status '{status}' with termination condition '{condition}'"
        )
    if "infeasible" in condition:
        raise RuntimeError("Solving status 'infeasible'")

    # Store the objective value from the model
    if hasattr(n.model, 'objective_value'):
        n.objective = n.model.objective_value
    elif hasattr(n.model, 'objective'):
        n.objective = n.model.objective.value
    else:
        logger.warning("Could not find objective value in model")

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
        solve_opts
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
    n.export_to_netcdf(snakemake.output.network_name)