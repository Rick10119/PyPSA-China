import logging

# Set up a logger for this module
logger = logging.getLogger(__name__)

# Import required libraries for plotting and data handling
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
# Import helper functions from local modules
from make_summary import assign_carriers
from plot_summary import preferred_order, rename_techs
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches

# Set the plotting style for all matplotlib plots
# This function applies a consistent style to all plots for better aesthetics
# and readability

def set_plot_style():
    plt.style.use(['classic', 'seaborn-v0_8-white',
                   {'axes.grid': False, 'grid.linestyle': '--', 'grid.color': u'0.6',
                    'hatch.color': 'white',
                    'patch.linewidth': 0.5,
                    'font.size': 12,
                    'legend.fontsize': 'medium',
                    'lines.linewidth': 1.5,
                    'pdf.fonttype': 42,
                    }])

# Assigns a 'location' attribute to each network component based on its name.
# This is useful for grouping and plotting components by their physical location.
def assign_location(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)
        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue
            names = ifind.index[ifind == i]
            c.df.loc[names, "location"] = names.str[:i]

# Helper function to process and order cost data for plotting.
# It ensures costs are grouped, ordered, and filtered according to the preferred order and available tech colors.
def get_costs(costs, tech_colors):
    costs = costs.groupby(costs.columns, axis=1).sum()
    costs.drop(list(costs.columns[(costs == 0.0).all()]), axis=1, inplace=True)
    new_columns = preferred_order.intersection(costs.columns).append(
        costs.columns.difference(preferred_order)
    )
    costs = costs[new_columns]
    for item in new_columns:
        if item not in tech_colors:
            logger.warning(f"{item} not in config/plotting/tech_colors")
    costs = costs.stack()  # .sort_index()
    to_drop = costs.index.levels[0].symmetric_difference(n.buses.index)
    if len(to_drop) != 0:
        #         logger.info(f"Dropping non-buses {to_drop.tolist()}")
        costs.drop(to_drop, level=0, inplace=True, axis=0, errors="ignore")
    # make sure they are removed from index
    costs.index = pd.MultiIndex.from_tuples(costs.index.values)
    return costs

# Returns the widths for lines and links for plotting, depending on whether
# total or only added grid capacity should be shown.
def get_link_widths(n, attr):
    if attr == 'total':
        line_widths = n.lines.s_nom_opt
        link_widths = n.links.p_nom_opt
        title = "total grid"
    else:
        line_widths = n.lines.s_nom_opt - n.lines.s_nom_min
        link_widths = n.links.p_nom_opt - n.links.p_nom_min
        title = "added grid"

    return line_widths, link_widths, title

# Main function to plot the cost map of the network.
# It visualizes the spatial distribution of system costs and grid expansion.
def plot_cost_map(
        network,
        opts,
        components=["generators", "links", "stores", "storage_units"],
):
    tech_colors = opts["tech_colors"]

    n = network.copy()
    assign_location(n)
    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    # Prepare dataframes to accumulate additional and nominal costs
    costs_add = pd.DataFrame(index=n.buses.index)
    costs_nom = pd.DataFrame(index=n.buses.index)

    # Loop through each component type and calculate costs
    for comp in components:
        df_c = getattr(n, comp)

        if df_c.empty:
            continue

        # Map carrier names to more readable names
        df_c["nice_group"] = df_c.carrier.map(rename_techs)

        # Choose the correct attribute for each component type
        attr = "e_nom_opt" if comp == "stores" else "p_nom_opt"
        attr2 = "e_nom" if comp == "stores" else "p_nom"

        # Calculate additional (added) capital costs
        costs_a = (
            (df_c.capital_cost * (df_c[attr] - df_c[attr2]))
                .groupby([df_c.location, df_c.nice_group])
                .sum()
                .unstack()
                .fillna(0.0)
        )
        costs_add = pd.concat([costs_add, costs_a], axis=1)

        # Calculate total (nominal) capital costs
        costs_n = (
            (df_c.capital_cost * (df_c[attr]))
                .groupby([df_c.location, df_c.nice_group])
                .sum()
                .unstack()
                .fillna(0.0)
        )

        costs_nom = pd.concat([costs_nom, costs_n], axis=1)

    # Process and order the cost data
    costs_add = get_costs(costs_add, tech_colors)
    costs_nom = get_costs(costs_nom, tech_colors)

    # Remove links with zero length (not real transmission lines)
    n.links.drop(n.links.index[n.links.length == 0],
                 inplace=True,
                 )

    # Only show carriers with significant costs (above threshold)
    threshold = 100e6  # 100 mEUR/a
    carriers = pd.concat([costs_add, costs_nom]).groupby(level=1).sum()
    carriers = carriers.where(carriers > threshold).dropna()
    carriers = list(carriers.index)

    # Set thresholds and factors for line widths and colors
    line_lower_threshold = 500.0
    line_upper_threshold = 1e4
    linewidth_factor = opts["map"]["linewidth_factor"]
    total_color = "rosybrown"
    added_color = "darkseagreen"

    # Prepare dataframe for bar plot (system cost breakdown)
    df = pd.DataFrame(index=carriers, columns=["total", "added"])
    df['total'] = costs_nom.groupby(level=1).sum()
    df['added'] = costs_add.groupby(level=1).sum()
    df = df.fillna(0)
    df = df / 1e9  # Convert to bEUR/a
    planning_horizon = int(snakemake.wildcards.planning_horizons)
    df = df / (1 + snakemake.config["costs"]["discountrate"]) ** (planning_horizon - 2020)

    # Order the index for plotting
    new_index = preferred_order.intersection(df.index).append(
        df.index.difference(preferred_order)
    )
    percent = round((df.sum()[1] / df.sum()[0]) * 100)

    # Create the figure and axes for the two maps
    fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={"projection": ccrs.PlateCarree()})
    fig.set_size_inches(opts["map"]["figsize"])

    bus_size_factor = opts["map"]["cost_size_factor"]

    # Plot the total grid (left map)
    line_widths, link_widths, title = get_link_widths(n, "total")

    line_widths = line_widths.clip(line_lower_threshold, line_upper_threshold)
    link_widths = link_widths.clip(line_lower_threshold, line_upper_threshold)

    line_widths = line_widths.replace(line_lower_threshold, 0)
    link_widths = link_widths.replace(line_lower_threshold, 0)

    n.plot(
        bus_sizes=costs_nom / bus_size_factor,
        bus_colors=tech_colors,
        line_colors=total_color,
        link_colors=total_color,
        line_widths=line_widths / linewidth_factor,
        link_widths=link_widths / linewidth_factor,
        ax=ax1,
        color_geomap=True,
        boundaries=opts["map"]["boundaries"]
    )

    # Add legend for line widths (total grid)
    sizes = [10, 5]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / linewidth_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.25, 1.03),
        frameon=False,
        labelspacing=0.8,
        handletextpad=1,
        title=title,
    )

    add_legend_lines(
        ax1, sizes, labels, patch_kw=dict(color=total_color), legend_kw=legend_kw
    )

    # Plot the added grid (right map)
    line_widths, link_widths, title = get_link_widths(n, "added")

    line_widths = line_widths.clip(line_lower_threshold, line_upper_threshold)
    link_widths = link_widths.clip(line_lower_threshold, line_upper_threshold)

    line_widths = line_widths.replace(line_lower_threshold, 0)
    link_widths = link_widths.replace(line_lower_threshold, 0)

    n.plot(
        bus_sizes=costs_add / bus_size_factor,
        bus_colors=tech_colors,
        line_colors=added_color,
        link_colors=added_color,
        line_widths=line_widths / linewidth_factor,
        link_widths=link_widths / linewidth_factor,
        ax=ax2,
        color_geomap=True,
        boundaries=opts["map"]["boundaries"]
    )

    # Add legend for bus sizes (system cost)
    sizes = [20, 10, 5]
    labels = [f"{s} bEUR/a" for s in sizes]
    sizes = [s / bus_size_factor * 1e9 for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.0, 1.03),
        labelspacing=0.8,
        frameon=False,
        handletextpad=0,
        title="system cost",
    )

    add_legend_circles(
        ax1,
        sizes,
        labels,
        srid=n.srid,
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw,
    )

    add_legend_circles(
        ax2,
        sizes,
        labels,
        srid=n.srid,
        patch_kw=dict(facecolor="lightgrey"),
        legend_kw=legend_kw,
    )

    # Add legend for line widths (added grid)
    sizes = [10, 5]
    labels = [f"{s} GW" for s in sizes]
    scale = 1e3 / linewidth_factor
    sizes = [s * scale for s in sizes]

    legend_kw = dict(
        loc="upper left",
        bbox_to_anchor=(0.25, 1.03),
        frameon=False,
        labelspacing=0.8,
        handletextpad=1,
        title=title,
    )

    add_legend_lines(
        ax2, sizes, labels, patch_kw=dict(color=added_color), legend_kw=legend_kw
    )

    # Add legend for technology colors
    legend_kw = dict(
        bbox_to_anchor=(1.42, 1.04),
        frameon=False,
    )

    colors = [tech_colors[c] for c in carriers] + [total_color, added_color]
    labels = carriers + ["HVDC or HVAC link", "HVDC or HVAC link"]

    add_legend_patches(
        ax2,
        colors,
        labels,
        legend_kw=legend_kw,
    )

    # Add a bar plot showing the system cost breakdown by technology
    ax3 = fig.add_axes([-0.09, 0.28, 0.09, 0.45])

    df.loc[new_index, df.columns].T.plot(
        kind="bar",
        ax=ax3,
        stacked=True,
        color=[tech_colors[i] for i in new_index],
    )
    ax3.legend().remove()
    ax3.set_ylabel("annualized system cost bEUR/a")
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation='horizontal')
    ax3.grid(axis="y")
    ax3.set_ylim([0, opts["costs_max"]])
    ax3.text(0.85, (df.sum()[1] + 15), str(percent) + "%", color='black')

    fig.tight_layout()

    # Save the figure to the output path specified by snakemake
    fig.savefig(snakemake.output.cost_map, transparent=True, bbox_inches="tight")

# Main script entry point
if __name__ == "__main__":
    # If not running from snakemake, mock the snakemake object for testing
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_network',
                                   opts='ll',
                                   topology ='current+Neighbor',
                                   pathway ='exponential175',
                                   planning_horizons="2020")
    # Set logging level from config
    logging.basicConfig(level=snakemake.config["logging"]["level"])

    set_plot_style()

    config = snakemake.config

    # Load the network from file
    n = pypsa.Network(snakemake.input.network)

    # Call the main plotting function
    plot_cost_map(
        n,
        opts=config["plotting"],
        components=["generators", "links", "stores", "storage_units"],
    )