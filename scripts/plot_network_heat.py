"""
Network Heat Map Plotting Script for PyPSA-China

This script creates visualizations of power system networks including:
- Geographic maps showing transmission lines, generators, and storage
- Energy mix pie charts
- Cost breakdown bar charts

The script is designed to work with PyPSA (Python for Power System Analysis) networks
and generates publication-ready figures for energy system analysis.
"""

import logging
from _helpers import (load_network_for_plots, aggregate_p, aggregate_costs, configure_logging)
from functions import pro_names

import pandas as pd
import numpy as np

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle, Ellipse
from matplotlib.legend_handler import HandlerPatch
to_rgba = mpl.colors.colorConverter.to_rgba

logger = logging.getLogger(__name__)

# Dictionary mapping component types to their optimization variable names
opt_name = {"Store": "e", "Line" : "s", "Transformer" : "s"}


def make_handler_map_to_scale_circles_as_in(ax, dont_resize_actively=False):
    """
    Creates a custom legend handler for circles that scales properly with the plot.
    
    This function ensures that circles in the legend maintain their relative sizes
    when the plot is resized or zoomed, making the legend more informative.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes object to attach the handler to
    dont_resize_actively : bool
        If True, the circles won't resize when the plot is resized
        
    Returns:
    --------
    dict : Handler map for matplotlib legend
    """
    fig = ax.get_figure()

    def axes2pt():
        """Convert axes coordinates to points for proper scaling"""
        return np.diff(ax.transData.transform([(0, 0), (1, 1)]), axis=0)[0] * (72. / fig.dpi)

    ellipses = []
    if not dont_resize_actively:
        def update_width_height(event):
            """Update circle sizes when plot is resized"""
            dist = axes2pt()
            for e, radius in ellipses: e.width, e.height = 2. * radius * dist

        # Connect resize events to update function
        fig.canvas.mpl_connect('resize_event', update_width_height)
        ax.callbacks.connect('xlim_changed', update_width_height)
        ax.callbacks.connect('ylim_changed', update_width_height)

    def legend_circle_handler(legend, orig_handle, xdescent, ydescent,
                              width, height, fontsize):
        """Create properly scaled circles for legend"""
        w, h = 2. * orig_handle.get_radius() * axes2pt()
        e = Ellipse(xy=(0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent), width=w, height=w)
        ellipses.append((e, orig_handle.get_radius()))
        return e

    return {Circle: HandlerPatch(patch_func=legend_circle_handler)}


def make_legend_circles_for(sizes, scale=1.0, **kw):
    """
    Create circle objects for legend with specified sizes.
    
    Parameters:
    -----------
    sizes : list
        List of sizes to create circles for
    scale : float
        Scaling factor for the circle sizes
    **kw : dict
        Additional keyword arguments for Circle objects
        
    Returns:
    --------
    list : List of Circle objects
    """
    return [Circle((0, 0), radius=(s / scale) ** 0.5, **kw) for s in sizes]


def set_plot_style():
    """
    Configure the matplotlib plotting style for consistent, publication-ready figures.
    
    Sets up a clean, professional appearance with appropriate line widths,
    font sizes, and grid styling.
    """
    plt.style.use(['default',
                   {'axes.grid': False, 'grid.linestyle': '--', 'grid.color': u'0.6',
                    'hatch.color': 'white',
                    'patch.linewidth': 0.5,
                    'font.size': 12,
                    'legend.fontsize': 'medium',
                    'lines.linewidth': 1.5,
                    'pdf.fonttype': 42,
                    'figure.facecolor': 'white',
                    'axes.facecolor': 'white',
                    'savefig.facecolor': 'white',
                    'savefig.edgecolor': 'white',
                    }])


def plot_opt_map(n, opts, ax=None, attribute='p_nom'):
    """
    Create a geographic map visualization of the power system network.
    
    This function plots:
    - Transmission lines (existing and expansion)
    - Power plants (generators) with size proportional to capacity
    - Storage facilities
    - Color-coded by technology type
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object containing the power system data
    opts : dict
        Configuration options for plotting (colors, scaling factors, etc.)
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, uses current axes.
    attribute : str
        Network attribute to visualize (currently only 'p_nom' supported)
        
    Returns:
    --------
    matplotlib.figure.Figure : The created figure object
    """
    if ax is None:
        ax = plt.gca()

    ## DATA PREPARATION
    # Define colors for existing vs expansion infrastructure
    line_colors = {'cur': "purple",
                   'exp': mpl.colors.rgb2hex(to_rgba("red", 0.7), True)}
    tech_colors = opts['tech_colors']

    if attribute == 'p_nom':
        # Aggregate generator capacities by bus and carrier (technology)
        # Exclude solar thermal and hydro inflow for clarity
        bus_sizes = pd.concat((n.generators.query('carrier != "solar thermal"' and 'carrier != "hydro_inflow"').groupby(['bus', 'carrier']).p_nom_opt.sum(),
                               n.links.query('carrier == ["gas-AC","coal-AC","stations-AC"]').groupby(['bus1', 'carrier']).p_nom_opt.sum()))
        bus_sizes.index.names = ['bus', 'carrier']
        bus_sizes = bus_sizes.groupby(['bus','carrier']).sum()
        
        # Filter out carriers that don't have colors defined
        available_colors = set(tech_colors.keys())
        available_carriers = bus_sizes.index.get_level_values('carrier').unique()
        missing_carriers = set(available_carriers) - available_colors
        
        if missing_carriers:
            logger.warning(f"Missing colors for carriers: {missing_carriers}. Filtering them out.")
            # Filter out carriers without colors
            bus_sizes = bus_sizes[bus_sizes.index.get_level_values('carrier').isin(available_colors)]
        
        # Check if we have any carriers left to plot
        if bus_sizes.empty:
            logger.warning("No carriers with defined colors found. Creating empty bus_sizes.")
            bus_sizes = pd.Series(dtype=float)
        else:
            logger.info(f"Plotting carriers: {bus_sizes.index.get_level_values('carrier').unique()}")
        
        # Get transmission line capacities (existing and expansion)
        line_widths_exp = n.lines.s_nom_opt  # Optimal (including expansion)
        line_widths_cur = n.lines.s_nom_min  # Existing minimum capacity
        
        # Get link capacities (transformers, converters, etc.)
        ac_ac_links = n.links.query('carrier == ["AC-AC"]').p_nom_opt
        other_links = n.links.query('carrier != ["AC-AC"]').p_min_pu
        link_widths_exp = pd.concat([ac_ac_links, other_links]) - n.links.p_nom_min
        link_widths_cur = n.links.p_nom_min
    else:
        raise 'plotting of {} has not been implemented yet'.format(attribute)

    # Create alpha masks for existing infrastructure visibility
    line_colors_with_alpha = \
        ((line_widths_cur / n.lines.s_nom > 1e-3)
         .map({True: line_colors['cur'], False: to_rgba(line_colors['cur'], 0.)}))
    link_colors_with_alpha = \
        ((link_widths_cur / n.links.p_nom > 1e-3)
         .map({True: line_colors['cur'], False: to_rgba(line_colors['cur'], 0.)}))

    ## FORMATTING
    # Get scaling factors from configuration
    linewidth_factor = opts['map']['linewidth_factor']
    bus_size_factor = opts['map']['bus_size_factor']

    ## PLOTTING
    # Plot expansion infrastructure (red lines, colored generators)
    n.plot(line_widths=line_widths_exp / linewidth_factor,
           link_widths=link_widths_exp / linewidth_factor,
           line_colors=line_colors['exp'],
           link_colors=line_colors['exp'],
           bus_sizes=bus_sizes / bus_size_factor,
           bus_colors=tech_colors,
           boundaries=map_boundaries,
           color_geomap=True, geomap=True,
           ax=ax)
    
    # Plot existing infrastructure (purple lines, transparent generators)
    n.plot(line_widths=line_widths_cur / linewidth_factor,
           link_widths=link_widths_cur / linewidth_factor,
           line_colors=line_colors_with_alpha,
           link_colors=link_colors_with_alpha,
           bus_sizes=0,
           boundaries=map_boundaries,
           color_geomap=True, geomap=True,
           ax=ax)
    
    # Set plot properties
    ax.set_aspect('equal')
    ax.axis('off')

    # Rasterize basemap for better performance
    # TODO : Check if this also works with cartopy
    for c in ax.collections[:2]: c.set_rasterized(True)

    ## LEGEND CREATION
    handles = []
    labels = []

    # Transmission line legend (expansion)
    for s in (50, 10):
        handles.append(plt.Line2D([0], [0], color=line_colors['exp'],
                                  linewidth=s * 1e3 / linewidth_factor))
        labels.append("{} GW".format(s))
    l1_1 = ax.legend(handles, labels,
                     loc="upper left", bbox_to_anchor=(0.24, 1.01),
                     frameon=False,
                     labelspacing=0.8, handletextpad=1.5,
                     title='Transmission Exist./Exp.             ')
    ax.add_artist(l1_1)

    # Transmission line legend (existing)
    handles = []
    labels = []
    for s in (50, 10):
        handles.append(plt.Line2D([0], [0], color=line_colors['cur'],
                                  linewidth=s * 1e3 / linewidth_factor))
        labels.append("/")
    l1_2 = ax.legend(handles, labels,
                     loc="upper left", bbox_to_anchor=(0.26, 1.01),
                     frameon=False,
                     labelspacing=0.8, handletextpad=0.5,
                     title=' ')
    ax.add_artist(l1_2)

    # Generation capacity legend (circles)
    handles = make_legend_circles_for([10e4, 5e4, 1e4], scale=bus_size_factor, facecolor="w")
    labels = ["{} GW".format(s) for s in (100, 50, 10)]
    l2 = ax.legend(handles, labels,
                   loc="upper left", bbox_to_anchor=(0.01, 1.01),
                   frameon=False, labelspacing=1.0,
                   title='Generation',
                   handler_map=make_handler_map_to_scale_circles_as_in(ax))
    ax.add_artist(l2)

    # Technology color legend
    techs = (bus_sizes.index.levels[1]).intersection(
        pd.Index(opts['vre_techs'] + opts['conv_techs'] + opts['storage_techs']))
    # Filter to only include technologies that have defined colors
    techs_with_colors = [t for t in techs if t in tech_colors]
    handles = []
    labels = []
    for t in techs_with_colors:
        handles.append(plt.Line2D([0], [0], color=tech_colors[t], marker='o', markersize=8, linewidth=0))
        labels.append(opts['nice_names'].get(t, t))
    l3 = ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.),  # bbox_to_anchor=(0.72, -0.05),
                   handletextpad=0., columnspacing=0.5, ncol=4, title='Technology')

    return fig


def plot_total_energy_pie(n, opts, ax=None):
    """
    Create a pie chart showing the energy mix by technology.
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object
    opts : dict
        Configuration options for plotting
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, uses current axes.
    """
    if ax is None: ax = plt.gca()

    ax.set_title('Energy per technology', fontdict=dict(fontsize="medium"))

    # Aggregate power by carrier and filter out small contributions
    e_primary = aggregate_p(n).drop('load', errors='ignore').loc[lambda s: s > 1]
    e_primary = e_primary.groupby('carrier').sum()
    
    # Filter out technologies that don't have defined colors
    tech_colors = opts['tech_colors']
    available_colors = set(tech_colors.keys())
    e_primary_filtered = e_primary[e_primary.index.isin(available_colors)]
    
    if e_primary_filtered.empty:
        logger.warning("No technologies with defined colors found for energy pie chart.")
        return
    
    # Create pie chart
    patches, texts, autotexts = ax.pie(e_primary_filtered,
                                       startangle=90,
                                       labels=e_primary_filtered.rename(opts['nice_names'].get('energy', {})).index,
                                       autopct='%.0f%%',
                                       shadow=False,
                                       colors=[tech_colors[tech] for tech in e_primary_filtered.index])
    
    # Remove labels for very small contributions (< 4% of total)
    for t1, t2, i in zip(texts, autotexts, e_primary_filtered.index):
        if e_primary_filtered.at[i] < 0.04 * e_primary_filtered.sum():
            t1.remove()
            t2.remove()


def plot_total_cost_bar(n, opts, ax=None):
    """
    Create a stacked bar chart showing system costs by technology.
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object
    opts : dict
        Configuration options for plotting
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, uses current axes.
    """
    if ax is None: ax = plt.gca()

    # Calculate total load for cost normalization
    total_load = (n.snapshot_weightings.generators * n.loads_t.p.sum(axis=1)).sum()
    tech_colors = opts['tech_colors']

    def split_costs(n):
        """
        Split costs into different categories: total, existing capital, new capital, and marginal.
        
        Returns:
        --------
        tuple : (total_costs, existing_capital, new_capital, marginal_costs)
        """
        costs = aggregate_costs(n).reset_index(level=0, drop=True)
        costs.index.rename(['cost','carrier'],inplace=True)
        costs = costs.groupby(['cost','carrier']).sum()
        costs_ex = aggregate_costs(n, existing_only=True).reset_index(level=0, drop=True)
        costs_ex.index.rename(['cost','carrier'],inplace=True)
        costs_ex = costs_ex.groupby(['cost','carrier']).sum()
        return (costs['capital'].add(costs['marginal'], fill_value=0.),
                costs_ex['capital'], costs['capital'] - costs_ex['capital'], costs['marginal'])

    costs, costs_cap_ex, costs_cap_new, costs_marg = split_costs(n)

    # Filter costs above threshold for clarity
    costs_graph = pd.DataFrame(dict(a=costs[costs > opts['costs_threshold']])).dropna()
    
    # Filter out technologies that don't have defined colors
    available_colors = set(tech_colors.keys())
    costs_graph_filtered = costs_graph[costs_graph.index.isin(available_colors)]
    
    if costs_graph_filtered.empty:
        logger.warning("No technologies with defined colors found for cost bar chart.")
        return

    bottom = np.array([0., 0.])
    texts = []

    # Create stacked bar chart
    for i, ind in enumerate(costs_graph_filtered.index):
        data = np.asarray(costs_graph_filtered.loc[ind]) / total_load
        ax.bar([0.5], data, bottom=bottom, color=tech_colors[ind],
               width=0.7, zorder=-1)
        bottom_sub = bottom
        bottom = bottom + data

        # Add sub-components for conventional technologies and transmission
        if ind in opts['conv_techs'] + ['AC line']:
            for c in [costs_cap_ex, costs_marg]:
                if ind in c:
                    data_sub = np.asarray([c.loc[ind]]) / total_load
                    ax.bar([0.5], data_sub, linewidth=0,
                           bottom=bottom_sub, color=tech_colors[ind],
                           width=0.7, zorder=-1, alpha=0.8)
                    bottom_sub += data_sub

        # Add labels for significant cost components
        if abs(data[-1]) < 5:
            continue

        text = ax.text(1.1, (bottom - 0.5 * data)[-1] - 3, opts['nice_names'].get(ind, ind))
        texts.append(text)

    # Format the plot
    ax.set_ylabel("Average system cost [Eur/MWh]")
    ax.set_ylim([0, opts.get('costs_avg', 80)])
    ax.set_xlim([0, 1])
    ax.set_xticklabels([])
    ax.grid(True, axis="y", color='k', linestyle='dotted')



if __name__ == "__main__":
    # Main execution block - creates the complete network visualization
    
    # Set up snakemake environment (for workflow integration)
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_network', opts='ll', planning_horizons=2020)
    configure_logging(snakemake)

    # Set up plotting style
    set_plot_style()

    # Extract configuration and parameters
    config, wildcards = snakemake.config, snakemake.wildcards

    map_figsize = config["plotting"]['map']['figsize']
    map_boundaries = config["plotting"]['map']['boundaries']

    cost_year = snakemake.wildcards.planning_horizons

    # Load the network data
    n = load_network_for_plots(snakemake.input.network, snakemake.input.tech_costs, config, cost_year)

    scenario_opts = wildcards.opts.split('-')

    # Create the main map figure
    fig, ax = plt.subplots(figsize=map_figsize, subplot_kw={"projection": ccrs.PlateCarree()})
    plot_opt_map(n, config["plotting"], ax=ax)

    # Save the map-only version
    fig.savefig(snakemake.output.only_map, dpi=150, bbox_inches='tight')

    # Add energy mix pie chart
    ax1 = fig.add_axes([-0.115, 0.625, 0.2, 0.2])
    plot_total_energy_pie(n, config["plotting"], ax=ax1)

    # Add cost breakdown bar chart
    ax2 = fig.add_axes([-0.075, 0.1, 0.1, 0.45])
    plot_total_cost_bar(n, config["plotting"], ax=ax2)

    # Save the complete figure with all components
    fig.savefig(snakemake.output.ext, transparent=True, bbox_inches='tight')
