"""
This script generates heatmap visualizations for energy storage technologies in the PyPSA-China model.
It creates two types of plots:
1. Heatmaps showing the temporal patterns of storage operation for H2, battery, and water storage
2. A time series plot showing the water tank storage level over time

The heatmaps show how storage technologies are used throughout the year, with hours on the y-axis
and days on the x-axis. The color intensity represents the normalized power output/input.
"""

from _helpers import configure_logging
import seaborn as sns
import pandas as pd
import pypsa
import matplotlib.pyplot as plt

def set_plot_style():
    """
    Sets up the plotting style for all matplotlib plots in this script.
    Uses a combination of classic and seaborn styles with custom modifications
    for better visualization quality.
    """
    plt.style.use(['classic', 'seaborn-v0_8-whitegrid',
                   {'axes.grid': False, 'grid.linestyle': '--', 'grid.color': u'0.6',
                    'hatch.color': 'white',
                    'patch.linewidth': 0.5,
                    'font.size': 12,
                    'legend.fontsize': 'medium',
                    'lines.linewidth': 1.5,
                    'pdf.fonttype': 42,
                    }])

def creat_df(n, tech):
    """
    Creates a DataFrame for heatmap visualization of a specific storage technology.
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object containing the simulation results
    tech : str
        The storage technology to analyze ('H2', 'battery', or 'water')
    
    Returns:
    --------
    tuple
        (summary DataFrame for heatmap, base power value in MW)
    """
    # Calculate the maximum power as the base value for normalization
    base = abs(n.stores_t.p.filter(like=tech).sum(axis=1)).max()
    # Normalize the power values by the base value
    df = (n.stores_t.p.filter(like=tech).sum(axis=1))/base
    df = df.to_frame()
    df.reset_index(inplace=True)
    renames = {0: 'p_store'}
    df.rename(columns=renames, inplace=True)
    
    # Convert timestamps to local time (Asia/Shanghai) and extract hour and day
    date = n.stores_t.p.filter(like='water').index
    date = date.tz_localize('utc')
    date = date.tz_convert("Asia/Shanghai")
    df['Hour'] = date.hour
    df['Day'] = date.strftime('%m-%d')
    
    # Create pivot table for heatmap visualization
    summary = pd.pivot_table(data=df,index='Hour',columns='Day',values='p_store')
    summary = summary.fillna(0)
    return summary, base

def plot_heatmap(n, config):
    """
    Generates heatmap plots for all storage technologies.
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object containing the simulation results
    config : dict
        Configuration dictionary containing plotting parameters
    """
    techs = ["H2", "battery", "water"]
    freq = config["freq"]
    planning_horizon = snakemake.wildcards.planning_horizons
    
    for tech in techs:
        fig, ax = plt.subplots(figsize=map_figsize)
        df, base = creat_df(n, tech)
        base = str(int(base / 1e3))  # Convert to GW for display
        
        # Create heatmap with coolwarm colormap and normalized values
        sns.heatmap(df, ax=ax, cmap='coolwarm', cbar_kws={'label': 'pu'},vmin=-1.0, vmax=1.0)
        ax.set_title(tech + ' heatmap with ' + freq + ' resolution in ' + planning_horizon + ' P_base = ' + base + ' GW')
        fig.savefig(snakemake.output[tech], dpi=150, bbox_inches='tight')

def plot_water_store(n):
    """
    Generates a time series plot of water tank storage levels.
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object containing the simulation results
    """
    planning_horizon = snakemake.wildcards.planning_horizons
    fig, ax = plt.subplots(figsize=map_figsize)
    
    # Calculate and plot the normalized water storage level over time
    (n.stores_t.e.filter(like='water').sum(axis=1) / n.stores.e_nom_opt.filter(like='water').sum()).plot(ax=ax)
    
    # Set y-axis limits and ticks for better visualization
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_title(" water tank storage in " + planning_horizon)
    fig.savefig(snakemake.output["water_store"], dpi=150, bbox_inches='tight')


if __name__ == "__main__":
    # Set up mock snakemake for testing if not running in snakemake
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_heatmap',
                                   opts='ll',
                                   topology ='current+Neighbor',
                                   pathway ='exponential175',
                                   planning_horizons="2020")
    configure_logging(snakemake)

    # Initialize plotting style and load configuration
    set_plot_style()
    config = snakemake.config

    # Get plotting parameters from config
    map_figsize = config["plotting"]['map']['figsize']
    map_boundaries = config["plotting"]['map']['boundaries']

    # Load the network and generate plots
    n = pypsa.Network(snakemake.input.network)
    plot_heatmap(n, config)
    plot_water_store(n)