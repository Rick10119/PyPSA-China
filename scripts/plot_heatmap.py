"""
This script generates heatmap visualizations for energy storage technologies and aluminum smelters in the PyPSA-China model.
It creates two types of plots:
1. Heatmaps showing the temporal patterns of storage operation for H2, battery, water storage, and aluminum smelters
2. A time series plot showing the water tank storage level over time

The heatmaps show how storage technologies and aluminum smelters are used throughout the year, with hours on the y-axis
and days on the x-axis. The color intensity represents the normalized power output/input.

This script now supports both national and province-specific heatmaps, saving them to different directories.
"""

from _helpers import configure_logging
import seaborn as sns
import pandas as pd
import pypsa
import matplotlib.pyplot as plt
import os

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

def creat_df(n, tech, province_filter=None):
    """
    Creates a DataFrame for heatmap visualization of a specific storage technology.
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object containing the simulation results
    tech : str
        The storage technology to analyze ('H2', 'battery', or 'water')
    province_filter : str, optional
        If provided, filter stores to only include those in the specified province
    
    Returns:
    --------
    tuple
        (summary DataFrame for heatmap, base power value in MW)
    """
    # Get stores for the specific technology
    stores = n.stores_t.p.filter(like=tech)
    
    # Apply province filter if specified
    if province_filter:
        # Filter stores to only include those in the specified province
        province_stores = stores.columns[stores.columns.str.contains(province_filter, case=False)]
        if len(province_stores) > 0:
            stores = stores[province_stores]
        else:
            print(f"Warning: No {tech} stores found in province {province_filter}")
            return pd.DataFrame(), 0
    
    # Calculate the maximum power as the base value for normalization
    base = abs(stores.sum(axis=1)).max()
    if base == 0:
        print(f"Warning: {tech} power is zero for {'province ' + province_filter if province_filter else 'national'}")
        return pd.DataFrame(), 0
    
    # Normalize the power values by the base value
    df = stores.sum(axis=1) / base
    df = df.to_frame()
    df.reset_index(inplace=True)
    renames = {0: 'p_store'}
    df.rename(columns=renames, inplace=True)
    
    # 数据已经是中国本地时间，直接使用时间戳，不需要时区转换
    date = n.stores_t.p.filter(like='water').index
    df['Hour'] = date.hour
    df['Day'] = date.strftime('%m-%d')
    
    # Create pivot table for heatmap visualization
    summary = pd.pivot_table(data=df, index='Hour', columns='Day', values='p_store')
    summary = summary.fillna(0).infer_objects(copy=False)
    return summary, base

def creat_aluminum_df(n, province_filter=None):
    """
    Creates a DataFrame for heatmap visualization of aluminum smelter operation.
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object containing the simulation results
    province_filter : str, optional
        If provided, filter aluminum smelters to only include those in the specified province
    
    Returns:
    --------
    tuple
        (summary DataFrame for heatmap, base power value in MW)
    """
    # Get aluminum smelter links (input power from electricity bus)
    aluminum_links = n.links_t.p0.filter(like='aluminum smelter')
    
    if aluminum_links.empty:
        # If no aluminum smelters found, return empty DataFrame
        print("Warning: No aluminum smelter links found in the network")
        return pd.DataFrame(), 0
    
    # Apply province filter if specified
    if province_filter:
        # Filter aluminum smelters to only include those in the specified province
        province_smelters = aluminum_links.columns[aluminum_links.columns.str.contains(province_filter, case=False)]
        if len(province_smelters) > 0:
            aluminum_links = aluminum_links[province_smelters]
        else:
            print(f"Warning: No aluminum smelters found in province {province_filter}")
            return pd.DataFrame(), 0
    
    # Calculate the maximum power as the base value for normalization
    base = abs(aluminum_links.sum(axis=1)).max()
    
    if base == 0:
        print(f"Warning: Aluminum smelter power is zero for {'province ' + province_filter if province_filter else 'national'}")
        return pd.DataFrame(), 0
    
    # Normalize the power values by the base value
    df = aluminum_links.sum(axis=1) / base
    df = df.to_frame()
    df.reset_index(inplace=True)
    renames = {0: 'p_smelter'}
    df.rename(columns=renames, inplace=True)
    
    # 数据已经是中国本地时间，直接使用时间戳，不需要时区转换
    date = aluminum_links.index
    df['Hour'] = date.hour
    df['Day'] = date.strftime('%m-%d')
    
    # Create pivot table for heatmap visualization
    summary = pd.pivot_table(data=df, index='Hour', columns='Day', values='p_smelter')
    summary = summary.fillna(0).infer_objects(copy=False)
    return summary, base

def plot_heatmap(n, config, output_dir, province_filter=None):
    """
    Generates heatmap plots for all storage technologies and aluminum smelters.
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object containing the simulation results
    config : dict
        Configuration dictionary containing plotting parameters
    output_dir : str
        Directory to save the heatmap plots
    province_filter : str, optional
        If provided, filter data to only include the specified province
    """
    techs = ["H2", "battery", "water"]
    freq = config["freq"]
    planning_horizon = snakemake.wildcards.planning_horizons
    
    # Create province-specific subdirectory if filtering by province
    if province_filter:
        province_dir = os.path.join(output_dir, f"province_{province_filter}")
        os.makedirs(province_dir, exist_ok=True)
        plot_title_suffix = f" in {province_filter}"
    else:
        province_dir = output_dir
        plot_title_suffix = " (National)"
    
    # Plot storage technologies
    for tech in techs:
        fig, ax = plt.subplots(figsize=map_figsize)
        df, base = creat_df(n, tech, province_filter)
        if not df.empty and base > 0:
            base = str(int(base / 1e3))  # Convert to GW for display
            
            # Create heatmap with coolwarm colormap and normalized values
            sns.heatmap(df, ax=ax, cmap='coolwarm', cbar_kws={'label': 'pu'}, vmin=-1.0, vmax=1.0)
            ax.set_title(tech + ' heatmap with ' + freq + ' resolution in ' + planning_horizon + 
                        plot_title_suffix + ' P_base = ' + base + ' GW')
            
            # Save to province-specific directory
            output_path = os.path.join(province_dir, f"{tech}-{snakemake.wildcards.opts}-{snakemake.wildcards.topology}-{snakemake.wildcards.pathway}-{planning_horizon}.png")
            fig.savefig(output_path, dpi=150, bbox_inches='tight')
            print(f"Saved {tech} heatmap to: {output_path}")
        plt.close()
    
    # Plot aluminum smelter heatmap only when aluminum is added and data exists
    if config.get("add_aluminum", False):
        # First check if aluminum data exists before creating any plots
        df, base = creat_aluminum_df(n, province_filter)
        
        # Only create and save the plot if we have valid aluminum data
        if not df.empty and base > 0:
            # Create figure with single heatmap
            fig, ax = plt.subplots(figsize=map_figsize)
            
            base = str(int(base / 1e3))  # Convert to GW for display
            
            # Create heatmap
            sns.heatmap(df, ax=ax, cmap='coolwarm', cbar_kws={'label': 'pu'}, vmin=0.0, vmax=1.0)
            ax.set_title('Aluminum Smelter heatmap with ' + freq + ' resolution in ' + planning_horizon + 
                        plot_title_suffix + ' P_base = ' + base + ' GW')
            
            # Save to province-specific directory
            output_path = os.path.join(province_dir, f"aluminum-{snakemake.wildcards.opts}-{snakemake.wildcards.topology}-{snakemake.wildcards.pathway}-{planning_horizon}.png")
            fig.savefig(output_path, dpi=150, bbox_inches='tight')
            print(f"Saved aluminum heatmap to: {output_path}")
            plt.close()
        else:
            print(f"Skipping aluminum heatmap for {'province ' + province_filter if province_filter else 'national'}: No data found or power is zero")

def plot_water_store(n, output_dir, province_filter=None):
    """
    Generates a time series plot of water tank storage levels.
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object containing the simulation results
    output_dir : str
        Directory to save the water storage plot
    province_filter : str, optional
        If provided, filter data to only include the specified province
    """
    planning_horizon = snakemake.wildcards.planning_horizons
    
    # Create province-specific subdirectory if filtering by province
    if province_filter:
        province_dir = os.path.join(output_dir, f"province_{province_filter}")
        os.makedirs(province_dir, exist_ok=True)
        plot_title_suffix = f" in {province_filter}"
    else:
        province_dir = output_dir
        plot_title_suffix = " (National)"
    
    fig, ax = plt.subplots(figsize=map_figsize)
    
    # Get water storage data
    water_stores = n.stores_t.e.filter(like='water')
    water_capacity = n.stores.e_nom_opt.filter(like='water')
    
    # Apply province filter if specified
    if province_filter:
        province_stores = water_stores.columns[water_stores.columns.str.contains(province_filter, case=False)]
        if len(province_stores) > 0:
            water_stores = water_stores[province_stores]
            water_capacity = water_capacity[water_capacity.index.isin(province_stores)]
        else:
            print(f"Warning: No water stores found in province {province_filter}")
            plt.close()
            return
    
    # Calculate and plot the normalized water storage level over time
    if not water_stores.empty and not water_capacity.empty:
        (water_stores.sum(axis=1) / water_capacity.sum()).plot(ax=ax)
        
        # Set y-axis limits and ticks for better visualization
        ax.set_ylim(0, 1.0)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_title("Water tank storage in " + planning_horizon + plot_title_suffix)
        
        # Save to province-specific directory
        output_path = os.path.join(province_dir, f"water_store-{snakemake.wildcards.opts}-{snakemake.wildcards.topology}-{snakemake.wildcards.pathway}-{planning_horizon}.png")
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Saved water storage plot to: {output_path}")
    else:
        print(f"No water storage data found for {'province ' + province_filter if province_filter else 'national'}")
    
    plt.close()

def save_to_snakemake_outputs(n, config):
    """
    Save plots to the exact paths expected by Snakefile outputs.
    This ensures compatibility with the existing workflow.
    """
    techs = ["H2", "battery", "water"]
    freq = config["freq"]
    planning_horizon = snakemake.wildcards.planning_horizons
    
    # Plot storage technologies to Snakefile output paths
    for tech in techs:
        fig, ax = plt.subplots(figsize=map_figsize)
        df, base = creat_df(n, tech)
        if not df.empty and base > 0:
            base = str(int(base / 1e3))  # Convert to GW for display
            
            # Create heatmap with coolwarm colormap and normalized values
            sns.heatmap(df, ax=ax, cmap='coolwarm', cbar_kws={'label': 'pu'}, vmin=-1.0, vmax=1.0)
            ax.set_title(tech + ' heatmap with ' + freq + ' resolution in ' + planning_horizon + ' P_base = ' + base + ' GW')
            
            # Save to Snakefile output path
            if hasattr(snakemake.output, tech):
                fig.savefig(snakemake.output[tech], dpi=150, bbox_inches='tight')
                print(f"Saved {tech} heatmap to Snakefile output: {snakemake.output[tech]}")
            else:
                print(f"Warning: No Snakefile output defined for {tech}")
        plt.close()
    
    # Plot aluminum smelter heatmap to Snakefile output path
    if config.get("add_aluminum", False):
        fig, ax = plt.subplots(figsize=map_figsize)
        df, base = creat_aluminum_df(n)
        if not df.empty and base > 0:
            base = str(int(base / 1e3))  # Convert to GW for display
            
            # Create heatmap
            sns.heatmap(df, ax=ax, cmap='coolwarm', cbar_kws={'label': 'pu'}, vmin=0.0, vmax=1.0)
            ax.set_title('Aluminum Smelter heatmap with ' + freq + ' resolution in ' + planning_horizon + ' P_base = ' + base + ' GW')
            
            # Save to Snakefile output path
            if hasattr(snakemake.output, "aluminum"):
                fig.savefig(snakemake.output["aluminum"], dpi=150, bbox_inches='tight')
                print(f"Saved aluminum heatmap to Snakefile output: {snakemake.output['aluminum']}")
            else:
                print("Warning: No Snakefile output defined for aluminum")
        else:
            print("Warning: No aluminum data found or power is zero")
        plt.close()
    
    # Plot water tank storage to Snakefile output path
    fig, ax = plt.subplots(figsize=map_figsize)
    water_stores = n.stores_t.e.filter(like='water')
    water_capacity = n.stores.e_nom_opt.filter(like='water')
    
    if not water_stores.empty and not water_capacity.empty:
        (water_stores.sum(axis=1) / water_capacity.sum()).plot(ax=ax)
        ax.set_ylim(0, 1.0)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_title("Water tank storage in " + planning_horizon)
        
        # Save to Snakefile output path
        if hasattr(snakemake.output, "water"):
            fig.savefig(snakemake.output["water"], dpi=150, bbox_inches='tight')
            print(f"Saved water tank plot to Snakefile output: {snakemake.output['water']}")
        else:
            print("Warning: No Snakefile output defined for water")
    plt.close()
    
    # Plot water store to Snakefile output path
    fig, ax = plt.subplots(figsize=map_figsize)
    if not water_stores.empty and not water_capacity.empty:
        (water_stores.sum(axis=1) / water_capacity.sum()).plot(ax=ax)
        ax.set_ylim(0, 1.0)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_title("Water store in " + planning_horizon)
        
        # Save to Snakefile output path
        if hasattr(snakemake.output, "water_store"):
            fig.savefig(snakemake.output["water_store"], dpi=150, bbox_inches='tight')
            print(f"Saved water store plot to Snakefile output: {snakemake.output['water_store']}")
        else:
            print("Warning: No Snakefile output defined for water_store")
    plt.close()

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
    
    # First, save to Snakefile output paths to ensure workflow compatibility
    print("Generating plots for Snakefile outputs...")
    save_to_snakemake_outputs(n, config)
    
    # Get base output directory from snakemake output for additional plots
    base_output_dir = os.path.dirname(snakemake.output["H2"])
    
    # Plot national heatmaps to additional directory
    print("\nGenerating national heatmaps to additional directory...")
    plot_heatmap(n, config, base_output_dir)
    plot_water_store(n, base_output_dir)
    
    # Plot province-specific heatmaps if single node mode is enabled
    if config.get("using_single_node", False) and config.get("single_node_province"):
        province = config["single_node_province"]
        print(f"\nGenerating province-specific heatmaps for {province}...")
        plot_heatmap(n, config, base_output_dir, province)
        plot_water_store(n, base_output_dir, province)