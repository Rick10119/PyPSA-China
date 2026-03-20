# SPDX-FileCopyrightText: : 2025 Ruike Lyu, rl8728@princeton.edu
"""
Heatmap comparison script for MMMU and NMMU scenarios.

This script is a refactored version of plot_heatmap.py to compare
two different configuration scenarios for storage technologies and
aluminum smelter operation.

Main features:
1. Read parameters from MMMU and NMMU config files
2. Load corresponding network results
3. Generate side-by-side comparison heatmaps
4. Support visualization of H2, battery, water storage and aluminum smelters

Scenario differences:
- MMMU: iterative_optimization: true, smelter_flexibility: mid, employment_transfer: unfavorable (U)
- NMMU: iterative_optimization: false, smelter_flexibility: non_constrained, employment_transfer: unfavorable (U)
"""

import yaml
import seaborn as sns
import pandas as pd
import pypsa
import matplotlib.pyplot as plt
import os
import argparse
from pathlib import Path

def set_plot_style():
    """
    Set global plotting style.
    """
    # Use Helvetica-like fonts
    plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'sans-serif']
    plt.rcParams['axes.unicode_minus'] = False
    
    plt.style.use(['classic', 'seaborn-v0_8-whitegrid',
                   {'axes.grid': False, 'grid.linestyle': '--', 'grid.color': u'0.6',
                    'hatch.color': 'white',
                    'patch.linewidth': 0.5,
                    'font.size': 22,
                    'legend.fontsize': 25,
                    'ytick.labelsize': 22,
                    'lines.linewidth': 1.5,
                    'pdf.fonttype': 42,
                    }])

def create_df(n, tech, province_filter=None):
    """
    Create a heatmap-ready dataframe for a given storage technology.
    
    Parameters:
    -----------
    n : pypsa.Network
        PyPSA Network with simulation results
    tech : str
        Storage technology to analyse ('H2', 'battery', or 'water')
    province_filter : str, optional
        If given, only include storage assets in the selected province
    
    Returns:
    --------
    tuple
        (pivoted dataframe for heatmap, base power in MW)
    """
    # Select storage assets for the given technology
    stores = n.stores_t.p.filter(like=tech)
    
    # Apply province filter if given
    if province_filter:
        province_stores = stores.columns[stores.columns.str.contains(province_filter, case=False)]
        if len(province_stores) > 0:
            stores = stores[province_stores]
        else:
            print(f"Warning: No {tech} storage found in province {province_filter}")
            return pd.DataFrame(), 0
    
    # Compute maximum power as base for normalization
    base = abs(stores.sum(axis=1)).max()
    if base == 0:
        print(f"Warning: {tech} power is zero {'in province ' + province_filter if province_filter else 'nationwide'}")
        return pd.DataFrame(), 0
    
    # Normalize by base
    df = stores.sum(axis=1) / base
    df = df.to_frame()
    df.reset_index(inplace=True)
    renames = {0: 'p_store'}
    df.rename(columns=renames, inplace=True)
    
    # Timestamps are already in China local time
    date = n.stores_t.p.filter(like='water').index
    df['Hour'] = date.hour
    df['Day'] = date.strftime('%m-%d')
    
    # Build pivot table for heatmap
    summary = pd.pivot_table(data=df, index='Hour', columns='Day', values='p_store')
    summary = summary.fillna(0).infer_objects(copy=False)
    return summary, base

def create_aluminum_df(n, province_filter=None):
    """
    Create a heatmap-ready dataframe for aluminum smelter operation.
    
    Parameters:
    -----------
    n : pypsa.Network
        PyPSA Network with simulation results
    province_filter : str, optional
        If given, only include aluminum smelters in the selected province
    
    Returns:
    --------
    tuple
        (pivoted dataframe for heatmap, base power in MW)
    """
    # Aluminum smelter links (power from the electricity bus)
    aluminum_links = n.links_t.p0.filter(like='aluminum smelter')
    
    if aluminum_links.empty:
        print("Warning: No aluminum smelter links found in network")
        return pd.DataFrame(), 0
    
    # Apply province filter if given
    if province_filter:
        province_smelters = aluminum_links.columns[aluminum_links.columns.str.contains(province_filter, case=False)]
        if len(province_smelters) > 0:
            aluminum_links = aluminum_links[province_smelters]
        else:
            print(f"Warning: No aluminum smelters found in province {province_filter}")
            return pd.DataFrame(), 0
    
    # Compute maximum power as base for normalization
    base = abs(aluminum_links.sum(axis=1)).max()
    
    if base == 0:
        print(f"Warning: Aluminum smelter power is zero {'in province ' + province_filter if province_filter else 'nationwide'}")
        return pd.DataFrame(), 0
    
    # Normalize by base
    df = aluminum_links.sum(axis=1) / base
    df = df.to_frame()
    df.reset_index(inplace=True)
    renames = {0: 'p_smelter'}
    df.rename(columns=renames, inplace=True)
    
    # Timestamps are already in China local time
    date = aluminum_links.index
    df['Hour'] = date.hour
    df['Day'] = date.strftime('%m-%d')
    
    # Build pivot table for heatmap
    summary = pd.pivot_table(data=df, index='Hour', columns='Day', values='p_smelter')
    summary = summary.fillna(0).infer_objects(copy=False)
    return summary, base

def get_aluminum_storage_daily_average(n, province_filter=None):
    """
    Compute daily average aluminum storage to overlay on the heatmap.
    
    Parameters:
    -----------
    n : pypsa.Network
        PyPSA Network with simulation results
    province_filter : str, optional
        If given, only include storage in the selected province
    
    Returns:
    --------
    tuple
        (daily average storage level, minimum storage level for normalization)
    """
    # Aluminum storage time series
    aluminum_stores = n.stores_t.e.filter(like='aluminum storage')
    
    if aluminum_stores.empty:
        print("Warning: No aluminum storage found in network")
        return pd.Series(), 0
    
    # Apply province filter if given
    if province_filter:
        province_stores = aluminum_stores.columns[aluminum_stores.columns.str.contains(province_filter, case=False)]
        if len(province_stores) > 0:
            aluminum_stores = aluminum_stores[province_stores]
        else:
            print(f"Warning: No aluminum storage found in province {province_filter}")
            return pd.Series(), 0
    
    # Aggregate total storage
    total_storage = aluminum_stores.sum(axis=1)
    
    if total_storage.empty:
        print(f"Warning: Aluminum storage data is empty {'in province ' + province_filter if province_filter else 'nationwide'}")
        return pd.Series(), 0
    
    # Minimum storage level for normalization
    min_storage = total_storage.min()
    
    # Daily average storage
    daily_avg = total_storage.groupby(total_storage.index.strftime('%m-%d')).mean()
    
    # Normalize by subtracting the minimum
    daily_avg_normalized = daily_avg - min_storage
    
    return daily_avg_normalized, min_storage

def plot_comparison_heatmap(n_mmm, n_nmm, config, output_dir, tech, province_filter=None):
    """
    Generate comparison heatmaps for MMMU and NMMU scenarios.
    
    Parameters:
    -----------
    n_mmm : pypsa.Network
        Network object for MMMU scenario
    n_nmm : pypsa.Network
        Network object for NMMU scenario
    config : dict
        Plotting configuration dictionary
    output_dir : str
        Directory to save the heatmap figure
    tech : str
        Technology to plot
    province_filter : str, optional
        If given, only include data in the selected province
    """
    freq = config["freq"]
    planning_horizon = "2050"
    
    # Create province-specific subdirectory when filtering by province
    if province_filter:
        province_dir = os.path.join(output_dir, f"province_{province_filter}")
        os.makedirs(province_dir, exist_ok=True)
        plot_title_suffix = f" in {province_filter}"
    else:
        province_dir = output_dir
        plot_title_suffix = " (National)"
    
    # Create vertically stacked subplots (approx. square overall aspect ratio)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12))
    
    # MMMU scenario
    if tech == "aluminum":
        df_mmm, base_mmm = create_aluminum_df(n_mmm, province_filter)
        daily_storage_mmm, min_storage_mmm = get_aluminum_storage_daily_average(n_mmm, province_filter)
    else:
        df_mmm, base_mmm = create_df(n_mmm, tech, province_filter)
        daily_storage_mmm, min_storage_mmm = None, None
    
    # NMMU scenario
    if tech == "aluminum":
        df_nmm, base_nmm = create_aluminum_df(n_nmm, province_filter)
        daily_storage_nmm, min_storage_nmm = get_aluminum_storage_daily_average(n_nmm, province_filter)
    else:
        df_nmm, base_nmm = create_df(n_nmm, tech, province_filter)
        daily_storage_nmm, min_storage_nmm = None, None
    
    # Plot MMMU heatmap
    if not df_mmm.empty and base_mmm > 0:
        base_mmm_display = str(int(base_mmm / 1e3))  # display in GW
        
        if tech == "aluminum":
            sns.heatmap(df_mmm, ax=ax1, cmap='coolwarm', cbar_kws={'label': 'pu', 'shrink': 0.8}, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(df_mmm, ax=ax1, cmap='coolwarm', cbar_kws={'label': 'pu', 'shrink': 0.8}, vmin=-1.0, vmax=1.0)
        
        ax1.set_title(f'Mid smelter flexibility (core scenario)')
        # Hide x-axis label (Day) for the first subplot
        ax1.set_xlabel('')
        # ax1.set_xticklabels([])
        # Keep y-axis labels upright
        ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0)
        
        # Overlay aluminum storage level (without extra labels)
        if tech == "aluminum" and not daily_storage_mmm.empty:
            day_columns = df_mmm.columns
            storage_values = []
            storage_positions = []
            
            for i, day in enumerate(day_columns):
                if day in daily_storage_mmm.index:
                    storage_values.append(daily_storage_mmm[day]/1e6)  # convert to Mt
                    storage_positions.append(i + 0.5)  # column center
            
            if storage_values:
                ax1_twin = ax1.twinx()
                # Outline in white first
                ax1_twin.plot(storage_positions, storage_values, 'w-', linewidth=3.5, zorder=1)
                # Then draw black line inside
                ax1_twin.plot(storage_positions, storage_values, 'k-', linewidth=2, zorder=2)
                ax1_twin.set_ylabel('Stored aluminum (Mt)', color='black')
                ax1_twin.tick_params(axis='y', labelcolor='black')
                ax1_twin.set_xlim(0, len(day_columns))
    
    # Plot NMMU heatmap
    if not df_nmm.empty and base_nmm > 0:
        base_nmm_display = str(int(base_nmm / 1e3))  # display in GW
        
        if tech == "aluminum":
            sns.heatmap(df_nmm, ax=ax2, cmap='coolwarm', cbar_kws={'label': 'pu', 'shrink': 0.8}, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(df_nmm, ax=ax2, cmap='coolwarm', cbar_kws={'label': 'pu', 'shrink': 0.8}, vmin=-1.0, vmax=1.0)
        
        ax2.set_title(f'Non-constrained smelter flexibility')
        # Add x-axis label for the second subplot
        ax2.set_xlabel('Day')
        # Keep y-axis labels upright
        ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0)
        
        # Overlay aluminum storage level
        if tech == "aluminum" and not daily_storage_nmm.empty:
            day_columns = df_nmm.columns
            storage_values = []
            storage_positions = []
            
            for i, day in enumerate(day_columns):
                if day in daily_storage_nmm.index:
                    storage_values.append(daily_storage_nmm[day]/1e6)  # convert to Mt
                    storage_positions.append(i + 0.5)  # column center
            
            if storage_values:
                ax2_twin = ax2.twinx()
                # Outline in white first
                ax2_twin.plot(storage_positions, storage_values, 'w-', linewidth=3.5, zorder=1)
                # Then draw black line inside
                ax2_twin.plot(storage_positions, storage_values, 'k-', linewidth=2, label='Stored aluminum', zorder=2)
                ax2_twin.set_ylabel('Stored aluminum (Mt)', color='black')
                ax2_twin.tick_params(axis='y', labelcolor='black')
                ax2_twin.legend(loc='lower right', bbox_to_anchor=(1.1, -0.5))
                ax2_twin.set_xlim(0, len(day_columns))
    
    # Adjust layout and leave more space for legend
    plt.tight_layout()
    # Further adjust spacing between subplots for right-side legend
    plt.subplots_adjust(right=2)
    
    # plt.show()
    
    # Save figure
    output_path = os.path.join(province_dir, f"smelter_flexibility_comparison.png")
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved {tech} comparison heatmap to: {output_path}")
    plt.close()

def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description='Compare heatmaps between MMMU and NMMU scenarios')
    parser.add_argument('--config-mmm', type=str, 
                      default='configs/config_MMMU_2050_15p.yaml',
                      help='MMMU configuration file path')
    parser.add_argument('--config-nmm', type=str,
                      default='configs/config_NMMU_2050_15p.yaml', 
                      help='NMMU configuration file path')
    parser.add_argument('--output-dir', type=str,
                       default='results/comparison_heatmaps',
                       help='Output directory')
    parser.add_argument('--province', type=str, default=None,
                       help='Province filter (optional)')
    parser.add_argument('--techs', nargs='+', 
                       default=['H2', 'battery', 'water', 'aluminum'],
                       help='List of technologies to compare')
    
    args = parser.parse_args()
    
    # Set drawing style
    set_plot_style()
    
    # Load config files
    print("Loading configuration files...")
    with open(args.config_mmm, 'r', encoding='utf-8') as f:
        config_mmm = yaml.safe_load(f)
    
    with open(args.config_nmm, 'r', encoding='utf-8') as f:
        config_nmm = yaml.safe_load(f)
    
    # Plotting parameters (currently unused, kept for compatibility)
    map_figsize = config_mmm["plotting"]['map']['figsize']
    
    # Build network paths from version strings
    mmm_version = config_mmm['version']
    nmm_version = config_nmm['version']
    
    def _pick_network_path(version: str) -> str:
        candidates = [
            f"results/version-{version}/postnetworks/positive/postnetwork-ll-current+FCG-linear2050-2050.nc",
            f"results/version-{version}/postnetworks/positive/postnetwork-ll-current+Neighbor-linear2050-2050.nc",
        ]
        for p in candidates:
            if os.path.exists(p):
                return p
        return candidates[0]

    mmm_network_path = _pick_network_path(mmm_version)
    nmm_network_path = _pick_network_path(nmm_version)
    
    # Check network files
    if not os.path.exists(mmm_network_path):
        print(f"Error: MMMU network file not found: {mmm_network_path}")
        return
    
    if not os.path.exists(nmm_network_path):
        print(f"Error: NMMU network file not found: {nmm_network_path}")
        return
    
    # Load networks
    print("Loading network data...")
    n_mmm = pypsa.Network(mmm_network_path)
    n_nmm = pypsa.Network(nmm_network_path)
    
    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # --- Additional output: monthly aluminum production & storage for MMMU ---
    try:
        # Aluminum smelter power (from electricity bus)
        smelter_p = n_mmm.links_t.p0.filter(like='aluminum smelter')
        if not smelter_p.empty:
            # Use absolute power to avoid sign ambiguity, then aggregate by calendar month
            smelter_power_total = smelter_p.abs().sum(axis=1)
            monthly = smelter_power_total.to_frame(name="smelter_power_total_MW")
            monthly["year_month"] = monthly.index.to_period("M").astype(str)
            monthly_stats = (
                monthly.groupby("year_month")["smelter_power_total_MW"]
                .agg(["mean", "sum"])
                .rename(columns={"mean": "smelter_power_MW_mean", "sum": "smelter_energy_MWh"})
            )
            # Convert electricity use to aluminum production using 13.4 MWh per tonne
            # (assuming 1-hour resolution so sum of MW over the month is MWh)
            monthly_stats["smelter_production_tonnes"] = (
                monthly_stats["smelter_energy_MWh"] / 13.4
            )
            monthly_stats["smelter_production_Mt"] = (
                monthly_stats["smelter_production_tonnes"] / 1e6
            )
        else:
            monthly_stats = None

        # Aluminum storage (state of charge)
        storage = n_mmm.stores_t.e.filter(like='aluminum storage')
        if not storage.empty:
            total_storage = storage.sum(axis=1)
            # Convert to Mt consistent with the overlay (divide by 1e6)
            storage_df = total_storage.to_frame(name="storage_Mt")
            storage_df["storage_Mt"] = storage_df["storage_Mt"] / 1e6
            storage_df["year_month"] = storage_df.index.to_period("M").astype(str)
            storage_monthly = (
                storage_df.groupby("year_month")["storage_Mt"]
                .agg(["mean", "min", "max"])
                .rename(
                    columns={
                        "mean": "storage_Mt_mean",
                        "min": "storage_Mt_min",
                        "max": "storage_Mt_max",
                    }
                )
            )
        else:
            storage_monthly = None

        if monthly_stats is not None or storage_monthly is not None:
            # Merge production and storage stats on year_month
            if monthly_stats is None:
                monthly_out = storage_monthly
            elif storage_monthly is None:
                monthly_out = monthly_stats
            else:
                monthly_out = monthly_stats.join(storage_monthly, how="outer")

            csv_path = os.path.join(
                args.output_dir, "mmmu_2050_monthly_aluminum_production_and_storage.csv"
            )
            monthly_out.to_csv(csv_path)
            print(f"Saved MMMU monthly aluminum production & storage to: {csv_path}")
    except Exception as e:
        print(f"Warning: failed to save MMMU monthly aluminum production/storage CSV: {e}")

    # Generate comparison heatmaps
    print("Generating comparison heatmaps...")
    for tech in args.techs:
        print(f"Processing {tech} technology...")
        plot_comparison_heatmap(n_mmm, n_nmm, config_mmm, args.output_dir, tech, args.province)
    
    print("All comparison heatmaps generated successfully!")

if __name__ == "__main__":
    main()
