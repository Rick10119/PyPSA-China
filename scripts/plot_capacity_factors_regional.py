"""
This script generates regional capacity factor plots for specific provinces (Xinjiang, Yunnan, InnerMongolia, Guangxi, Shandong)
in the PyPSA-China model. It creates visualizations showing how capacity factors vary by month for different technologies.
"""

from _helpers import configure_logging
import seaborn as sns
import pandas as pd
import pypsa
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

def set_plot_style():
    """
    Sets up the plotting style for all matplotlib plots in this script.
    """
    plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'sans-serif']
    plt.rcParams['axes.unicode_minus'] = False
    plt.style.use(['classic', 'seaborn-v0_8-whitegrid',
                   {'axes.grid': False, 'grid.linestyle': '--', 'grid.color': u'0.6',
                    'hatch.color': 'white',
                    'patch.linewidth': 0.5,
                    'font.size': 18,
                    'legend.fontsize': 'large',
                    'lines.linewidth': 1.5,
                    'pdf.fonttype': 42,
                    }])

def filter_network_by_province(n, target_province=None):
    """
    Filter the network to include only components from a specific province.
    """
    if target_province is None:
        return n
    
    print(f"正在过滤网络，只保留 {target_province} 省份的组件...")
    
    # Create a copy of the network to avoid modifying the original
    n_filtered = n.copy()
    
    # Find buses in the target province
    province_buses = n_filtered.buses[n_filtered.buses.index.str.contains(target_province, case=False)].index
    
    if len(province_buses) == 0:
        print(f"警告：未找到 {target_province} 省份的节点")
        return n_filtered
    
    print(f"找到 {len(province_buses)} 个 {target_province} 省份的节点")
    
    # Remove generators not in the target province
    non_province_generators = n_filtered.generators[~n_filtered.generators.bus.isin(province_buses)].index
    if len(non_province_generators) > 0:
        n_filtered.mremove("Generator", non_province_generators)
    
    # Remove loads not in the target province
    non_province_loads = n_filtered.loads[~n_filtered.loads.bus.isin(province_buses)].index
    if len(non_province_loads) > 0:
        n_filtered.mremove("Load", non_province_loads)
    
    # Remove storage units not in the target province
    non_province_storage = n_filtered.storage_units[~n_filtered.storage_units.bus.isin(province_buses)].index
    if len(non_province_storage) > 0:
        n_filtered.mremove("StorageUnit", non_province_storage)
    
    # Remove stores not in the target province
    non_province_stores = n_filtered.stores[~n_filtered.stores.bus.isin(province_buses)].index
    if len(non_province_stores) > 0:
        n_filtered.mremove("Store", non_province_stores)
    
    # Remove links not connected to the target province
    non_province_links = n_filtered.links[~(n_filtered.links.bus0.isin(province_buses) | n_filtered.links.bus1.isin(province_buses))].index
    if len(non_province_links) > 0:
        n_filtered.mremove("Link", non_province_links)
    
    # Remove lines not connected to the target province
    non_province_lines = n_filtered.lines[~(n_filtered.lines.bus0.isin(province_buses) | n_filtered.lines.bus1.isin(province_buses))].index
    if len(non_province_lines) > 0:
        n_filtered.mremove("Line", non_province_lines)
    
    # Finally remove non-province buses
    non_province_buses = n_filtered.buses[~n_filtered.buses.index.isin(province_buses)].index
    if len(non_province_buses) > 0:
        n_filtered.mremove("Bus", non_province_buses)
    
    return n_filtered

def calculate_monthly_capacity_factors(n):
    """
    Calculate monthly average capacity factors for all generators and power-producing links in the network.
    """
    if not hasattr(n, 'generators_t') or not hasattr(n.generators_t, 'p'):
        print("Warning: No generator time series data found")
        return {}
    
    # Get generator power output and calculate actual maximum power output
    gen_power = n.generators_t.p
    gen_max_power = gen_power.max()
    
    # Get link power output and calculate actual maximum power output
    link_power = pd.DataFrame()
    link_max_power = pd.Series(dtype=float)
    
    if hasattr(n, 'links_t') and hasattr(n.links_t, 'p0'):
        # Filter links that produce electricity (bus1 is electricity bus)
        elec_links = n.links[n.links.bus1.isin(n.buses[n.buses.carrier == 'AC'].index)]
        if not elec_links.empty:
            link_power = n.links_t.p0[elec_links.index].copy()
            link_max_power = link_power.max()
        
        # Also include aluminum smelters
        aluminum_smelters = n.links[n.links.carrier == 'aluminum']
        if not aluminum_smelters.empty:
            for link in aluminum_smelters.index:
                if link in n.links_t.p0.columns:
                    if link not in link_power.columns:
                        link_power[link] = n.links_t.p0[link]
                    if link not in link_max_power.index:
                        link_max_power[link] = n.links_t.p0[link].max()

    # Define technology groups
    tech_groups = {
        'Hydro': ['hydro', 'hydroelectricity'],
        'Nuclear': ['nuclear'],
        'Coal': ['coal cc', 'CHP coal', 'coal power plant'],
        'Gas': ['OCGT gas', 'CHP gas'],
        'Wind': ['onwind', 'offwind', 'wind'],
        'Solar': ['solar', 'solar pv', 'pv'],
        'Aluminum': ['aluminum', 'smelter'],
        'Other': []
    }
    
    monthly_cf = {}
    
    for group_name, carriers in tech_groups.items():
        # Find generators and links belonging to this group
        if carriers:
            group_generators = []
            group_links = []
            
            for carrier in carriers:
                # Generators matching
                exact_matches = n.generators[n.generators.carrier == carrier].index.tolist()
                exact_matches = [gen for gen in exact_matches if 'fuel' not in gen.lower()]
                group_generators.extend(exact_matches)
                
                if not exact_matches:
                    partial_matches = n.generators[n.generators.carrier.str.contains(carrier, case=False, na=False)].index.tolist()
                    partial_matches = [gen for gen in partial_matches if 'fuel' not in gen.lower()]
                    group_generators.extend(partial_matches)
                
                # Links matching
                if not link_power.empty:
                    exact_link_matches = n.links[n.links.carrier == carrier].index.tolist()
                    exact_link_matches = [link for link in exact_link_matches if link in link_power.columns]
                    group_links.extend(exact_link_matches)
                    
                    if not exact_link_matches:
                        partial_link_matches = n.links[n.links.carrier.str.contains(carrier, case=False, na=False)].index.tolist()
                        partial_link_matches = [link for link in partial_link_matches if link in link_power.columns]
                        group_links.extend(partial_link_matches)
                    
                    if carrier in ['aluminum', 'smelter']:
                        smelter_links = [link for link in link_power.columns if 'smelter' in link.lower()]
                        group_links.extend(smelter_links)
        else:
            # 'Other' group logic
            all_used_generators = set()
            all_used_links = set()
            for carriers_list in tech_groups.values():
                if carriers_list:
                    for carrier in carriers_list:
                        all_used_generators.update(n.generators[n.generators.carrier == carrier].index.tolist())
                        all_used_generators.update(n.generators[n.generators.carrier.str.contains(carrier, case=False, na=False)].index.tolist())
                        if not link_power.empty:
                            all_used_links.update(n.links[n.links.carrier == carrier].index.tolist())
                            all_used_links.update(n.links[n.links.carrier.str.contains(carrier, case=False, na=False)].index.tolist())
            
            group_generators = [gen for gen in n.generators.index if gen not in all_used_generators and 'fuel' not in gen.lower()]
            group_links = [link for link in link_power.columns if link not in all_used_links and 'smelter' not in link.lower()]
        
        # Calculate CF
        total_power = pd.Series(0, index=gen_power.index)
        total_max_power = 0
        
        if group_generators:
            total_power += gen_power[group_generators].sum(axis=1)
            total_max_power += gen_max_power[group_generators].sum()
        
        if group_links:
            total_power += link_power[group_links].sum(axis=1)
            total_max_power += link_max_power[group_links].sum()
        
        if total_max_power > 0:
            cf = total_power / total_max_power
            cf_df = cf.to_frame('cf')
            cf_df['month'] = cf_df.index.month
            monthly_cf[group_name] = cf_df.groupby('month')['cf'].mean()
    
    return monthly_cf

def calculate_monthly_load_factors(n):
    """
    Calculate monthly average load factors.
    """
    monthly_load = {}
    if hasattr(n, 'loads_t') and hasattr(n.loads_t, 'p_set'):
        # Electricity Load
        elec_loads = n.loads_t.p_set.filter(regex='^(?!.*(heat|aluminum)).*$', axis=1)
        if not elec_loads.empty:
            total_elec_load = elec_loads.sum(axis=1)
            if total_elec_load.max() > 0:
                load_df = (total_elec_load / total_elec_load.max()).to_frame('load_factor')
                load_df['month'] = load_df.index.month
                monthly_load['Electricity Load'] = load_df.groupby('month')['load_factor'].mean()
        
        # Heating Load
        heat_loads = n.loads_t.p_set.filter(like='heat')
        if not heat_loads.empty:
            total_heat_load = heat_loads.sum(axis=1)
            if total_heat_load.max() > 0:
                load_df = (total_heat_load / total_heat_load.max()).to_frame('load_factor')
                load_df['month'] = load_df.index.month
                monthly_load['Heating Load'] = load_df.groupby('month')['load_factor'].mean()
                
    return monthly_load

def save_monthly_data_to_csv(monthly_cf, monthly_load, planning_horizon, target_province):
    """
    Save results to CSV.
    """
    output_dir = "results/monthly_capacity_factors_regional"
    os.makedirs(output_dir, exist_ok=True)
    
    csv_filename = f"{output_dir}/monthly_cf_{planning_horizon}_{target_province}.csv"
    
    all_data = {}
    for tech, cf_data in monthly_cf.items():
        if not cf_data.empty: all_data[f"{tech}_CF_Avg"] = cf_data
    for load_type, load_data in monthly_load.items():
        if not load_data.empty: all_data[f"{load_type}_Load_Factor"] = load_data
        
    if all_data:
        df = pd.DataFrame(all_data)
        df.index.name = 'Month'
        month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        df['Month_Name'] = [month_names[i-1] for i in df.index]
        df = df[['Month_Name'] + [col for col in df.columns if col != 'Month_Name']]
        df.to_csv(csv_filename)
        print(f"Data saved to {csv_filename}")

def plot_province_capacity_factors(n, planning_horizon, target_province):
    """
    Generate plot for a specific province.
    """
    n_prov = filter_network_by_province(n, target_province)
    
    monthly_cf = calculate_monthly_capacity_factors(n_prov)
    monthly_load = calculate_monthly_load_factors(n_prov)
    
    if not monthly_cf and not monthly_load:
        print(f"Warning: No data for {target_province}")
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9))
    colors = {
        'Hydro': '#0000FF',      # Blue
        'Nuclear': '#800080',    # Purple
        'Coal': '#000000',       # Black
        'Gas': '#FF0000',        # Red
        'Wind': '#00FF00',       # Green
        'Solar': '#FFD700',      # Gold
        'Aluminum': '#FF69B4',   # Hot pink
        'Other': '#808080'       # Gray
    }
    load_colors = {
        'Electricity Load': '#1f77b4',    # Blue
        'Heating Load': '#ff7f0e',        # Orange
        'Aluminum Load': '#2ca02c'        # Green
    }
    display_labels = {
        'Heating Load': 'Heating demand',
        'Aluminum Load': 'Aluminum smelter',
        'Electricity Load': 'Electricity load'
    }

    load_types_upper = ['Electricity Load', 'Heating Load']
    for load_type in monthly_load.keys():
        if load_type in load_types_upper:
            months = monthly_load[load_type].index
            values = monthly_load[load_type].values
            color = load_colors.get(load_type, '#000000')
            display_label = display_labels.get(load_type, load_type)
            ax1.plot(months, values, 's--', color=color,
                    linewidth=4, markersize=6, label=display_label)

    if 'Aluminum' in monthly_cf:
        months = monthly_cf['Aluminum'].index
        values = monthly_cf['Aluminum'].values
        color = colors.get('Aluminum', '#000000')
        ax1.plot(months, values, 'o-', color=color,
                linewidth=4, markersize=6, label='Aluminum smelter')

    tech_types_lower = ['Hydro', 'Coal', 'Gas', 'Wind', 'Solar']
    for tech in monthly_cf.keys():
        if tech in tech_types_lower:
            months = monthly_cf[tech].index
            values = monthly_cf[tech].values
            color = colors.get(tech, '#000000')
            ax2.plot(months, values, 'o-', color=color,
                    linewidth=4, markersize=6, label=tech)

    ax1.set_ylabel('Capacity Factor', fontsize=30)
    ax1.set_title('Monthly Load & Smelter Capacity Factors', fontsize=30)
    ax1.set_xlim(1.0, 12.0)
    ax1.set_ylim(-0.005, 1.005)
    ax1.set_xticks(range(1, 13))
    ax1.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                         'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], fontsize=30)
    ax1.tick_params(axis='y', labelsize=30)
    ax1.grid(True, alpha=0.3)

    ax2.set_ylabel('Capacity Factor', fontsize=30)
    ax2.set_xlabel('Month', fontsize=30)
    ax2.set_title('Monthly Generation Capacity Factors', fontsize=30)
    ax2.set_xlim(1.0, 12.0)
    ax2.set_ylim(-0.0025, 0.8025)
    ax2.set_xticks(range(1, 13))
    ax2.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                         'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'], fontsize=30)
    ax2.tick_params(axis='y', labelsize=30)
    ax2.grid(True, alpha=0.3)

    ax1.legend(loc='center left', bbox_to_anchor=(1.04, 0.5),
               ncol=1, fontsize=30, borderaxespad=0.)
    legend_order = ['Hydro', 'Wind', 'Solar', 'Gas', 'Coal']
    handles, labels = ax2.get_legend_handles_labels()
    label_handle_map = dict(zip(labels, handles))
    ordered_handles = [label_handle_map[l] for l in legend_order if l in label_handle_map]
    ordered_labels = [l for l in legend_order if l in label_handle_map]
    ax2.legend(ordered_handles, ordered_labels, loc='center left',
               bbox_to_anchor=(1.04, 0.5), ncol=1, fontsize=30, borderaxespad=0.)
    
    plt.tight_layout()
    plt.subplots_adjust(right=2)
    
    output_dir = "results/plots/capacity_factors_regional"
    os.makedirs(output_dir, exist_ok=True)
    plot_filename = f"{output_dir}/cf_{planning_horizon}_{target_province}.png"
    fig.savefig(plot_filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {plot_filename}")
    
    save_monthly_data_to_csv(monthly_cf, monthly_load, planning_horizon, target_province)

if __name__ == "__main__":
    # Hardcoded scenario as requested by the user
    scenario = "version-0120.1H.1-MMMU-2050-10p"
    planning_horizon = "2050"
    network_path = f"results/{scenario}/postnetworks/positive/postnetwork-ll-current+FCG-linear2050-2050.nc"
    
    if 'snakemake' not in globals():
        from types import SimpleNamespace
        class MockSnakemake:
            def __init__(self):
                self.input = SimpleNamespace(network=network_path)
                self.wildcards = SimpleNamespace(planning_horizons=planning_horizon)
                self.config = {}
                self.rule = 'plot_capacity_factors_regional'
                self.log = self.MockLog()
            
            class MockLog:
                def __init__(self):
                    self.python = "logs/plot_capacity_factors_regional.log"
                def get(self, key, default):
                    val = getattr(self, key, None)
                    return val if val is not None else default
                def __getitem__(self, key):
                    return self.python if key == 0 or key == 'python' else None
                def __iter__(self):
                    return iter([self.python])
                def __len__(self):
                    return 1
                def __bool__(self):
                    return True
        
        snakemake = MockSnakemake()
    else:
        planning_horizon = snakemake.wildcards.planning_horizons  # type: ignore[name-defined]
        network_path = snakemake.input.network  # type: ignore[name-defined]
    
    configure_logging(snakemake)
    set_plot_style()
    
    print(f"\nScenario: {scenario}")
    print(f"Loading network from: {network_path}")
    
    if not os.path.exists(network_path):
        print(f"Error: Network file not found at {network_path}")
        # Try relative path if absolute fails, though results/ is usually in root
        if not os.path.exists(network_path):
            import sys
            print("Current working directory:", os.getcwd())
            sys.exit(1)

    n = pypsa.Network(network_path)
    
    # List of provinces for regional analysis as requested by the user
    provinces = ["Xinjiang", "Yunnan", "InnerMongolia", "Guangxi", "Shandong"]
    
    for province in provinces:
        print(f"\nProcessing province: {province}")
        plot_province_capacity_factors(n, planning_horizon, province)

    print("\nRegional analysis completed.")
