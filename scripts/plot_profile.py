"""
This script generates weekly system operation plots showing various resource outputs during heating and non-heating periods.
It creates two types of plots:
1. Weekly system operation plots showing renewable generation, conventional generation, storage technologies, and load
2. Heating system comparison plots showing heating demand and supply during heating vs non-heating periods

The plots show how different energy resources operate throughout one week, with time on the x-axis and power on the y-axis.
For storage technologies, a positive value indicates the process of discharging, whereas a negative value signifies the charging process.
"""

from _helpers import configure_logging
import seaborn as sns
import pandas as pd
import pypsa
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta

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

def plot_weekly_system_operation(n, config):
    """
    Generates weekly system operation plots showing various resource outputs.
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object containing the simulation results
    config : dict
        Configuration dictionary containing plotting parameters
    """
    planning_horizon = snakemake.wildcards.planning_horizons
    
    # Define time periods for analysis
    # Heating period typically ends around March 15th in China
    # We'll analyze the first week after heating ends (around April 1st) and a week during heating (around January 15th)
    
    # Get the full time index
    time_index = n.stores_t.p.index
    
    # Find April 1st week (first week after heating ends)
    april_start = None
    for date in time_index:
        if date.month == 4 and date.day <= 7:
            april_start = date
            break
    
    # Find January 15th week (during heating period)
    january_start = None
    for date in time_index:
        if date.month == 1 and 10 <= date.day <= 17:
            january_start = date
            break
    
    if april_start is None or january_start is None:
        print("Warning: Could not find appropriate dates for weekly analysis")
        return
    
    # Define one week (168 hours)
    week_hours = 168
    
    # Get data for both periods
    periods = {
        'Non-heating period (April 1st week)': april_start,
        'Heating period (January 15th week)': january_start
    }
    
    for period_name, start_date in periods.items():
        # Get one week of data
        end_date = start_date + timedelta(hours=week_hours-1)
        mask = (time_index >= start_date) & (time_index <= end_date)
        week_data = time_index[mask]
        
        if len(week_data) < week_hours:
            print(f"Warning: Insufficient data for {period_name}")
            continue
        
        # Create subplots for different resource types
        fig, axes = plt.subplots(4, 1, figsize=(15, 12))
        fig.suptitle(f'System Operation Throughout 1 Week - {period_name} ({planning_horizon})', fontsize=16)
        
        # 1. Renewable Generation
        ax1 = axes[0]
        if hasattr(n, 'generators_t') and hasattr(n.generators_t, 'p'):
            # Solar PV
            solar_pv = n.generators_t.p.filter(like='solar').sum(axis=1)[week_data]
            ax1.plot(week_data, solar_pv, label='Solar PV', color='orange', linewidth=2)
            
            # Wind (onshore + offshore)
            wind_onshore = n.generators_t.p.filter(like='onwind').sum(axis=1)[week_data]
            wind_offshore = n.generators_t.p.filter(like='offwind').sum(axis=1)[week_data]
            ax1.plot(week_data, wind_onshore, label='Wind Onshore', color='lightblue', linewidth=2)
            ax1.plot(week_data, wind_offshore, label='Wind Offshore', color='blue', linewidth=2)
            
            # Hydro
            hydro = n.generators_t.p.filter(like='hydro').sum(axis=1)[week_data]
            ax1.plot(week_data, hydro, label='Hydro', color='darkblue', linewidth=2)
        
        ax1.set_ylabel('Power (MW)')
        ax1.set_title('Renewable Generation')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Conventional Generation
        ax2 = axes[1]
        if hasattr(n, 'generators_t') and hasattr(n.generators_t, 'p'):
            # Coal
            coal = n.generators_t.p.filter(like='coal').sum(axis=1)[week_data]
            ax2.plot(week_data, coal, label='Coal', color='black', linewidth=2)
            
            # Gas
            gas = n.generators_t.p.filter(like='gas').sum(axis=1)[week_data]
            ax2.plot(week_data, gas, label='Gas', color='red', linewidth=2)
            
            # Nuclear
            nuclear = n.generators_t.p.filter(like='nuclear').sum(axis=1)[week_data]
            ax2.plot(week_data, nuclear, label='Nuclear', color='purple', linewidth=2)
            
            # Biomass
            biomass = n.generators_t.p.filter(like='biomass').sum(axis=1)[week_data]
            ax2.plot(week_data, biomass, label='Biomass', color='green', linewidth=2)
        
        ax2.set_ylabel('Power (MW)')
        ax2.set_title('Conventional Generation')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Storage Technologies
        ax3 = axes[2]
        if hasattr(n, 'stores_t') and hasattr(n.stores_t, 'p'):
            # Battery storage (positive = discharging, negative = charging)
            battery = n.stores_t.p.filter(like='battery').sum(axis=1)[week_data]
            ax3.plot(week_data, battery, label='Battery Storage', color='green', linewidth=2)
            
            # Hydrogen storage
            h2 = n.stores_t.p.filter(like='H2').sum(axis=1)[week_data]
            ax3.plot(week_data, h2, label='H2 Storage', color='cyan', linewidth=2)
            
            # Water storage
            water = n.stores_t.p.filter(like='water').sum(axis=1)[week_data]
            ax3.plot(week_data, water, label='Water Storage', color='blue', linewidth=2)
        
        ax3.set_ylabel('Power (MW)')
        ax3.set_title('Storage Technologies (Positive = Discharging, Negative = Charging)')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        ax3.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        
        # 4. Load and Net Load
        ax4 = axes[3]
        if hasattr(n, 'loads_t') and hasattr(n.loads_t, 'p'):
            # Total load
            total_load = n.loads_t.p.sum(axis=1)[week_data]
            ax4.plot(week_data, total_load, label='Total Load', color='red', linewidth=2)
            
            # Net load (load minus renewables)
            if hasattr(n, 'generators_t') and hasattr(n.generators_t, 'p'):
                renewables = (n.generators_t.p.filter(like='solar').sum(axis=1) + 
                            n.generators_t.p.filter(like='wind').sum(axis=1) + 
                            n.generators_t.p.filter(like='hydro').sum(axis=1))[week_data]
                net_load = total_load - renewables
                ax4.plot(week_data, net_load, label='Net Load', color='orange', linewidth=2)
        
        ax4.set_ylabel('Power (MW)')
        ax4.set_title('Load and Net Load')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        # Format x-axis
        for ax in axes:
            ax.set_xlim(week_data[0], week_data[-1])
            # Show only every 24 hours (daily)
            ax.set_xticks(week_data[::24])
            ax.set_xticklabels([d.strftime('%m-%d') for d in week_data[::24]], rotation=45)
        
        # Save the plot
        period_suffix = 'non_heating' if 'April' in period_name else 'heating'
        output_file = snakemake.output[f"weekly_operation_{period_suffix}"]
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()

def plot_heating_comparison(n, config):
    """
    Generates a comparison plot showing heating demand and supply during heating vs non-heating periods.
    
    Parameters:
    -----------
    n : pypsa.Network
        The PyPSA network object containing the simulation results
    config : dict
        Configuration dictionary containing plotting parameters
    """
    planning_horizon = snakemake.wildcards.planning_horizons
    
    # Get the full time index
    time_index = n.stores_t.p.index
    
    # Find representative weeks
    april_start = None
    january_start = None
    
    for date in time_index:
        if date.month == 4 and date.day <= 7:
            april_start = date
        if date.month == 1 and 10 <= date.day <= 17:
            january_start = date
        if april_start and january_start:
            break
    
    if april_start is None or january_start is None:
        print("Warning: Could not find appropriate dates for heating comparison")
        return
    
    # Get one week of data for each period
    week_hours = 168
    
    fig, axes = plt.subplots(2, 1, figsize=(15, 10))
    fig.suptitle(f'Heating System Comparison - {planning_horizon}', fontsize=16)
    
    periods = [
        (april_start, 'Non-heating period (April)', axes[0]),
        (january_start, 'Heating period (January)', axes[1])
    ]
    
    for start_date, period_name, ax in periods:
        end_date = start_date + timedelta(hours=week_hours-1)
        mask = (time_index >= start_date) & (time_index <= end_date)
        week_data = time_index[mask]
        
        if len(week_data) < week_hours:
            continue
        
        # Heating demand (if available)
        if hasattr(n, 'loads_t') and hasattr(n.loads_t, 'p'):
            heating_load = n.loads_t.p.filter(like='heat').sum(axis=1)[week_data]
            ax.plot(week_data, heating_load, label='Heating Demand', color='red', linewidth=2)
        
        # Heat supply from different sources
        if hasattr(n, 'links_t') and hasattr(n.links_t, 'p1'):
            # CHP heat
            chp_heat = n.links_t.p1.filter(like='CHP').sum(axis=1)[week_data]
            ax.plot(week_data, chp_heat, label='CHP Heat', color='orange', linewidth=2)
            
            # Heat pump heat
            heat_pump = n.links_t.p1.filter(like='heat pump').sum(axis=1)[week_data]
            ax.plot(week_data, heat_pump, label='Heat Pump', color='blue', linewidth=2)
            
            # Boiler heat
            boiler = n.links_t.p1.filter(like='boiler').sum(axis=1)[week_data]
            ax.plot(week_data, boiler, label='Boiler', color='green', linewidth=2)
        
        # Water storage for heat
        if hasattr(n, 'stores_t') and hasattr(n.stores_t, 'p'):
            water_heat = n.stores_t.p.filter(like='water').sum(axis=1)[week_data]
            ax.plot(week_data, water_heat, label='Water Storage (Heat)', color='cyan', linewidth=2)
        
        ax.set_ylabel('Power (MW)')
        ax.set_title(period_name)
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(week_data[0], week_data[-1])
        ax.set_xticks(week_data[::24])
        ax.set_xticklabels([d.strftime('%m-%d') for d in week_data[::24]], rotation=45)
    
    # Save the plot
    output_file = snakemake.output["heating_comparison"]
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    # Set up mock snakemake for testing if not running in snakemake
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_profile',
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
    
    # Generate weekly system operation plots
    plot_weekly_system_operation(n, config)
    plot_heating_comparison(n, config)