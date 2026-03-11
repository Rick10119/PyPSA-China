#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simplified plotting script for the MMMU-2050 scenario cost analysis.

The script reads the pre-processed file `mmmu_2050_M_detailed_data.csv` and
produces a figure with:
  - x-axis: aluminum smelting capacity (5p–100p, converted to Mt/year)
  - left y-axis: cost savings / increases (billion CNY)

It displays electricity-system cost savings, aluminum operational cost changes,
and net cost savings.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import logging

# Use fonts that render Latin characters cleanly across platforms
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def load_detailed_data(csv_path):
    """
    Load the detailed MMMU-2050 data from CSV.

    Parameters
    ----------
    csv_path : str or Path
        Path to the CSV file.

    Returns
    -------
    pd.DataFrame
        Data frame with scenario metrics.
    """
    try:
        df = pd.read_csv(csv_path)
        logger.info(f"Successfully loaded data file: {csv_path}")
        logger.info(f"Data shape: {df.shape}")
        return df
    except Exception as e:
        logger.error(f"Failed to load data file: {str(e)}")
        return None

def plot_mmm_2050_from_csv(csv_path, output_dir=None):
    """
    Plot MMMU-2050 cost and emissions metrics from a prepared CSV file.

    Parameters
    ----------
    csv_path : str or Path
        Path to the CSV file created by the analysis workflow.
    output_dir : str or Path, optional
        Output directory for the figure (defaults to the CSV parent directory).
    """
    # Load data
    df = load_detailed_data(csv_path)
    if df is None:
        return
    
    # Configure output directory
    if output_dir is None:
        output_dir = Path(csv_path).parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    # Prepare data arrays
    x = df['Capacity_Value_Mt'].values
    power_savings = df['Power_Cost_Changes_Billion_CNY'].values
    aluminum_changes = df['Aluminum_Cost_Changes_Billion_CNY'].values
    net_savings = df['Net_Cost_Savings_Billion_CNY'].values
    
    # Set bar width
    bar_width = 150
    
    # Use slightly offset x-positions so electricity and aluminum bars do not overlap
    x_power = [pos - bar_width/6 for pos in x]      # electricity-system costs shifted left
    x_aluminum = [pos + bar_width/6 for pos in x]   # aluminum operational costs shifted right
    
    # Plot electricity-system cost savings (bottom bars)
    bars1 = ax.bar(x_power, power_savings, bar_width*0.8, color='#1f77b4', alpha=0.8, 
                   label='Power System Cost Savings')
    
    # Plot aluminum operational cost changes (stacked and slightly offset)
    bars2 = ax.bar(x_aluminum, aluminum_changes, bar_width*0.8, bottom=power_savings, 
                   color='#ff7f0e', alpha=0.8, label='Aluminum Operation Cost Increase')
    
    # Plot net cost savings curve (black line at original x positions)
    ax.plot(x, net_savings, 'k-', linewidth=3, label='Net Cost Savings', marker='o', markersize=8, zorder=20)
    
    # Find the point with maximum net savings
    max_saving_index = np.argmax(net_savings)
    max_saving_value = net_savings[max_saving_index]
    max_saving_capacity = x[max_saving_index]
    
    # Mark the maximum net-savings point with a star
    ax.plot(max_saving_capacity, max_saving_value, 'r*', markersize=15, 
            label=f'Highest Net Savings: {max_saving_value:.0f}B CNY', zorder=30)
    
    # Add numeric label for the maximum net-savings point
    ax.annotate(f'{max_saving_value:.0f}B',
                xy=(max_saving_capacity, max_saving_value),
                xytext=(0, 20),
                textcoords="offset points",
                ha='center', va='bottom', 
                fontsize=20, weight='bold', color='red',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Axis labels
    ax.set_xlabel('Aluminum Smelting Capacity (Mt/Year)', fontsize=20)
    ax.set_ylabel('Cost Savings/Increase (Billion CNY)', fontsize=20, color='blue')
    
    # Add light grid
    ax.grid(True, alpha=0.3, axis='y')
    
    # Configure x-axis ticks and labels
    ax.set_xticks(x)
    ax.set_xticklabels([f'{cap/100:.0f}' for cap in x], fontsize=20)
    
    # Configure left y-axis tick labels
    y_ticks = ax.get_yticks()
    y_tick_labels = [f'{tick:.0f}' for tick in y_ticks]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels, fontsize=20, color='blue')
    
    # Build legend elements
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor='#1f77b4', alpha=0.8, label='Electricity System Cost Savings'),
        plt.Rectangle((0,0),1,1, facecolor='#ff7f0e', alpha=0.8, label='Smelter Operational Cost Increase'),
        plt.Line2D([0], [0], color='black', linewidth=3, marker='o', markersize=8, label='Net Benefit'),
        plt.Line2D([0], [0], marker='*', color='red', markersize=15, linestyle='', label='Highest Net Benefit'),
    ]
    
    # Add legend
    ax.legend(handles=legend_elements, loc='lower left', fontsize=20)
    
    plt.tight_layout()
    
    # Save figure
    plot_file = output_dir / "mmmu_2050_analysis.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"MMMU-2050 scenario analysis figure saved to: {plot_file}")
    
    
    return fig, ax

def main():
    """CLI entry point."""
    # Default CSV file path
    default_csv_path = "results/mmmu_2050_analysis/mmmu_2050_M_detailed_data.csv"
    
    # Check that the CSV file exists
    csv_path = Path(default_csv_path)
    if not csv_path.exists():
        # Backward-compatible fallback to the old naming convention
        legacy_csv_path = Path("results/mmm_2050_analysis/mmm_2050_M_detailed_data.csv")
        if legacy_csv_path.exists():
            logger.warning(f"Default CSV not found; falling back to legacy path: {legacy_csv_path}")
            csv_path = legacy_csv_path
        else:
            logger.error(f"CSV file not found: {csv_path}")
            logger.info("Please ensure that `mmmu_2050_M_detailed_data.csv` (or legacy `mmm_2050_M_detailed_data.csv`) is available.")
            return
    
    logger.info("Start plotting MMMU-2050 scenario analysis figure.")
    logger.info(f"Data file: {csv_path}")
    
    # Generate the plot
    try:
        fig, ax = plot_mmm_2050_from_csv(csv_path)
        logger.info("Plotting complete.")
        
        # Optionally show the figure in interactive environments
        # plt.show()
        
    except Exception as e:
        logger.error(f"Error while plotting MMMU-2050 analysis figure: {str(e)}")
        raise

if __name__ == "__main__":
    main()
