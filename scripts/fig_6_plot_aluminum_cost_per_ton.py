#!/usr/bin/env python3
"""
Visualization script for the cost composition of primary aluminum smelting.

It shows how the levelized cost per tonne is decomposed into cost components
for 2020 and several 2050 capacity-ratio cases.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
import os

# Use Helvetica-like fonts for clear English labels
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def create_aluminum_cost_bar_chart():
    """Create a stacked bar chart of aluminum cost composition."""
    
    # Data preparation
    categories = ['Raw Materials', 'Labor', 'Fixed o&m', 'Restart', 'Depreciation', 'Retirement loss', 'Storage', 'Electricity']
    costs_2020 = [8451.2, 150, 1900, 0.01, 300, 0, 1, 6250]
    costs_2050_5p = [8451.2, 159.889114, 2171.893745, 769.1493874, 347.0, 809.8, 2.4, 4529.662849]
    costs_2050_20p = [8451.2, 159.1901429, 2984.318766, 299.1122492, 462.1, 692.7, 15.6, 2476.284987]
    costs_2050_100p = [8451.2, 160.7041192, 7326.478149, 707.9584908, 1156.812339, 0, 25.2, 797.0111764]
    
    # Create DataFrame with cost breakdown per tonne of primary aluminum
    df = pd.DataFrame({
        'Category': categories,
        '2020 cost (CNY/tonne)': costs_2020,
        '2050_5p cost (CNY/tonne)': costs_2050_5p,
        '2050_20p cost (CNY/tonne)': costs_2050_20p,
        '2050_100p cost (CNY/tonne)': costs_2050_100p
    })
    
    # Compute total levelized cost per tonne for each case
    total_2020 = sum(costs_2020)
    total_2050_5p = sum(costs_2050_5p)
    total_2050_20p = sum(costs_2050_20p)
    total_2050_100p = sum(costs_2050_100p)
    
    print(f"Total cost in 2020: {total_2020:.2f} CNY/tonne")
    print(f"Total cost in 2050_5p: {total_2050_5p:.2f} CNY/tonne")
    print(f"Total cost in 2050_20p: {total_2050_20p:.2f} CNY/tonne")
    print(f"Total cost in 2050_100p: {total_2050_100p:.2f} CNY/tonne")
    
    # Create a single subplot for the stacked bar chart
    fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    
    # Define colors chosen to loosely reflect the meaning of each category
    colors = [
        '#8B4513',  # Raw materials – brown (mineral/raw material)
        '#2E8B57',  # Labor – dark green (people/labor)
        '#4169E1',  # Fixed O&M – blue (stable/technical)
        '#FF4500',  # Restart – orange-red (restart/dynamics)
        '#696969',  # Depreciation – gray (time/aging)
        '#DC143C',  # Retirement loss – crimson (loss/shutdown)
        '#9370DB',  # Storage – purple (storage/space)
        '#FFD700'   # Electricity – gold (power/energy)
    ]
    
    # Prepare data – each scenario corresponds to one stacked bar
    scenarios = ['2020', 
                '2050\n5% overcapacity', 
                '2050\n36% overcapacity', 
                '2050\nNo-discommissioning']
    scenario_costs = [costs_2020, costs_2050_5p, costs_2050_20p, costs_2050_100p]
    scenario_totals = [total_2020, total_2050_5p, total_2050_20p, total_2050_100p]
    
    # Set x-axis positions
    x = np.arange(len(scenarios))
    width = 0.6
    
    # Build stacked bars for each cost component
    bottom = np.zeros(len(scenarios))
    
    for i, (category, color) in enumerate(zip(categories, colors)):
        # Extract the cost of the current category in each scenario
        category_costs = [costs[i] for costs in scenario_costs]
        
        # Plot stacked bars for this category
        bars = ax.bar(x, category_costs, width, bottom=bottom, 
                     label=category, color=color, alpha=0.8, 
                     edgecolor='black', linewidth=0.5)
        
        # Add labels for mid-sized components to keep the figure readable
        for j, (bar, cost) in enumerate(zip(bars, category_costs)):
            if cost > 1 and cost < 8000:  # Only show labels for intermediate-sized components
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2, 
                           bar.get_y() + height/2, 
                           f'{cost:.0f}', ha='center', va='center', 
                           fontsize=16, fontweight='bold', color='black')
        
        # Update stacked bar baseline for the next component
        bottom += category_costs
    
    # Axis labels and title
    # ax.set_xlabel('Scenarios', fontsize=20, fontweight='bold')
    ax.set_ylabel('Levelized cost (CNY/tonne)', fontsize=20, fontweight='bold')
    # ax.set_title('Aluminum Cost Composition Comparison\n(Stacked Bar Chart)', 
    #             fontsize=16, fontweight='bold', pad=20)
    
    # X-axis labels
    ax.set_xticks(x)
    ax.set_xticklabels(scenarios, fontsize=20)
    
    # Y-axis tick labels
    ax.set_yticks(np.arange(8000, 20000, 2000))
    ax.set_yticklabels([f'{int(tick)}' for tick in np.arange(8000, 20000, 2000)], fontsize=20)
    
    # Add total-cost labels at the top of each bar
    for i, total in enumerate(scenario_totals):
        ax.text(i, total + 200, f'Total: {total:.0f}', 
               ha='center', va='bottom', fontsize=20, fontweight='bold')
    
    # Legend in the middle, two rows, reversed order so electricity appears on top
    handles, labels = ax.get_legend_handles_labels()
    # Reverse order so the legend matches the visual stacking (top component listed first)
    handles = handles[::-1]
    labels = labels[::-1]
    ax.legend(handles, labels, bbox_to_anchor=(0.5, -0.1), loc='upper center', 
              ncol=4, fontsize=20, frameon=False)
    
    # Light grid on y-axis
    ax.grid(True, alpha=0.3, axis='y')
    
    # Fix y-axis range for consistent comparison across cases
    ax.set_ylim(8000, 20000)
    
    # Tighten layout
    plt.tight_layout()
    
    # Save figure
    output_path = 'results/aluminum_cost_composition_2020_2050_stacked_bar.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Image saved to: {output_path}")
    
    # Show figure
    plt.show()
    
    return df

def create_detailed_comparison_table():
    """Print a detailed comparison table of the cost breakdown."""
    
    # Data preparation
    categories = ['Raw Materials', 'Labor', 'Fixed o&m', 'Restart', 'Depreciation', 'Retirement loss', 'Storage', 'Electricity']
    costs_2020 = [8451.2, 150, 1900, 0, 300, 0, 1, 6250]
    costs_2050_5p = [8451.2, 159.889114, 2171.893745, 769.1493874, 347.0, 809.8, 2.4, 4529.662849]
    costs_2050_20p = [8451.2, 159.1901429, 2984.318766, 299.1122492, 462.1, 692.7, 15.6, 2476.284987]
    costs_2050_100p = [8451.2, 160.7041192, 7326.478149, 707.9584908, 1156.812339, 0, 25.2, 797.0111764]
    
    # Create DataFrame with cost breakdown per tonne
    df = pd.DataFrame({
        'Category': categories,
        '2020 cost (CNY/tonne)': costs_2020,
        '2050_5p cost (CNY/tonne)': costs_2050_5p,
        '2050_20p cost (CNY/tonne)': costs_2050_20p,
        '2050_100p cost (CNY/tonne)': costs_2050_100p
    })
    
    # Compute total levelized cost per tonne for each case
    total_2020 = df['2020 cost (CNY/tonne)'].sum()
    total_2050_5p = df['2050_5p cost (CNY/tonne)'].sum()
    total_2050_20p = df['2050_20p cost (CNY/tonne)'].sum()
    total_2050_100p = df['2050_100p cost (CNY/tonne)'].sum()
    
    print("\n=== Detailed comparison of aluminum cost composition ===")
    print(df.to_string(index=False, float_format='%.2f'))
    print("\nTotal cost comparison (CNY/tonne):")
    print(f"2020:       {total_2020:.2f}")
    print(f"2050_5p:    {total_2050_5p:.2f} (Δ: {total_2050_5p - total_2020:.2f})")
    print(f"2050_20p:   {total_2050_20p:.2f} (Δ: {total_2050_20p - total_2020:.2f})")
    print(f"2050_100p:  {total_2050_100p:.2f} (Δ: {total_2050_100p - total_2020:.2f})")
    
    return df

def main():
    """CLI entry point for cost-composition visualization."""
    print("Creating aluminum cost composition bar chart...")
    
    # Create the bar chart
    df = create_aluminum_cost_bar_chart()
    
    # Create the detailed comparison table
    create_detailed_comparison_table()
    
    print("\nVisualization completed!")

if __name__ == "__main__":
    main()
