# test_plot_module.py
# 测试绘图模块的功能

import sys
import os

# 添加scripts目录到Python路径
sys.path.append(os.path.join(os.path.dirname(__file__)))

from plot_capacity_expansion import (
    plot_results, 
    analyze_ramp_constraints, 
    plot_capacity_comparison,
    plot_network_summary,
    plot_time_series_analysis,
    create_summary_plots
)

def test_plot_module():
    """测试绘图模块的导入和基本功能"""
    print("Testing plot module import...")
    
    # 检查所有函数是否都能正常导入
    functions = [
        plot_results,
        analyze_ramp_constraints,
        plot_capacity_comparison,
        plot_network_summary,
        plot_time_series_analysis,
        create_summary_plots
    ]
    
    print(f"Successfully imported {len(functions)} plotting functions:")
    for func in functions:
        print(f"  - {func.__name__}: {func.__doc__}")
    
    print("\nPlot module test complete!")
    print("All plotting functions have been successfully separated into the plot_capacity_expansion.py module.")

if __name__ == "__main__":
    test_plot_module() 