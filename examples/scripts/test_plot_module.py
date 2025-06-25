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
    print("测试绘图模块导入...")
    
    # 检查所有函数是否都能正常导入
    functions = [
        plot_results,
        analyze_ramp_constraints,
        plot_capacity_comparison,
        plot_network_summary,
        plot_time_series_analysis,
        create_summary_plots
    ]
    
    print(f"成功导入 {len(functions)} 个绘图函数:")
    for func in functions:
        print(f"  - {func.__name__}: {func.__doc__}")
    
    print("\n绘图模块测试完成！")
    print("所有绘图函数已成功从主文件中分离到 plot_capacity_expansion.py 模块中。")

if __name__ == "__main__":
    test_plot_module() 