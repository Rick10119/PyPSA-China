#!/usr/bin/env python3
"""
从现有的all_plot_data_costs.csv文件读取数据并生成可视化图表
只保留画图功能，不改变CSV文件
"""

import logging
import os
import sys
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def load_plot_data(csv_path):
    """
    从CSV文件加载绘图数据
    
    Parameters:
    -----------
    csv_path : str or Path
        CSV文件路径
        
    Returns:
    --------
    pd.DataFrame
        加载的数据
    """
    try:
        df = pd.read_csv(csv_path)
        logger.info(f"成功加载数据，共 {len(df)} 行")
        logger.info(f"数据列: {list(df.columns)}")
        logger.info(f"唯一场景数: {df['Scenario_Code'].nunique()}")
        logger.info(f"唯一分类数: {df['Category'].nunique()}")
        return df
    except Exception as e:
        logger.error(f"加载CSV文件 {csv_path} 时出错: {str(e)}")
        return None

def generate_scenario_plots_from_csv(df, output_dir):
    """
    从CSV数据生成场景对比图表
    
    Parameters:
    -----------
    df : pd.DataFrame
        绘图数据
    output_dir : Path
        输出目录
    """
    # 创建输出目录
    plots_dir = output_dir / "scenario_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"图表输出目录: {plots_dir}")
    
    # 定义场景代码到描述的映射
    scenario_descriptions = {
        'L': 'Low',
        'M': 'Mid', 
        'H': 'High',
        'N': 'Non-constrained'
    }
    
    # 定义需求级别和市场级别
    demand_levels = ['L', 'M', 'H']
    market_levels = ['L', 'M', 'H']
    flexibility_levels = ['L', 'M', 'H', 'N']
    
    # 硬编码成本分类颜色
    category_colors = {
        # variable cost-non-renewable
        "Non-renewable operation": "#ff8c00",  # coal color
        
        # capital-non-renewable
        "Non-renewable investment": "#545454",  # nuclear color
        
        # Heating electrification
        "Heating electrification": "#bf13a0",  # heat pump color
        
        # capital-renewable
        "Renewable investment": "#2fb537",  # Renewable investment color

        # transmission lines
        "Transmission lines": "#6c9459",
        
        # batteries
        "Batteries": "#ace37f",  # battery color
        
        # long-duration storages
        "Long-duration storages": "#235ebc",
        
        # carbon capture
        "carbon capture": "#f29dae",  # CO2 color
        
        # synthetic fuels
        "Synthetic fuels": "#9850ad",  # Sabatier color
    }
    
    # 定义资源分类的优先级顺序，用于在正负号相同时进行排序
    category_priority = {
        "Renewable investment": 1,
        "Non-renewable operation": 5,
        "Non-renewable investment": 6,
        "Transmission lines": 7,
        "Batteries": 8,
        "Long-duration storages": 9,
        "Heating electrification": 10,
    }
    
    # 按照优先级对所有分类进行统一排序
    def sort_key(category):
        priority = category_priority.get(category, 999)
        return priority
    
    # 获取所有分类并排序
    all_categories = sorted(df['Category'].unique(), key=sort_key)
    
    # 为所有分类统一分配颜色
    global_category_colors = {}
    for cat_idx, category in enumerate(all_categories):
        if category in category_colors:
            color = category_colors[category]
        else:
            # 使用默认颜色映射
            color = plt.cm.tab20(cat_idx % 20)
        global_category_colors[category] = color
    
    # 为每个market-flexibility组合创建子图
    # 布局：3行（Market: L,M,H）x 4列（Flexibility: L,M,H,N）
    fig, axes = plt.subplots(3, 4, figsize=(24, 16), sharey=True)
    
    # 使用gridspec来精确控制行间距
    from matplotlib import gridspec
    gs = gridspec.GridSpec(3, 4, figure=fig, height_ratios=[1, 1, 1], hspace=0.1, wspace=0.1)
    
    # 重新分配axes
    axes = []
    for i in range(3):
        row_axes = []
        for j in range(4):
            ax = fig.add_subplot(gs[i, j])
            row_axes.append(ax)
        axes.append(row_axes)
    axes = np.array(axes)
    
    # 修复：正确的循环逻辑
    # i对应Market (0=L, 1=M, 2=H)
    # j对应Flexibility (0=L, 1=M, 2=H, 3=N)
    for i, market in enumerate(market_levels):
        for j, flex in enumerate(flexibility_levels):
            ax = axes[i, j]  # 正确的索引：axes[market_index, flexibility_index]
            
            # 从CSV数据中筛选当前market-flexibility组合的数据
            current_data = df[
                (df['Market'] == market) & 
                (df['Flexibility'] == flex)
            ]
            
            if not current_data.empty:
                # 获取该market-flexibility组合下所有demand级别的数据
                demand_data_dict = {}
                for demand in demand_levels:
                    demand_data = current_data[current_data['Demand'] == demand]
                    if not demand_data.empty:
                        demand_data_dict[demand] = {}
                        for _, row in demand_data.iterrows():
                            demand_data_dict[demand][row['Category']] = row['Value (CNY)']
                
                if demand_data_dict:
                    # 获取所有分类
                    all_categories = set()
                    for demand_data in demand_data_dict.values():
                        all_categories.update(demand_data.keys())
                    
                    if all_categories:
                        categories = list(all_categories)
                        demand_names = list(demand_data_dict.keys())
                        
                        # 创建堆叠柱状图
                        x_pos = np.arange(len(demand_names))
                        width = 0.6  # 减少柱子宽度，让图表更美观
                        
                        # 分析每个分类在不同demand级别下的正负号情况
                        category_signs = {}
                        for category in categories:
                            signs = []
                            for demand in demand_names:
                                demand_data = demand_data_dict.get(demand, {})
                                value = demand_data.get(category, 0)
                                if value < 0:  # 成本减少
                                    signs.append(-1)
                                elif value > 0:  # 成本增加
                                    signs.append(1)
                                else:  # 无变化
                                    signs.append(0)
                            category_signs[category] = signs
                        
                        # 按照正负号一致性进行排序：先按主要变化方向排序，再按优先级排序
                        def sort_key(category):
                            signs = category_signs[category]
                            # 计算主要变化方向（-1表示主要减少，1表示主要增加，0表示无变化）
                            if all(s == 0 for s in signs):
                                main_direction = 0
                            else:
                                # 计算加权平均方向
                                total_change = sum(abs(s) for s in signs)
                                if total_change == 0:
                                    main_direction = 0
                                else:
                                    weighted_direction = sum(s * abs(s) for s in signs) / total_change
                                    main_direction = -1 if weighted_direction < -0.1 else (1 if weighted_direction > 0.1 else 0)
                            
                            # 优先级（数字越小优先级越高）
                            priority = category_priority.get(category, 999)
                            
                            # 返回排序键：(主要方向, 优先级)
                            # 主要方向：-1(减少) -> 0(无变化) -> 1(增加)
                            return (main_direction, priority)
                        
                        # 按照新的排序逻辑重新排列分类
                        categories = sorted(categories, key=sort_key)
                        
                        # 使用预定义的全局颜色
                        colors = []
                        for category in categories:
                            color = global_category_colors.get(category, plt.cm.tab20(len(colors) % 20))
                            colors.append(color)
                        
                        # 准备堆叠数据 - 为每个demand级别创建完整的数组
                        positive_changes = []
                        negative_changes = []
                        
                        for demand in demand_names:
                            demand_data = demand_data_dict.get(demand, {})
                            demand_positive = []
                            demand_negative = []
                            
                            for category in categories:
                                value = demand_data.get(category, 0)
                                # 调整方向：成本减少（负值）显示在上方，成本增加（正值）显示在下方
                                if value < 0:  # 成本减少，显示在上方
                                    demand_positive.append(abs(value))  # 取绝对值
                                    demand_negative.append(0)
                                elif value > 0:  # 成本增加，显示在下方
                                    demand_positive.append(0)
                                    demand_negative.append(value)
                                else:
                                    demand_positive.append(0)
                                    demand_negative.append(0)
                            
                            positive_changes.append(demand_positive)
                            negative_changes.append(demand_negative)
                        
                        # 用于跟踪已经添加到legend的分类
                        added_to_legend = set()
                        
                        # 绘制正值堆叠（成本减少，在横轴上面）
                        bottom_positive = np.zeros(len(x_pos))
                        for cat_idx, category in enumerate(categories):
                            category_values = [changes[cat_idx] for changes in positive_changes]
                            if any(val > 0 for val in category_values):
                                # 只在第一次遇到该分类时添加到legend
                                label = category if category not in added_to_legend else ""
                                ax.bar(x_pos, category_values, width, 
                                       bottom=bottom_positive,
                                       color=colors[cat_idx], 
                                       alpha=0.8,
                                       label=label)
                                bottom_positive += np.array(category_values)
                                added_to_legend.add(category)
                        
                        # 绘制负值堆叠（成本增加，在横轴下面）
                        bottom_negative = np.zeros(len(x_pos))
                        for cat_idx, category in enumerate(categories):
                            category_values = [changes[cat_idx] for changes in negative_changes]
                            if any(val > 0 for val in category_values):
                                # 负值部分也添加到legend，但只在第一次遇到该分类时添加
                                label = category if category not in added_to_legend else ""
                                ax.bar(x_pos, -np.array(category_values), width,  # 使用负值
                                       bottom=bottom_negative,
                                       color=colors[cat_idx], 
                                       alpha=0.8,
                                       label=label)
                                bottom_negative += np.array(category_values)
                                added_to_legend.add(category)
                        
                        # 计算并显示净值（每个柱子的总高度）
                        net_values = []
                        for i, demand in enumerate(demand_names):
                            demand_data = demand_data_dict.get(demand, {})
                            total_value = sum(demand_data.values())
                            net_values.append(total_value)
                        
                        # 在堆叠图顶部显示净值
                        for i, (x, net_value) in enumerate(zip(x_pos, net_values)):
                            # 计算文本位置：正值堆叠的顶部
                            text_y = bottom_positive[i] + 5e9  # 在堆叠图顶部稍微上方
                            
                            # 格式化净值显示
                            value_text = f'{-net_value/1e9:.0f}'
                            
                            # 添加净值文本
                            ax.text(x, text_y, value_text, 
                                   ha='center', va='bottom' if net_value >= 0 else 'top',
                                   fontsize=10, fontweight='bold')
                        
                        # 设置y轴刻度为10为间隔，统一范围
                        y_min, y_max = -20e9, 100e9
                        y_ticks = np.arange(y_min, y_max + 10e9, 10e9)
                        y_tick_labels = [f'{int(tick/1e9)}' for tick in y_ticks]
                        
                        # 强制设置y轴范围和刻度，避免matplotlib自动调整
                        ax.set_ylim(y_min, y_max)
                        ax.set_yticks(y_ticks)
                        ax.set_yticklabels(y_tick_labels, fontsize=12)
                        
                        # 完全禁用自动刻度调整，确保只显示指定的刻度
                        ax.yaxis.set_major_locator(plt.FixedLocator(y_ticks))
                        ax.yaxis.set_minor_locator(plt.NullLocator())
                        
                        # 设置x轴标签为Demand: L, M, H (只有第一个显示Demand)
                        ax.set_xticks(x_pos)
                        x_labels = []
                        for i, demand in enumerate(demand_names):
                            if i == 0:
                                x_labels.append(f'Demand: {demand}')
                            else:
                                x_labels.append(f'{demand}')
                        ax.set_xticklabels(x_labels, fontsize=14)
                        
                        # 禁用所有自动刻度调整
                        ax.tick_params(axis='both', which='both', left=False, right=False, top=False, bottom=False)
                        ax.tick_params(axis='x', which='minor', bottom=False)
                        ax.tick_params(axis='y', which='minor', left=False)
                        
                        # 确保y轴只显示主要刻度，不显示次要刻度
                        ax.tick_params(axis='y', which='major', left=True)
                        ax.tick_params(axis='y', which='minor', left=False)
                        
                        # 只在最左边的子图显示y轴标签
                        if j == 0:
                            ax.set_ylabel('Cost Change (Billion CNY)', fontsize=14)
                        else:
                            ax.set_ylabel('')
                        
                        # 添加零线
                        ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1.5)
                        
                        # 设置标题和标签
                        # 每列顶部显示Flexibility标签
                        if i == 0:  # 第一行显示flexibility标签
                            ax.set_title(f'Flexibility: {scenario_descriptions[flex]}', 
                                       fontsize=14, fontweight='bold', pad=10)
                        # 每行左侧显示Market标签
                        if j == 0:  # 第一列显示market标签
                            ax.text(-0.2, 0.5, f'Market: {scenario_descriptions[market]}', 
                                   fontsize=14, fontweight='bold', rotation=90, 
                                   ha='center', va='center', transform=ax.transAxes)
                        
                        # 第三列和第四列的y轴标签向左移动
                        if j >= 2:  # 第三列和第四列
                            ax.tick_params(axis='y', pad=12.5)  # 增加标签与轴的距离
                    else:
                        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', 
                               transform=ax.transAxes, fontsize=10)
                else:
                    ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', 
                           transform=ax.transAxes, fontsize=10)
            else:
                ax.text(0.5, 0.5, f'No data for\nMarket:{market}, Flexibility:{flex}', ha='center', va='center', 
                       transform=ax.transAxes, fontsize=10)
    
    # 创建图例
    if not df.empty and global_category_colors:
        # 按照预定义的顺序创建图例，确保与图表中的颜色一致
        legend_elements = []
        for category in all_categories:
            if category in global_category_colors:
                color = global_category_colors[category]
                legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=color, 
                                                   label=category, alpha=0.8))
        
        # 在图表下方添加统一图例，横向排列
        fig.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.05),
                   ncol=min(len(legend_elements), 5), 
                   fontsize=15)
        
        logger.info(f"图例包含 {len(legend_elements)} 个分类")
    
    plt.tight_layout()
    
    # 保存图表
    plot_file = plots_dir / "scenario_comparison_costs.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight', pad_inches=0.3)
    logger.info(f"Scenario comparison plot saved to: {plot_file}")
    
    plt.show()
    
    plt.close()

def generate_summary_table_from_csv(df, output_dir):
    """
    从CSV数据生成场景对比摘要表格
    
    Parameters:
    -----------
    df : pd.DataFrame
        绘图数据
    output_dir : Path
        输出目录
    """
    # 创建输出目录
    tables_dir = output_dir / "summary_tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    
    # 准备表格数据
    table_data = []
    
    for scenario_code in df['Scenario_Code'].unique():
        scenario_data = df[df['Scenario_Code'] == scenario_code]
        
        if not scenario_data.empty:
            # 获取场景信息
            flexibility = scenario_data['Flexibility'].iloc[0]
            demand = scenario_data['Demand'].iloc[0]
            market = scenario_data['Market'].iloc[0]
            
            # 计算各分类的变化量
            category_changes = {}
            for _, row in scenario_data.iterrows():
                category = row['Category']
                value = row['Value (Billion CNY)']  # 已经是十亿人民币单位
                category_changes[category] = value
            
            # 找出变化最大的前3个分类
            top_changes = sorted(category_changes.items(), key=lambda x: abs(x[1]), reverse=True)[:3]
            
            row = {
                'Scenario': scenario_code,
                'Flexibility': flexibility,
                'Demand': demand,
                'Market': market,
                'Top Category 1': f"{top_changes[0][0]}: {top_changes[0][1]:.2f}B" if len(top_changes) > 0 else "N/A",
                'Top Category 2': f"{top_changes[1][0]}: {top_changes[1][1]:.2f}B" if len(top_changes) > 1 else "N/A",
                'Top Category 3': f"{top_changes[2][0]}: {top_changes[2][1]:.2f}B" if len(top_changes) > 2 else "N/A"
            }
            table_data.append(row)
    
    if table_data:
        # 创建DataFrame并排序
        df_table = pd.DataFrame(table_data)
        df_table = df_table.sort_values(['Demand', 'Market'])
        
        # 保存表格
        table_file = tables_dir / "scenario_summary_costs.csv"
        df_table.to_csv(table_file, index=False)
        logger.info(f"Summary table saved to: {table_file}")
        
        # 打印表格
        print(f"\n=== Costs Summary Table ===")
        print(df_table.to_string(index=False))
    else:
        logger.warning("No data available for summary table")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='从CSV文件读取数据并生成可视化图表')
    parser.add_argument('--csv-file', 
                       default='results/scenario_analysis/scenario_plots/all_plot_data_costs.csv',
                       help='CSV文件路径')
    parser.add_argument('--output', default='results/scenario_analysis', help='输出目录')
    
    args = parser.parse_args()
    
    # 设置日志级别
    log_level = logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    
    logger.info(f"开始从CSV文件生成图表")
    logger.info(f"CSV文件: {args.csv_file}")
    logger.info(f"输出目录: {args.output}")
    
    # 检查CSV文件是否存在
    if not Path(args.csv_file).exists():
        logger.error(f"CSV文件不存在: {args.csv_file}")
        return
    
    # 加载数据
    df = load_plot_data(args.csv_file)
    if df is None:
        logger.error("无法加载数据，退出")
        return
    
    # 创建输出目录
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 生成场景对比图表
    logger.info("生成场景对比图表...")
    generate_scenario_plots_from_csv(df, output_path)
    
    # 生成摘要表格
    logger.info("生成摘要表格...")
    generate_summary_table_from_csv(df, output_path)
    
    logger.info("分析完成！")
    logger.info(f"结果保存在: {output_path}")

if __name__ == "__main__":
    main()
