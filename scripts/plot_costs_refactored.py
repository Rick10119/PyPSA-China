#!/usr/bin/env python3
"""
重构的成本分析图表生成脚本
从CSV文件读取数据并生成清晰的场景对比图表
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import logging

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def load_data(csv_path):
    """加载CSV数据"""
    try:
        df = pd.read_csv(csv_path)
        logger.info(f"成功加载数据，共 {len(df)} 行")
        return df
    except Exception as e:
        logger.error(f"加载CSV文件失败: {e}")
        return None

def get_scenario_descriptions():
    """获取场景描述映射"""
    return {
        'L': 'Low',
        'M': 'Mid', 
        'H': 'High',
        'N': 'Non-constrained'
    }

def get_category_colors():
    """获取分类颜色映射"""
    return {
        "Non-renewable operation": "#ff8c00",
        "Non-renewable investment": "#545454",
        "Heating electrification": "#bf13a0",
        "Renewable investment": "#2fb537",
        "Transmission lines": "#6c9459",
        "Batteries": "#ace37f",
        "Long-duration storages": "#235ebc",
    }

def create_cost_comparison_plot(df, output_dir):
    """创建成本对比图表"""
    
    # 创建输出目录
    plots_dir = output_dir / "scenario_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # 获取场景描述和颜色
    scenario_descriptions = get_scenario_descriptions()
    category_colors = get_category_colors()
    
    # 定义维度
    market_levels = ['L', 'M', 'H']
    flexibility_levels = ['L', 'M', 'H', 'N']
    demand_levels = ['L', 'M', 'H']
    
    # 创建图表：3行（Market）x 4列（Flexibility）
    fig, axes = plt.subplots(3, 4, figsize=(20, 15), sharey=True)
    # fig.suptitle('Cost Changes by Market and Flexibility Scenarios', fontsize=16, fontweight='bold')
    
    # 为每个Market-Flexibility组合创建子图
    for market_idx, market in enumerate(market_levels):
        for flex_idx, flexibility in enumerate(flexibility_levels):
            ax = axes[market_idx, flex_idx]
            
            # 筛选当前组合的数据
            current_data = df[
                (df['Market'] == market) & 
                (df['Flexibility'] == flexibility)
            ]
            
            if current_data.empty:
                ax.text(0.5, 0.5, f'No data\nMarket: {market}\nFlex: {flexibility}', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=10)
                continue
            
            # 按Demand级别组织数据
            demand_data = {}
            for demand in demand_levels:  # demand_levels = ['L', 'M', 'H']
                demand_subset = current_data[current_data['Demand'] == demand]
                if not demand_subset.empty:
                    demand_data[demand] = dict(zip(demand_subset['Category'], demand_subset['Value (CNY)']))
            
            if not demand_data:
                ax.text(0.5, 0.5, 'No demand data', ha='center', va='center', 
                       transform=ax.transAxes, fontsize=10)
                continue
            
            # 获取所有分类
            all_categories = set()
            for data in demand_data.values():
                all_categories.update(data.keys())
            all_categories = sorted(list(all_categories))
            
            # 准备堆叠柱状图数据 - 按照L, M, H的顺序
            demand_names = [d for d in demand_levels if d in demand_data]  # 保持L, M, H的顺序
            x_pos = np.arange(len(demand_names))
            width = 0.6
            
            # 分离正负值
            positive_data = []
            negative_data = []
            
            for demand in demand_names:
                pos_values = []
                neg_values = []
                for category in all_categories:
                    value = demand_data[demand].get(category, 0)
                    if value < 0:  # 成本减少，显示在上方
                        pos_values.append(abs(value))
                        neg_values.append(0)
                    elif value > 0:  # 成本增加，显示在下方
                        pos_values.append(0)
                        neg_values.append(value)
                    else:
                        pos_values.append(0)
                        neg_values.append(0)
                positive_data.append(pos_values)
                negative_data.append(neg_values)
            
            # 绘制堆叠柱状图
            bottom_pos = np.zeros(len(x_pos))
            bottom_neg = np.zeros(len(x_pos))
            
            for cat_idx, category in enumerate(all_categories):
                color = category_colors.get(category, plt.cm.tab20(cat_idx % 20))
                
                # 绘制正值（成本减少）
                pos_values = [data[cat_idx] for data in positive_data]
                if any(v > 0 for v in pos_values):
                    ax.bar(x_pos, pos_values, width, bottom=bottom_pos, 
                          color=color, alpha=0.8, label=category)
                    bottom_pos += np.array(pos_values)
                
                # 绘制负值（成本增加）
                neg_values = [data[cat_idx] for data in negative_data]
                if any(v > 0 for v in neg_values):
                    ax.bar(x_pos, -np.array(neg_values), width, bottom=bottom_neg, 
                          color=color, alpha=0.8)
                    bottom_neg += np.array(neg_values)
            
            # 计算并显示净值（每个柱子的总高度）
            net_values = []
            for i, demand in enumerate(demand_names):
                demand_data_dict = demand_data.get(demand, {})
                total_value = sum(demand_data_dict.values())
                net_values.append(total_value)
            
            # 在堆叠图顶部显示净值
            for i, (x, net_value) in enumerate(zip(x_pos, net_values)):
                # 计算文本位置：正值堆叠的顶部
                text_y = bottom_pos[i] + 4e9  # 在堆叠图顶部稍微上方
                
                # 格式化净值显示
                value_text = f'{-net_value/1e9:.0f}'
                
                # 添加净值文本
                ax.text(x, text_y, value_text, 
                       ha='center', va='bottom' if net_value >= 0 else 'top',
                       fontsize=12, fontweight='bold')
            
            # 设置x轴标签
            ax.set_xticks(x_pos)
            x_labels = [f'Demand: {d}' if i == 0 else d for i, d in enumerate(demand_names)]
            ax.set_xticklabels(x_labels, fontsize=15)
            
            # 设置y轴
            ax.set_ylim(-20e9, 100e9)
            y_ticks = np.arange(-20e9, 101e9, 20e9)
            y_labels = [f'{int(tick/1e9)}' for tick in y_ticks]
            ax.set_yticks(y_ticks)
            ax.set_yticklabels(y_labels, fontsize=15)
            
            # 添加零线
            ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1)
            
            # 设置标题和标签
            if market_idx == 0:  # 第一行显示Flexibility标签
                ax.set_title(f'Flexibility: {scenario_descriptions[flexibility]}', 
                           fontsize=15, fontweight='bold')
            
            if flex_idx == 0:  # 第一列显示Market标签
                ax.set_ylabel(f'Market: {scenario_descriptions[market]}\nCost Change (Billion CNY)', 
                            fontsize=15, fontweight='bold')
            else:
                ax.set_ylabel('')
    
    # 创建图例 - 按照预定义的颜色映射顺序
    legend_elements = []
    # 按照get_category_colors()中定义的顺序创建图例
    ordered_categories = [
        "Renewable investment",
        "Non-renewable operation",
        "Non-renewable investment", 
        "Transmission lines",
        "Batteries",
        "Long-duration storages",
        "Heating electrification"
    ]
    
    for category in ordered_categories:
        if category in all_categories:  # 只添加实际存在的分类
            color = category_colors.get(category, plt.cm.tab20(len(legend_elements) % 20))
            legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=color, 
                                               label=category, alpha=0.8))
    
    fig.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.05),
               ncol=min(len(legend_elements), 4), fontsize=15)
    
    plt.tight_layout()
    
    # 保存图表
    plot_file = plots_dir / "cost_comparison_refactored.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight', pad_inches=0.3)
    logger.info(f"图表已保存到: {plot_file}")
    
    plt.show()
    plt.close()

def create_summary_table(df, output_dir):
    """创建摘要表格"""
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
                value = row['Value (Billion CNY)']
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
        table_file = tables_dir / "cost_summary_refactored.csv"
        df_table.to_csv(table_file, index=False)
        logger.info(f"摘要表格已保存到: {table_file}")
        
        # 打印表格
        print(f"\n=== 成本摘要表格 ===")
        print(df_table.to_string(index=False))
    else:
        logger.warning("没有数据可用于摘要表格")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='重构的成本分析图表生成脚本')
    parser.add_argument('--csv-file', 
                       default='results/scenario_analysis/scenario_plots/all_plot_data_costs.csv',
                       help='CSV文件路径')
    parser.add_argument('--output', default='results/scenario_analysis', help='输出目录')
    
    args = parser.parse_args()
    
    logger.info(f"开始生成成本分析图表")
    logger.info(f"CSV文件: {args.csv_file}")
    logger.info(f"输出目录: {args.output}")
    
    # 检查CSV文件是否存在
    if not Path(args.csv_file).exists():
        logger.error(f"CSV文件不存在: {args.csv_file}")
        return
    
    # 加载数据
    df = load_data(args.csv_file)
    if df is None:
        logger.error("无法加载数据，退出")
        return
    
    # 创建输出目录
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 生成图表和表格
    logger.info("生成成本对比图表...")
    create_cost_comparison_plot(df, output_path)
    
    logger.info("生成摘要表格...")
    create_summary_table(df, output_path)
    
    logger.info("分析完成！")
    logger.info(f"结果保存在: {output_path}")

if __name__ == "__main__":
    main()
