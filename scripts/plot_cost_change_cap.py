# SPDX-FileCopyrightText: : 2022 The PyPSA-China Authors
#
# SPDX-License-Identifier: MIT

"""
对比不同容量比例的结果summary的cost和capacity差异（仅2050年）
生成相对于基准版本的成本减少量图表
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
import yaml

logger = logging.getLogger(__name__)

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def load_config(config_path):
    """
    加载配置文件
    
    Parameters:
    -----------
    config_path : str or Path
        配置文件路径
        
    Returns:
    --------
    dict
        配置内容
    """
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        logger.error(f"加载配置文件 {config_path} 时出错: {str(e)}")
        return None

def get_base_version_from_config(config_path):
    """
    从配置文件中获取基准版本号
    
    Parameters:
    -----------
    config_path : str or Path
        配置文件路径
        
    Returns:
    --------
    str
        基准版本号（包含55p容量比例后缀）
    """
    config = load_config(config_path)
    if config is None:
        logger.warning(f"无法加载配置文件 {config_path}，使用默认版本号 0723.8H.5-55p")
        return "0723.8H.5-55p"
    
    version = config.get('version', '0723.8H.5')
    
    # 构建基准版本号（版本号-55p）
    if '-' in version and version.endswith('p'):
        # 如果版本号已经包含容量比例后缀，先移除它
        base_version_without_suffix = version.rsplit('-', 1)[0]
        base_version = f"{base_version_without_suffix}-55p"
        logger.info(f"从配置文件 {config_path} 读取到版本号: {version}，基准版本: {base_version}")
    else:
        # 如果版本号不包含容量比例后缀，直接添加-55p
        base_version = f"{version}-55p"
        logger.info(f"从配置文件 {config_path} 读取到版本号: {version}，基准版本: {base_version}")
    
    return base_version

def collect_all_capacity_data(capacity_ratios, file_type, base_version):
    """
    从所有容量比例目录中收集2050年数据，自适应跳过不存在的年份数据
    
    Parameters:
    -----------
    capacity_ratios : list
        容量比例列表，如 ['0723.8H.5', '60p', '70p', '80p', '90p', '100p']
    file_type : str
        文件类型
    base_version : str
        基准版本号
        
    Returns:
    --------
    dict
        按容量比例组织的数据
    """
    capacity_data = {}
    year = 2050
    
    for ratio in capacity_ratios:
        if ratio == base_version:
            # 基准版本（55p）直接使用
            version_name = ratio
        else:
            # 其他容量比例需要从基准版本中提取基础版本号，然后添加对应的容量比例
            if base_version.endswith('-55p'):
                base_version_without_suffix = base_version[:-4]  # 移除 '-55p'
                version_name = f"{base_version_without_suffix}-{ratio}"
            else:
                version_name = f"{base_version}-{ratio}"
        
        # 构建目录路径
        dir_pattern = f"postnetwork-ll-current+Neighbor-linear2050-{year}"
        version_dir = Path(f"results/version-{version_name}/summary/postnetworks/positive/{dir_pattern}")
        
        if version_dir.exists():
            # 读取costs文件
            costs_file = version_dir / f"{file_type}.csv"
            
            if costs_file.exists():
                try:
                    # 读取数据
                    df = load_single_csv_file(costs_file)
                    
                    if df is not None:
                        capacity_data[ratio] = df
                        logger.info(f"成功加载 {version_name} 的 {file_type} 数据")
                    else:
                        logger.warning(f"无法加载 {version_name} 的 {file_type} 数据，将使用空数据")
                        # 创建空的DataFrame，确保数据结构一致
                        capacity_data[ratio] = pd.DataFrame()
                        
                except Exception as e:
                    logger.warning(f"加载 {version_name} 的 {file_type} 数据时出错: {str(e)}，将使用空数据")
                    # 创建空的DataFrame，确保数据结构一致
                    capacity_data[ratio] = pd.DataFrame()
            else:
                logger.warning(f"在 {version_name} 目录中未找到 {file_type}.csv 文件，将使用空数据")
                # 创建空的DataFrame，确保数据结构一致
                capacity_data[ratio] = pd.DataFrame()
        else:
            logger.warning(f"版本 {version_name} 的目录不存在，将使用空数据")
            # 创建空的DataFrame，确保数据结构一致
            capacity_data[ratio] = pd.DataFrame()
    
    return capacity_data

def load_single_csv_file(file_path):
    """
    加载单个CSV文件
    
    Parameters:
    -----------
    file_path : Path
        CSV文件路径
        
    Returns:
    --------
    pd.DataFrame or None
        加载的数据
    """
    try:
        # 首先尝试读取文件的前几行来判断结构
        with open(file_path, 'r') as f:
            first_lines = [f.readline().strip() for _ in range(5)]
        
        # 检查是否有多级索引（通过检查前几行是否有空列）
        has_multiindex = any(',' in line and line.split(',')[1] == '' for line in first_lines)
        
        if has_multiindex:
            # 多级索引文件 - 读取时保留所有列
            df = pd.read_csv(file_path, header=None)
            
            # 处理多级索引结构
            if len(df.columns) >= 4:
                # 设置多级索引：前两列作为索引，第三列作为技术名称
                df.set_index([0, 1, 2], inplace=True)
                # 重命名最后一列为数值列名
                df.columns = [df.columns[0]]
                # 将数值列转换为数值类型
                df[df.columns[0]] = pd.to_numeric(df[df.columns[0]], errors='coerce')
            else:
                # 如果列数不足，使用默认的多级索引
                df = pd.read_csv(file_path, index_col=[0, 1])
                numeric_col = df.columns[0]
                df[numeric_col] = pd.to_numeric(df[numeric_col], errors='coerce')
        else:
            # 普通单级索引文件
            df = pd.read_csv(file_path, index_col=0)
            # 尝试将所有列转换为数值类型
            for col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        return df
        
    except Exception as e:
        logger.warning(f"加载 {file_path} 时出错: {str(e)}")
        return None

def generate_capacity_reduction_plot(capacity_data, output_dir, base_version):
    """
    生成相对于基准版本的成本减少量图表
    
    Parameters:
    -----------
    capacity_data : dict
        按容量比例组织的数据
    output_dir : Path
        输出目录
    base_version : str
        基准版本号
    """
    # 创建输出目录
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # 简化的成本分类映射 - 分为两个大类
    cost_category_mapping = {
        # 电力系统成本 - 包含所有现有的电力相关成本
        ('marginal', 'coal'): 'Power System Cost',
        ('marginal', 'coal power plant'): 'Power System Cost',
        ('marginal', 'coal cc'): 'Power System Cost',
        ('marginal', 'gas'): 'Power System Cost',
        ('marginal', 'nuclear'): 'Power System Cost',
        ('marginal', 'CHP coal'): 'Power System Cost',
        ('marginal', 'CHP gas'): 'Power System Cost',
        ('marginal', 'OCGT gas'): 'Power System Cost',
        ('marginal', 'coal boiler'): 'Power System Cost',
        ('marginal', 'gas boiler'): 'Power System Cost',
        ('capital', 'coal'): 'Power System Cost',
        ('capital', 'coal power plant'): 'Power System Cost',
        ('capital', 'coal cc'): 'Power System Cost',
        ('capital', 'gas'): 'Power System Cost',
        ('capital', 'nuclear'): 'Power System Cost',
        ('capital', 'CHP coal'): 'Power System Cost',
        ('capital', 'CHP gas'): 'Power System Cost',
        ('capital', 'OCGT gas'): 'Power System Cost',
        ('capital', 'coal boiler'): 'Power System Cost',
        ('capital', 'gas boiler'): 'Power System Cost',
        ('capital', 'heat pump'): 'Power System Cost',
        ('capital', 'resistive heater'): 'Power System Cost',
        ('capital', 'hydro_inflow'): 'Power System Cost',
        ('capital', 'hydroelectricity'): 'Power System Cost',
        ('capital', 'offwind'): 'Power System Cost',
        ('capital', 'onwind'): 'Power System Cost',
        ('capital', 'solar'): 'Power System Cost',
        ('capital', 'solar thermal'): 'Power System Cost',
        ('capital', 'biomass'): 'Power System Cost',
        ('capital', 'biogas'): 'Power System Cost',
        ('capital', 'AC'): 'Power System Cost',
        ('capital', 'stations'): 'Power System Cost',
        ('capital', 'battery'): 'Power System Cost',
        ('capital', 'battery discharger'): 'Power System Cost',
        ('marginal', 'battery'): 'Power System Cost',
        ('marginal', 'battery discharger'): 'Power System Cost',
        ('capital', 'PHS'): 'Power System Cost',
        ('capital', 'water tanks'): 'Power System Cost',
        ('capital', 'H2'): 'Power System Cost',
        ('capital', 'H2 CHP'): 'Power System Cost',
        ('marginal', 'PHS'): 'Power System Cost',
        ('marginal', 'water tanks'): 'Power System Cost',
        ('marginal', 'H2'): 'Power System Cost',
        ('marginal', 'H2 CHP'): 'Power System Cost',
        ('capital', 'CO2 capture'): 'Power System Cost',
        ('marginal', 'CO2 capture'): 'Power System Cost',
        ('capital', 'Sabatier'): 'Power System Cost',
        ('marginal', 'Sabatier'): 'Power System Cost',
        ('capital', 'CO2'): 'Power System Cost',
        ('marginal', 'CO2'): 'Power System Cost',
        
        # 电解铝运行成本 - 专门针对电解铝的成本
        ('marginal', 'aluminum'): 'Aluminum Operation Cost',
        ('capital', 'aluminum'): 'Aluminum Operation Cost',
        ('marginal', 'aluminum smelter'): 'Aluminum Operation Cost',
        ('capital', 'aluminum smelter'): 'Aluminum Operation Cost',
        ('shutdown', 'aluminum'): 'Aluminum Operation Cost',
        ('standby', 'aluminum'): 'Aluminum Operation Cost',
        ('startup', 'aluminum'): 'Aluminum Operation Cost',
    }
    
    # 欧元到人民币转换率
    EUR_TO_CNY = 7.8
    
    # 排除不需要展示的分类
    exclude_categories = {
        'synthetic fuels',
        'marginal - renewable',
        'marginal - heat pump',
        'marginal - resistive heater',
        'marginal - onwind',
        'marginal - offwind',
        'marginal - solar',
        'marginal - solar thermal',
        'marginal - biomass',
        'marginal - biogas',
        'marginal - H2',
        'marginal - H2 CHP',
    }
    
    # 检查是否有基准版本数据
    if base_version not in capacity_data:
        logger.error(f"缺少{base_version}基准数据，无法计算变化量")
        return
    
    baseline_data = capacity_data[base_version]
    
    # 计算相对于基准版本的变化量
    change_data = {}
    comparison_ratios = ['60p', '70p', '80p', '90p', '100p']
    net_changes = []
    
    # 首先收集所有容量比例中的所有索引，确保不遗漏任何成本项
    all_indices = set()
    all_indices.update(baseline_data.index)
    for ratio in comparison_ratios:
        if ratio in capacity_data:
            all_indices.update(capacity_data[ratio].index)
    
    print(f"调试信息：总共找到 {len(all_indices)} 个成本项")
    print(f"调试信息：基准数据中有 {len(baseline_data.index)} 个成本项")
    
    # 检查电解铝相关的成本项
    aluminum_indices = []
    for idx in all_indices:
        if len(idx) >= 3:
            component_type, cost_type, carrier = idx[0], idx[1], idx[2]
            if isinstance(carrier, str) and 'aluminum' in carrier.lower():
                aluminum_indices.append(idx)
    
    print(f"调试信息：找到 {len(aluminum_indices)} 个电解铝相关成本项")
    for idx in aluminum_indices:
        print(f"  电解铝成本项: {idx}")
    
    for ratio in comparison_ratios:
        if ratio not in capacity_data:
            logger.warning(f"缺少 {ratio} 数据，跳过")
            continue
            
        comparison_data = capacity_data[ratio]
        change_data[ratio] = {}
        
        # 按成本分类组织变化量
        category_changes = {}
        ratio_net_change = 0
        
        # 计算每个成本分类的变化量 - 遍历所有索引
        for idx in all_indices:
            if len(idx) >= 3:
                component_type, cost_type, carrier = idx[0], idx[1], idx[2]
                
                # 使用分类映射
                category_key = (cost_type, carrier)
                category_name = cost_category_mapping.get(category_key, 'Power System Cost')  # 默认为电力系统成本
                
                if category_name not in exclude_categories:
                    # 获取基准值（从基准版本数据中读取）
                    baseline_value = 0
                    if idx in baseline_data.index:
                        baseline_value = baseline_data.loc[idx].iloc[0]
                        if pd.isna(baseline_value):
                            baseline_value = 0
                    
                    # 获取对比值 - 如果数据为空或缺失，设为0
                    comparison_value = 0
                    if not comparison_data.empty and idx in comparison_data.index:
                        comparison_value = comparison_data.loc[idx].iloc[0]
                        if pd.isna(comparison_value):
                            comparison_value = 0
                    
                    # 计算变化量（baseline - comparison，减少为正，增加为负）
                    change = baseline_value - comparison_value
                    
                    if category_name not in category_changes:
                        category_changes[category_name] = 0
                    category_changes[category_name] += change
                    ratio_net_change += change
        
        # 转换为人民币
        for category, change in category_changes.items():
            category_changes[category] = change * EUR_TO_CNY
        
        change_data[ratio] = category_changes
        net_changes.append(ratio_net_change * EUR_TO_CNY)
    
    # 收集所有成本分类
    all_categories = set()
    for ratio_data in change_data.values():
        all_categories.update(ratio_data.keys())
    
    # 过滤掉总和为0的分类
    filtered_categories = {}
    for category in all_categories:
        total_change = sum(change_data[ratio].get(category, 0) for ratio in comparison_ratios if ratio in change_data)
        if abs(total_change) > 1e6:  # 只保留变化量绝对值大于1M CNY的分类
            filtered_categories[category] = {}
            for ratio in comparison_ratios:
                if ratio in change_data:
                    filtered_categories[category][ratio] = change_data[ratio].get(category, 0)
                else:
                    filtered_categories[category][ratio] = 0
    
    if not filtered_categories:
        logger.warning("没有找到有效的成本变化数据")
        return
    
    # 创建单个图表
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # 准备数据
    x = np.arange(len(comparison_ratios))
    x_labels = [f"{ratio.replace('p', '')}%" for ratio in comparison_ratios]
    
    # 绘制两个成本类别的柱状图
    bar_width = 0.35
    opacity = 0.8
    
    # 电力系统成本
    if 'Power System Cost' in filtered_categories:
        power_values = [filtered_categories['Power System Cost'].get(ratio, 0) for ratio in comparison_ratios]
        bars1 = ax.bar(x - bar_width/2, power_values, bar_width, 
                      color='#1f77b4', alpha=opacity, label='Power System Cost')
    
    # 电解铝运行成本
    if 'Aluminum Operation Cost' in filtered_categories:
        aluminum_values = [filtered_categories['Aluminum Operation Cost'].get(ratio, 0) for ratio in comparison_ratios]
        bars2 = ax.bar(x + bar_width/2, aluminum_values, bar_width, 
                      color='#ff7f0e', alpha=opacity, label='Aluminum Operation Cost')
    
    # 绘制净变化量的黑线
    ax.plot(x, net_changes, 'k-', linewidth=3, label='Net Reduction', marker='o', markersize=10, zorder=20)
    
    # 添加数值标签
    for i, net_change in enumerate(net_changes):
        if abs(net_change) > 1e3:  # 显示大于1k的变化
            ax.annotate(f'{net_change/1e9:.1f}B',
                        xy=(i, net_change),
                        xytext=(0, 10 if net_change > 0 else -20),
                        textcoords="offset points",
                        ha='center', va='bottom' if net_change > 0 else 'top', 
                        fontsize=10, weight='bold', color='black',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # 设置图表属性
    ax.set_xlabel('Capacity Ratio (vs 55%)')
    ax.set_ylabel('Cost Reduction (Billion CNY)')
    ax.set_title('Cost Reduction by Category for Different Capacity Ratios vs 55% (2050)')
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    
    # 设置y轴刻度标签为十亿人民币单位
    y_ticks = ax.get_yticks()
    y_tick_labels = [f'{tick/1e9:.1f}B' for tick in y_ticks]
    ax.set_yticklabels(y_tick_labels)
    
    plt.tight_layout()
    
    # 保存图表
    plot_file = plots_dir / f"cost_reduction_by_capacity_ratio_vs_{base_version}_2050.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    logger.info(f"Capacity ratio cost reduction plot saved to: {plot_file}")
    
    # 保存数据
    data_rows = []
    for category, ratios_data in filtered_categories.items():
        for ratio in comparison_ratios:
            # 如果某个容量比例的数据缺失，设为0
            cost_reduction = ratios_data.get(ratio, 0)
            data_rows.append({
                'Cost Category': category,
                'Capacity Ratio': ratio.replace('p', '%'),
                'Cost Reduction (CNY)': cost_reduction,
                'Cost Reduction (Billion CNY)': cost_reduction / 1e9
            })
    
    # 添加净变化量数据
    for i, ratio in enumerate(comparison_ratios):
        if i < len(net_changes):
            data_rows.append({
                'Cost Category': 'Net Reduction',
                'Capacity Ratio': ratio.replace('p', '%'),
                'Cost Reduction (CNY)': net_changes[i],
                'Cost Reduction (Billion CNY)': net_changes[i] / 1e9
            })
    
    data_df = pd.DataFrame(data_rows)
    data_file = plots_dir / f"cost_reduction_by_capacity_ratio_vs_{base_version}_2050.csv"
    data_df.to_csv(data_file, index=False)
    logger.info(f"Capacity ratio cost reduction data saved to: {data_file}")
    
    plt.close()
    
    # 打印调试信息
    print(f"\n=== 容量比例成本减少量分析（相对于55%，人民币）===")
    for i, ratio in enumerate(comparison_ratios):
        if i < len(net_changes):
            print(f"{ratio.replace('p', '%')} vs 55%: {net_changes[i]/1e9:.3f}B CNY")
        else:
            print(f"{ratio.replace('p', '%')} vs 55%: 数据缺失")
    
    print(f"\n=== 各分类减少量（100% vs 55%）===")
    for category in filtered_categories:
        change_100p = filtered_categories[category].get('100p', 0)
        if abs(change_100p) > 1e6:  # 只显示变化量大于1M的分类
            print(f"{category}: {change_100p/1e9:.3f}B CNY")
    
    # 添加缺失数据统计
    print(f"\n=== 数据完整性统计 ===")
    for ratio in capacity_ratios:
        if ratio in capacity_data:
            if capacity_data[ratio].empty:
                print(f"{ratio}: 数据缺失（空DataFrame）")
            else:
                print(f"{ratio}: 数据完整（{len(capacity_data[ratio].index)} 个成本项）")
        else:
            print(f"{ratio}: 数据缺失（未找到）")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='生成不同容量比例相对于基准版本的成本减少量图表（仅2050年）')
    parser.add_argument('--output', default='results/comparison_results', help='输出目录')
    parser.add_argument('--verbose', '-v', action='store_true', help='详细输出')
    parser.add_argument('--results-dir', default='results', help='结果目录路径 (默认: results)')
    parser.add_argument('--config', default='config.yaml', help='配置文件路径 (默认: config.yaml)')
    
    args = parser.parse_args()
    
    # 设置日志级别
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    
    # 从配置文件读取基准版本号
    base_version = get_base_version_from_config(args.config)
    
    # 定义容量比例（基准版本为55p，其他为对比版本）
    capacity_ratios = [base_version, '60p', '70p', '80p', '90p', '100p']
    
    logger.info(f"开始收集容量比例 {capacity_ratios} 的2050年数据")
    logger.info(f"基准版本: {base_version} (55%)")
    
    # 收集所有容量比例的数据
    capacity_data = collect_all_capacity_data(capacity_ratios, 'costs', base_version)
    
    if not capacity_data:
        logger.error("没有找到任何容量比例的数据")
        return
    
    logger.info(f"成功收集到 {len(capacity_data)} 个容量比例的数据")
    
    # 生成容量比例成本减少量图表
    generate_capacity_reduction_plot(capacity_data, Path(args.output), base_version)
    
    logger.info("容量比例成本减少量图表生成完成！")

if __name__ == "__main__":
    main() 