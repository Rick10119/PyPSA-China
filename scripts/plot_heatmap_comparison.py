"""
比较 MMM 和 NMM 场景的热力图可视化脚本

该脚本基于 plot_heatmap.py 重写，用于比较两个不同配置场景的储能技术和铝冶炼厂运行模式。
主要功能：
1. 从 MMM 和 NMM 配置文件中读取参数
2. 加载对应的网络数据
3. 生成并排比较的热力图
4. 支持 H2、电池、水储能和铝冶炼厂的可视化

主要差异：
- MMM: iterative_optimization: true, smelter_flexibility: mid
- NMM: iterative_optimization: false, smelter_flexibility: non_constrained
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
    设置绘图样式
    """
    # 设置字体为Helvetica
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
    为特定储能技术创建热力图数据框
    
    Parameters:
    -----------
    n : pypsa.Network
        包含仿真结果的 PyPSA 网络对象
    tech : str
        要分析的储能技术 ('H2', 'battery', 或 'water')
    province_filter : str, optional
        如果提供，则只包含指定省份的储能设备
    
    Returns:
    --------
    tuple
        (用于热力图的数据框摘要, 基础功率值 MW)
    """
    # 获取特定技术的储能设备
    stores = n.stores_t.p.filter(like=tech)
    
    # 如果指定了省份过滤器，则应用过滤
    if province_filter:
        province_stores = stores.columns[stores.columns.str.contains(province_filter, case=False)]
        if len(province_stores) > 0:
            stores = stores[province_stores]
        else:
            print(f"警告: 在省份 {province_filter} 中未找到 {tech} 储能设备")
            return pd.DataFrame(), 0
    
    # 计算最大功率作为归一化的基础值
    base = abs(stores.sum(axis=1)).max()
    if base == 0:
        print(f"警告: {tech} 功率为零 {'省份 ' + province_filter if province_filter else '全国'}")
        return pd.DataFrame(), 0
    
    # 通过基础值归一化功率值
    df = stores.sum(axis=1) / base
    df = df.to_frame()
    df.reset_index(inplace=True)
    renames = {0: 'p_store'}
    df.rename(columns=renames, inplace=True)
    
    # 数据已经是中国本地时间，直接使用时间戳
    date = n.stores_t.p.filter(like='water').index
    df['Hour'] = date.hour
    df['Day'] = date.strftime('%m-%d')
    
    # 创建用于热力图可视化的透视表
    summary = pd.pivot_table(data=df, index='Hour', columns='Day', values='p_store')
    summary = summary.fillna(0).infer_objects(copy=False)
    return summary, base

def create_aluminum_df(n, province_filter=None):
    """
    为铝冶炼厂运行创建热力图数据框
    
    Parameters:
    -----------
    n : pypsa.Network
        包含仿真结果的 PyPSA 网络对象
    province_filter : str, optional
        如果提供，则只包含指定省份的铝冶炼厂
    
    Returns:
    --------
    tuple
        (用于热力图的数据框摘要, 基础功率值 MW)
    """
    # 获取铝冶炼厂链接（从电力总线的输入功率）
    aluminum_links = n.links_t.p0.filter(like='aluminum smelter')
    
    if aluminum_links.empty:
        print("警告: 在网络中未找到铝冶炼厂链接")
        return pd.DataFrame(), 0
    
    # 如果指定了省份过滤器，则应用过滤
    if province_filter:
        province_smelters = aluminum_links.columns[aluminum_links.columns.str.contains(province_filter, case=False)]
        if len(province_smelters) > 0:
            aluminum_links = aluminum_links[province_smelters]
        else:
            print(f"警告: 在省份 {province_filter} 中未找到铝冶炼厂")
            return pd.DataFrame(), 0
    
    # 计算最大功率作为归一化的基础值
    base = abs(aluminum_links.sum(axis=1)).max()
    
    if base == 0:
        print(f"警告: 铝冶炼厂功率为零 {'省份 ' + province_filter if province_filter else '全国'}")
        return pd.DataFrame(), 0
    
    # 通过基础值归一化功率值
    df = aluminum_links.sum(axis=1) / base
    df = df.to_frame()
    df.reset_index(inplace=True)
    renames = {0: 'p_smelter'}
    df.rename(columns=renames, inplace=True)
    
    # 数据已经是中国本地时间，直接使用时间戳
    date = aluminum_links.index
    df['Hour'] = date.hour
    df['Day'] = date.strftime('%m-%d')
    
    # 创建用于热力图可视化的透视表
    summary = pd.pivot_table(data=df, index='Hour', columns='Day', values='p_smelter')
    summary = summary.fillna(0).infer_objects(copy=False)
    return summary, base

def get_aluminum_storage_daily_average(n, province_filter=None):
    """
    获取铝储能的日平均值用于热力图叠加
    
    Parameters:
    -----------
    n : pypsa.Network
        包含仿真结果的 PyPSA 网络对象
    province_filter : str, optional
        如果提供，则只包含指定省份的铝储能
    
    Returns:
    --------
    tuple
        (日平均储能水平, 最小储能水平用于归一化)
    """
    # 获取铝储能数据
    aluminum_stores = n.stores_t.e.filter(like='aluminum storage')
    
    if aluminum_stores.empty:
        print("警告: 在网络中未找到铝储能")
        return pd.Series(), 0
    
    # 如果指定了省份过滤器，则应用过滤
    if province_filter:
        province_stores = aluminum_stores.columns[aluminum_stores.columns.str.contains(province_filter, case=False)]
        if len(province_stores) > 0:
            aluminum_stores = aluminum_stores[province_stores]
        else:
            print(f"警告: 在省份 {province_filter} 中未找到铝储能")
            return pd.Series(), 0
    
    # 计算所有储能的总额
    total_storage = aluminum_stores.sum(axis=1)
    
    if total_storage.empty:
        print(f"警告: 铝储能数据为空 {'省份 ' + province_filter if province_filter else '全国'}")
        return pd.Series(), 0
    
    # 获取最小储能水平用于归一化
    min_storage = total_storage.min()
    
    # 计算日平均储能水平
    daily_avg = total_storage.groupby(total_storage.index.strftime('%m-%d')).mean()
    
    # 减去最小值进行归一化
    daily_avg_normalized = daily_avg - min_storage
    
    return daily_avg_normalized, min_storage

def plot_comparison_heatmap(n_mmm, n_nmm, config, output_dir, tech, province_filter=None):
    """
    生成 MMM 和 NMM 场景的比较热力图
    
    Parameters:
    -----------
    n_mmm : pypsa.Network
        MMM 场景的网络对象
    n_nmm : pypsa.Network
        NMM 场景的网络对象
    config : dict
        包含绘图参数的配置字典
    output_dir : str
        保存热力图绘图的目录
    tech : str
        要绘制的技术类型
    province_filter : str, optional
        如果提供，则只包含指定省份的数据
    """
    freq = config["freq"]
    planning_horizon = "2050"
    
    # 创建省份特定子目录（如果按省份过滤）
    if province_filter:
        province_dir = os.path.join(output_dir, f"province_{province_filter}")
        os.makedirs(province_dir, exist_ok=True)
        plot_title_suffix = f" in {province_filter}"
    else:
        province_dir = output_dir
        plot_title_suffix = " (National)"
    
    # 创建上下排列的子图，每个图长宽比9:4，整体约1:1
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12))
    
    # 处理 MMM 场景
    if tech == "aluminum":
        df_mmm, base_mmm = create_aluminum_df(n_mmm, province_filter)
        daily_storage_mmm, min_storage_mmm = get_aluminum_storage_daily_average(n_mmm, province_filter)
    else:
        df_mmm, base_mmm = create_df(n_mmm, tech, province_filter)
        daily_storage_mmm, min_storage_mmm = None, None
    
    # 处理 NMM 场景
    if tech == "aluminum":
        df_nmm, base_nmm = create_aluminum_df(n_nmm, province_filter)
        daily_storage_nmm, min_storage_nmm = get_aluminum_storage_daily_average(n_nmm, province_filter)
    else:
        df_nmm, base_nmm = create_df(n_nmm, tech, province_filter)
        daily_storage_nmm, min_storage_nmm = None, None
    
    # 绘制 MMM 热力图
    if not df_mmm.empty and base_mmm > 0:
        base_mmm_display = str(int(base_mmm / 1e3))  # 转换为 GW 显示
        
        if tech == "aluminum":
            sns.heatmap(df_mmm, ax=ax1, cmap='coolwarm', cbar_kws={'label': 'pu', 'shrink': 0.8}, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(df_mmm, ax=ax1, cmap='coolwarm', cbar_kws={'label': 'pu', 'shrink': 0.8}, vmin=-1.0, vmax=1.0)
        
        ax1.set_title(f'Mid smelter flexibility')
        # 隐藏第一个子图的x轴标签（Day）
        ax1.set_xlabel('')
        # ax1.set_xticklabels([])
        # 设置y轴标签正常站立（不旋转）
        ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0)
        
        # 为铝冶炼厂添加储能水平叠加（不显示标签）
        if tech == "aluminum" and not daily_storage_mmm.empty:
            day_columns = df_mmm.columns
            storage_values = []
            storage_positions = []
            
            for i, day in enumerate(day_columns):
                if day in daily_storage_mmm.index:
                    storage_values.append(daily_storage_mmm[day]/1e6)  # 转换为万吨
                    storage_positions.append(i + 0.5)  # 列的中心
            
            if storage_values:
                ax1_twin = ax1.twinx()
                # 先画白线作为外边框
                ax1_twin.plot(storage_positions, storage_values, 'w-', linewidth=3.5, zorder=1)
                # 再画黑线在内部
                ax1_twin.plot(storage_positions, storage_values, 'k-', linewidth=2, zorder=2)
                ax1_twin.set_ylabel('Stored aluminum (Mt)', color='black')
                ax1_twin.tick_params(axis='y', labelcolor='black')
                ax1_twin.set_xlim(0, len(day_columns))
    
    # 绘制 NMM 热力图
    if not df_nmm.empty and base_nmm > 0:
        base_nmm_display = str(int(base_nmm / 1e3))  # 转换为 GW 显示
        
        if tech == "aluminum":
            sns.heatmap(df_nmm, ax=ax2, cmap='coolwarm', cbar_kws={'label': 'pu', 'shrink': 0.8}, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(df_nmm, ax=ax2, cmap='coolwarm', cbar_kws={'label': 'pu', 'shrink': 0.8}, vmin=-1.0, vmax=1.0)
        
        ax2.set_title(f'Non-constrained smelter flexibility')
        # 为第二个子图添加x轴标签
        ax2.set_xlabel('Day')
        # 设置y轴标签正常站立（不旋转）
        ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0)
        
        # 为铝冶炼厂添加储能水平叠加
        if tech == "aluminum" and not daily_storage_nmm.empty:
            day_columns = df_nmm.columns
            storage_values = []
            storage_positions = []
            
            for i, day in enumerate(day_columns):
                if day in daily_storage_nmm.index:
                    storage_values.append(daily_storage_nmm[day]/1e6)  # 转换为万吨
                    storage_positions.append(i + 0.5)  # 列的中心
            
            if storage_values:
                ax2_twin = ax2.twinx()
                # 先画白线作为外边框
                ax2_twin.plot(storage_positions, storage_values, 'w-', linewidth=3.5, zorder=1)
                # 再画黑线在内部
                ax2_twin.plot(storage_positions, storage_values, 'k-', linewidth=2, label='Stored aluminum', zorder=2)
                ax2_twin.set_ylabel('Stored aluminum (Mt)', color='black')
                ax2_twin.tick_params(axis='y', labelcolor='black')
                ax2_twin.legend(loc='lower right', bbox_to_anchor=(1.1, -0.5))
                ax2_twin.set_xlim(0, len(day_columns))
    
    # 调整布局，为图例留出更多空间
    plt.tight_layout()
    # 进一步调整子图间距，为右侧图例留出更多空间
    plt.subplots_adjust(right=2)
    
    # plt.show()
    
    # 保存图片
    output_path = os.path.join(province_dir, f"smelter_flexibility_comparison.png")
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved {tech} comparison heatmap to: {output_path}")
    plt.close()

def main():
    """
    主函数
    """
    parser = argparse.ArgumentParser(description='Compare heatmaps between MMM and NMM scenarios')
    parser.add_argument('--config-mmm', type=str, 
                       default='configs/config_MMM_2050_20p.yaml',
                       help='MMM configuration file path')
    parser.add_argument('--config-nmm', type=str,
                       default='configs/config_NMM_2050_20p.yaml', 
                       help='NMM configuration file path')
    parser.add_argument('--output-dir', type=str,
                       default='results/comparison_heatmaps',
                       help='Output directory')
    parser.add_argument('--province', type=str, default=None,
                       help='Province filter (optional)')
    parser.add_argument('--techs', nargs='+', 
                       default=['H2', 'battery', 'water', 'aluminum'],
                       help='List of technologies to compare')
    
    args = parser.parse_args()
    
    # 设置绘图样式
    set_plot_style()
    
    # 加载配置文件
    print("Loading configuration files...")
    with open(args.config_mmm, 'r', encoding='utf-8') as f:
        config_mmm = yaml.safe_load(f)
    
    with open(args.config_nmm, 'r', encoding='utf-8') as f:
        config_nmm = yaml.safe_load(f)
    
    # 获取绘图参数
    map_figsize = config_mmm["plotting"]['map']['figsize']
    
    # 构建网络文件路径
    mmm_version = config_mmm['version']
    nmm_version = config_nmm['version']
    
    mmm_network_path = f"results/version-{mmm_version}/postnetworks/positive/postnetwork-ll-current+Neighbor-linear2050-2050.nc"
    nmm_network_path = f"results/version-{nmm_version}/postnetworks/positive/postnetwork-ll-current+Neighbor-linear2050-2050.nc"
    
    # 检查网络文件是否存在
    if not os.path.exists(mmm_network_path):
        print(f"Error: MMM network file not found: {mmm_network_path}")
        return
    
    if not os.path.exists(nmm_network_path):
        print(f"Error: NMM network file not found: {nmm_network_path}")
        return
    
    # 加载网络
    print("Loading network data...")
    n_mmm = pypsa.Network(mmm_network_path)
    n_nmm = pypsa.Network(nmm_network_path)
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 生成比较热力图
    print("Generating comparison heatmaps...")
    for tech in args.techs:
        print(f"Processing {tech} technology...")
        plot_comparison_heatmap(n_mmm, n_nmm, config_mmm, args.output_dir, tech, args.province)
    
    print("All comparison heatmaps generated successfully!")

if __name__ == "__main__":
    main()
