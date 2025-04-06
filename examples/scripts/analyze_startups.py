import pandas as pd
import matplotlib.pyplot as plt

def count_link_startups(n, link_name='smelter'):
    """统计链接组件的启动次数"""
    if hasattr(n, 'links_t') and hasattr(n.links_t, 'status'):
        # 获取状态时间序列
        status = n.links_t.status[link_name]
        
        # 计算启动事件：当状态从0变为1时
        startups = (status.diff() > 0).astype(int)
        
        # 总启动次数
        total_startups = startups.sum()
        
        # 每月启动次数
        monthly_startups = startups.resample('M').sum()
        
        return {
            'total_startups': total_startups,
            'monthly_startups': monthly_startups,
            'startup_series': startups
        }
    else:
        print("警告：网络中没有状态变量。可能没有启用整数优化或没有记录状态。")
        return None

def plot_link_startups(startup_data, link_name='smelter'):
    """绘制链接启动次数图表"""
    if startup_data is None:
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # 绘制启动事件时间序列
    ax1.step(startup_data['startup_series'].index, 
             startup_data['startup_series'].values, 
             where='post', 
             label='Startup Events')
    ax1.set_title(f'{link_name} Startup Events')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Startup (1=Yes, 0=No)')
    ax1.legend()
    
    # 绘制每月启动次数
    monthly = startup_data['monthly_startups']
    ax2.bar(range(len(monthly)), monthly.values)
    ax2.set_title(f'{link_name} Monthly Startup Count')
    ax2.set_xlabel('Month')
    ax2.set_ylabel('Number of Startups')
    ax2.set_xticks(range(len(monthly)))
    ax2.set_xticklabels([d.strftime('%Y-%m') for d in monthly.index])
    ax2.tick_params(axis='x', rotation=45)
    
    # plt.tight_layout()
    # plt.savefig(f'examples/results/{link_name}_startups.png', dpi=300, bbox_inches='tight')
    # plt.show()
    
    # 打印总启动次数
    print(f"\n{link_name} Total Startups: {startup_data['total_startups']}")
    print(f"\n{link_name} Monthly Startups:")
    for date, count in zip(monthly.index, monthly.values):
        print(f"{date.strftime('%Y-%m')}: {count}")

# 主函数中调用这些函数
def analyze_startups(n):
    """分析网络中链接的启动情况"""
    # 分析铝电解槽启动情况
    smelter_startups = count_link_startups(n, 'smelter')
    plot_link_startups(smelter_startups, 'smelter')
    
    # 你还可以分析其他链接组件
    # other_link_startups = count_link_startups(n, 'other_link_name')
    # plot_link_startups(other_link_startups, 'other_link_name') 