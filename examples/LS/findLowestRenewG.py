import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 读取CSV文件
df = pd.read_csv('examples/data/time-series-lecture-2.csv', index_col=0, parse_dates=True)

# 计算总新能源出力
df['total_renewables'] = df['onwind'] + df['offwind'] + df['solar']

# 将数据按天分组
daily_sum = df['total_renewables'].resample('D').sum()

# 计算连续7天的滚动总和
rolling_7days = daily_sum.rolling(window=7).sum()

# 找出总和最小的7天窗口
min_index = rolling_7days.idxmin()
start_date = min_index - pd.Timedelta(days=6)

# 获取最低七天的日期范围
lowest_7days = daily_sum.loc[start_date:min_index]

# 设置中文字体以便正确显示中文
plt.rcParams['font.sans-serif'] = ['SimHei']  # 设置中文字体


print(f"新能源出力最低的连续七天:")
print(f"开始日期: {start_date.date()}")
print(f"结束日期: {min_index.date()}")
print(f"总出力: {rolling_7days.loc[min_index]:.2f}")

# 展示这7天的每日出力
print("\n每日出力详情:")
for date, value in lowest_7days.items():
    print(f"{date.date()}: {value:.2f}")

# 可视化这7天的出力情况
plt.figure(figsize=(12, 6))
plt.bar(lowest_7days.index.strftime('%Y-%m-%d'), lowest_7days.values)
plt.title('新能源出力最低的连续七天')
plt.xlabel('日期')
plt.ylabel('总出力')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('examples/LS/results/lowest_7days.png')
plt.show()

# 提取这7天的小时级数据，可用于后续分析
lowest_hours_data = df.loc[start_date:min_index + pd.Timedelta(days=1) - pd.Timedelta(seconds=1)]

# 保存这7天的小时级数据
lowest_hours_data.to_csv('examples/LS/results/lowest_7days_hourly.csv')

print(f"\n已将详细数据保存到 'examples/LS/results/lowest_7days_hourly.csv'")