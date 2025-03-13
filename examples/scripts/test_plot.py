import matplotlib.pyplot as plt
import numpy as np

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# 数据
years = ['2025', '2030', '2035', '2040']
new_storage_no_al = [54, 63, 83, 120]
hydrogen_no_al = [51, 54, 54, 53]
system_cost_no_al = [91, 81, 75, 70]

new_storage_with_al = [39, 48, 54, 96]
hydrogen_with_al = [46, 48, 50, 49]
system_cost_with_al = [82, 73, 68, 63]

# 设置柱状图的位置
x = np.arange(len(years))  # 年份的数量
width = 0.25  # 柱的宽度

# 创建图形和轴
fig, ax1 = plt.subplots(figsize=(12, 6))

# 绘制不考虑电解铝响应的柱状图
ax1.bar(x - width, new_storage_no_al, width, label='新型储能装机 (GW) - 不考虑电解铝响应', color='lightblue')
ax1.bar(x, hydrogen_no_al, width, label='电制氢装机 (GW) - 不考虑电解铝响应', color='lightgreen')
ax1.bar(x + width, system_cost_no_al, width, label='系统成本 (B €/a) - 不考虑电解铝响应', color='salmon')

# 绘制考虑电解铝响应的柱状图
ax1.bar(x - width, new_storage_with_al, width, label='新型储能装机 (GW) - 考虑电解铝响应', color='blue', alpha=0.5)
ax1.bar(x, hydrogen_with_al, width, label='电制氢装机 (GW) - 考虑电解铝响应', color='green', alpha=0.5)
ax1.bar(x + width, system_cost_with_al, width, label='系统成本 (B €/a) - 考虑电解铝响应', color='red', alpha=0.5)

# 计算减少比例并标注
for i in range(len(years)):
    storage_reduction = (new_storage_no_al[i] - new_storage_with_al[i]) / new_storage_no_al[i] * 100
    hydrogen_reduction = (hydrogen_no_al[i] - hydrogen_with_al[i]) / hydrogen_no_al[i] * 100
    system_cost_reduction = (system_cost_no_al[i] - system_cost_with_al[i]) / system_cost_no_al[i] * 100
    ax1.text(x[i] - width, new_storage_with_al[i] + 2, f'{storage_reduction:.1f}%', ha='center', color='black', fontsize=12)
    ax1.text(x[i], hydrogen_with_al[i] + 2, f'{hydrogen_reduction:.1f}%', ha='center', color='black', fontsize=12)
    ax1.text(x[i] + width, system_cost_with_al[i] + 2, f'{system_cost_reduction:.1f}%', ha='center', color='black', fontsize=12)

# 添加标签和标题
ax1.set_xlabel('价格场景', fontsize=16)
ax1.set_ylabel('容量 (GW) / 成本 (B €/a)', fontsize=16)
ax1.set_title('不同场景下的新型储能装机、电制氢装机和系统成本对比', fontsize=16)
ax1.set_xticks(x)
ax1.set_xticklabels(years, fontsize=12)

# 调整图例位置
ax1.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)

# 显示图形
plt.tight_layout()
plt.show()
