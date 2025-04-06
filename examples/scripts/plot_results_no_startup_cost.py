# 画出基础的结果（不考虑启停成本）
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 创建数据框
data = pd.DataFrame({
    'Excess_Rate': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50],
    'Al_Capacity': [5.50, 5.55, 5.61, 5.67, 5.72, 5.78, 5.83, 5.89, 5.94, 6.00, 6.05, 6.11, 6.16, 6.21, 6.27, 6.32, 6.38, 6.43, 6.49, 6.54, 6.60, 7.15, 7.70, 8.25],
    'System_Cost': [90993.16, 90035.52, 89092.96, 88221.07, 87393.41, 86574.35, 85842.13, 85220.64, 84661.17, 84159.36, 83790.10, 83516.74, 83243.23, 83010.33, 82818.41, 82659.95, 82508.05, 82357.35, 82219.89, 82152.94, 82152.94, 82152.94, 82152.94, 82152.94],
    'Wind': [370.47, 380.93, 390.87, 387.47, 382.24, 379.37, 374.96, 368.93, 361.08, 352.59, 359.54, 361.75, 365.28, 357.51, 345.60, 341.44, 338.45, 337.80, 341.59, 342.69, 342.69, 342.69, 342.69, 342.69],
    'Solar': [225.87, 219.54, 214.68, 219.86, 225.31, 227.75, 227.35, 221.56, 218.08, 216.15, 198.83, 198.35, 197.00, 203.66, 210.71, 211.78, 211.53, 208.58, 204.23, 196.63, 196.63, 196.63, 196.63, 196.63],
    'Battery': [54.45, 53.49, 52.22, 58.65, 65.56, 68.34, 67.49, 59.56, 55.46, 52.08, 39.50, 39.09, 37.95, 42.94, 51.65, 51.87, 51.32, 48.78, 39.61, 39.25, 39.25, 39.25, 39.25, 39.25],
    'H2_Store': [50.74, 48.45, 47.16, 45.22, 43.12, 41.97, 41.48, 42.48, 43.42, 44.33, 44.89, 44.68, 44.40, 44.58, 45.01, 45.31, 45.59, 45.86, 46.00, 46.27, 46.27, 46.27, 46.27, 46.27]
})

# 创建图表
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

# 1. System Cost vs Excess Rate
ax1.plot(data['Excess_Rate'], data['System_Cost']/1000, 'b-', linewidth=2)
# 设置x轴文字大小
ax1.set_xlabel('Excess Rate (%)', fontsize=15)
ax1.set_ylabel('System Cost (M€/a)', fontsize=15)
ax1.set_title('System Cost vs Excess Rate', fontsize=15)
ax1.grid(True)

# 2. Generation Capacity vs Excess Rate
ax2.plot(data['Excess_Rate'], data['Wind'], 'b-', label='Wind', linewidth=2)
ax2.plot(data['Excess_Rate'], data['Solar'], 'r-', label='Solar', linewidth=2)
ax2.set_xlabel('Excess Rate (%)', fontsize=15)
ax2.set_ylabel('Capacity (GW)', fontsize=15)
ax2.set_title('Generation Capacity vs Excess Rate', fontsize=15)
ax2.legend(fontsize=15)
ax2.grid(True)
ax2.set_xlim(0, max(data['Excess_Rate']))
ax2.set_ylim(0, max(max(data['Wind']), max(data['Solar'])) * 1.1)

# 3. Storage Capacity vs Excess Rate
ax3.plot(data['Excess_Rate'], data['Battery'], 'g-', label='Battery', linewidth=2)
ax3.plot(data['Excess_Rate'], data['H2_Store'], 'r-', label='H2 Storage', linewidth=2)
ax3.set_xlabel('Excess Rate (%)', fontsize=15)
ax3.set_ylabel('Capacity (GW)', fontsize=15)
ax3.set_title('Storage Capacity vs Excess Rate', fontsize=15)
ax3.legend(fontsize=15)
ax3.grid(True)

# 4. Aluminum Smelter Capacity vs Excess Rate
ax4.plot(data['Excess_Rate'], data['Al_Capacity'], 'purple', linewidth=2)
ax4.set_xlabel('Excess Rate (%)', fontsize=15)
ax4.set_ylabel('Capacity (GW)', fontsize=15)
ax4.set_title('Aluminum Smelter Capacity vs Excess Rate', fontsize=15)
ax4.grid(True)

# 调整布局
plt.tight_layout()

# 保存图表
plt.savefig('examples/results/system_analysis.png', dpi=300, bbox_inches='tight')

# 显示图表
plt.show()

# 打印关键分析结果
print("\nKey Analysis Results:")
print(f"Minimum System Cost: {data['System_Cost'].min()/1000:.2f} M€/a at {data['Excess_Rate'][data['System_Cost'].argmin()]}% excess rate")
print(f"Maximum Wind Capacity: {data['Wind'].max():.2f} GW at {data['Excess_Rate'][data['Wind'].argmax()]}% excess rate")
print(f"Maximum Battery Capacity: {data['Battery'].max():.2f} GW at {data['Excess_Rate'][data['Battery'].argmax()]}% excess rate") 