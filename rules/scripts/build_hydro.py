# 示例代码生成 hydro_p_nom.h5
import pandas as pd

# 创建数据
hydro_p_nom = pd.Series({
    'Beijing': 0,
    'Tianjin': 0,
    'Hebei': 1000,  # MW
    'Sichuan': 80000,
    'Yunnan': 70000,
    # ... 其他省份
}, name='p_nom')

# 保存为h5文件
hydro_p_nom.to_hdf('data/p_nom/hydro_p_nom.h5', key='hydro_p_nom')

# 示例代码生成 hydro_p_max_pu.h5
import pandas as pd
import numpy as np

# 创建时间索引 (8760小时)
index = pd.date_range('2023-01-01', '2023-12-31 23:00', freq='H')

# 创建各省份的出力曲线 (0-1之间的值)
provinces = ['Beijing', 'Tianjin', 'Hebei', 'Sichuan', 'Yunnan', # ... 其他省份]
data = {}
for province in provinces:
    # 这里用随机数模拟，实际应该用真实的水电出力数据
    # 可以考虑季节性变化，丰水期出力高，枯水期出力低
    seasonal = np.sin(np.linspace(0, 2*np.pi, len(index))) * 0.3 + 0.5
    noise = np.random.normal(0, 0.1, len(index))
    data[province] = np.clip(seasonal + noise, 0, 1)

hydro_p_max_pu = pd.DataFrame(data, index=index)

# 保存为h5文件
hydro_p_max_pu.to_hdf('data/p_nom/hydro_p_max_pu.h5', key='hydro_p_max_pu')