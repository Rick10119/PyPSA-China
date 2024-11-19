import pandas as pd
import numpy as np
import os

# 创建目录
os.makedirs('data/existing_infrastructure', exist_ok=True)

# 省份列表
provinces = [
    "Beijing", "Tianjin", "Hebei", "Shanxi", "InnerMongolia",
    "Liaoning", "Jilin", "Heilongjiang", "Shanghai", "Jiangsu",
    "Zhejiang", "Anhui", "Fujian", "Jiangxi", "Shandong",
    "Henan", "Hubei", "Hunan", "Guangdong", "Guangxi",
    "Hainan", "Chongqing", "Sichuan", "Guizhou", "Yunnan", "Tibet",
    "Shaanxi", "Gansu", "Qinghai", "Ningxia", "Xinjiang"
]

# 初始数据（部分真实数据，部分为占位数据）
capacities = {
    'coal': {  # 煤电装机
        "Beijing": 2730, "Tianjin": 11470, "Hebei": 44140,
        "Inner Mongolia": 94300, "Shanxi": 67890
    },
    'nuclear': {  # 核电装机
        "Liaoning": 2167, "Jiangsu": 4840, "Zhejiang": 9140,
        "Fujian": 4090, "Guangdong": 13340, "Guangxi": 2170
    },
    'CHP coal': {  # 煤电热电联产（估算值）
        "Beijing": 5000, "Tianjin": 6000, "Hebei": 8000,
        "Heilongjiang": 10000, "Jilin": 8000
    },
    'CHP gas': {  # 燃气热电联产（估算值）
        "Beijing": 3000, "Shanghai": 4000, "Guangdong": 3500
    },
    'OCGT': {  # 开式燃气轮机（占位数据）
        "Beijing": 1000, "Shanghai": 1500, "Guangdong": 2000
    },
    'solar': {  # 光伏发电
        "Qinghai": 15800, "Xinjiang": 13600, "Gansu": 12300,
        "Inner Mongolia": 14200, "Ningxia": 9500
    },
    'solar thermal': {  # 太阳能集热器（占位数据）
        "Beijing": 500, "Tianjin": 400, "Hebei": 600
    },
    'onwind': {  # 陆上风电
        "Inner Mongolia": 37300, "Xinjiang": 21050, "Hebei": 15430
    },
    'offwind': {  # 海上风电（部分真实数据）
        "Jiangsu": 4700, "Guangdong": 2000, "Fujian": 1500
    },
    'coal boiler': {  # 煤锅炉（估算值）
        "Beijing": 2000, "Tianjin": 2500, "Hebei": 3000
    },
    'ground heat pump': {  # 地源热泵（占位数据）
        "Beijing": 1000, "Tianjin": 800, "Hebei": 600
    }
}

# 创建并保存所有CSV文件
for tech, data in capacities.items():
    df = pd.DataFrame(index=provinces)
    df['2020'] = pd.Series(data)
    df = df.fillna(0)  # 将缺失值填充为0
    df.index.name = 'Region'
    filename = f'data/existing_infrastructure/{tech} capacity.csv'
    df.to_csv(filename)
    print(f"Created {filename}")