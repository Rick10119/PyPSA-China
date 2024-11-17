import geopandas as gpd
import os

# 检查所有相关的geojson文件
files = [
    "data/resources/regions_offshore.geojson",
    "data/resources/regions_offshore_province.geojson",
]

for file in files:
    if os.path.exists(file):
        print(f"\n检查文件: {file}")
        print(f"文件大小: {os.path.getsize(file) / (1024*1024):.2f} MB")
        
        # 读取并检查内容
        gdf = gpd.read_file(file)
        print(f"特征数量: {len(gdf)}")
        if 'province' in gdf.columns:
            print("省份分布:")
            print(gdf['province'].value_counts())