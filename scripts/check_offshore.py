import geopandas as gpd
import atlite
import xarray as xr
import numpy as np
import pandas as pd

print("1. 检查 offshore_province_shapes:")
EEZ_province_shp = gpd.read_file("data/resources/regions_offshore_province.geojson")
EEZ_province_shp = EEZ_province_shp.set_index('province')
print("列名:", EEZ_province_shp.columns.tolist())
print("索引:", EEZ_province_shp.index.tolist())

print("\n2. 检查 provinces_shp:")
provinces_shp = gpd.read_file("data/province_shapes/CHN_adm1.shp")[['NAME_1', 'geometry']]
provinces_shp.replace(to_replace={'Nei Mongol': 'InnerMongolia',
                              'Xinjiang Uygur': 'Xinjiang',
                              'Ningxia Hui': 'Ningxia',
                              'Xizang': 'Tibet'}, inplace=True)
print("替换后的省份名称:", provinces_shp['NAME_1'].tolist())

print("\n3. 检查 offwind_pro_names:")
offwind_pro_names = np.array(['Fujian', 'Guangdong', 'Guangxi', 'Hainan', 'Hebei',
                           'Jiangsu', 'Liaoning', 'Shandong', 'Shanghai', 'Tianjin', 'Zhejiang'],
                           dtype=str)
print("沿海省份:", offwind_pro_names.tolist())

print("\n4. 检查 cutout:")
try:
    cutout = atlite.Cutout("cutouts/China-2020.nc")
    cutout.prepare()
    print("Cutout 加载成功")
    print("时间范围:", cutout.coords['time'].values[[0, -1]])
except Exception as e:
    print("Cutout 加载失败:", str(e)) 