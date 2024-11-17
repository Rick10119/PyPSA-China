import geopandas as gpd
import json
from shapely.geometry import Polygon, MultiPolygon, mapping

# 读取现有的 offshore regions
offshore = gpd.read_file(snakemake.input.offshore_shapes)

# 定义沿海省份及其大致范围
province_bounds = {
    "Liaoning": [120.0, 38.5, 124.5, 41.0],
    "Hebei": [118.0, 38.0, 120.0, 40.0],
    "Shandong": [119.0, 35.0, 122.5, 38.0],
    "Jiangsu": [120.0, 31.5, 123.0, 35.0],
    "Shanghai": [121.0, 30.5, 123.0, 31.5],
    "Zhejiang": [120.0, 27.0, 123.0, 30.5],
    "Fujian": [118.5, 23.5, 120.5, 27.0],
    "Guangdong": [110.0, 20.0, 117.0, 23.5],
    "Guangxi": [107.0, 20.0, 110.0, 22.0],
    "Hainan": [108.0, 18.0, 111.0, 20.0]
}

# 创建新的 GeoDataFrame
features = []
for province, bounds in province_bounds.items():
    # 获取该省份范围内的 offshore 区域
    minx, miny, maxx, maxy = bounds
    mask = offshore.geometry.intersects(Polygon([
        [minx, miny], [maxx, miny], 
        [maxx, maxy], [minx, maxy], [minx, miny]
    ]))
    
    if any(mask):
        province_offshore = offshore[mask].copy()
        province_offshore['province'] = province
        province_offshore['name'] = f"{province}_Offshore"
        
        # 将每个特征转换为可序列化的字典
        for _, row in province_offshore.iterrows():
            feature = {
                "type": "Feature",
                "properties": {
                    "province": row['province'],
                    "name": row['name']
                },
                "geometry": mapping(row['geometry'])  # 使用 mapping 函数转换几何对象
            }
            features.append(feature)

# 创建新的 GeoJSON
offshore_province = {
    "type": "FeatureCollection",
    "crs": {
        "type": "name",
        "properties": {
            "name": "urn:ogc:def:crs:OGC:1.3:CRS84"
        }
    },
    "features": features
}

# 保存文件
with open(snakemake.output.offshore_province, 'w') as f:
    json.dump(offshore_province, f, indent=2)

print(f"已创建文件: {snakemake.output.offshore_province}")