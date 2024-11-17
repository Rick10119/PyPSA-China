import geopandas as gpd
import json
from shapely.geometry import Polygon, MultiPolygon, mapping
from shapely.ops import unary_union
import os

def simplify_geometry(geom, tolerance=0.01):
    """简化几何形状，减少点的数量"""
    return geom.simplify(tolerance, preserve_topology=True)

# 读取现有的 offshore regions
offshore = gpd.read_file(snakemake.input.offshore_shapes)

# 定义沿海省份及其大致范围
province_bounds = {
    "Liaoning": [120.0, 38.5, 124.5, 41.0],
    "Hebei": [118.0, 38.0, 120.0, 40.0],
    "Tianjin": [117.5, 38.5, 118.5, 39.5],
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
        # 获取该省份的海上区域
        province_offshore = offshore[mask].copy()
        # 合并重叠的几何形状
        merged_geom = unary_union(province_offshore.geometry)
        # 简化几何形状
        simplified_geom = simplify_geometry(merged_geom)
        
        # 创建特征
        feature = {
            "type": "Feature",
            "properties": {
                "bus": province,  # 将 "province" 改为 "bus"
                "name": f"{province}_Offshore"
            },
            "geometry": mapping(simplified_geom)
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
    json.dump(offshore_province, f)

# 验证文件大小
print(f"原始文件大小: {os.path.getsize(snakemake.input.offshore_shapes) / (1024*1024):.2f} MB")
print(f"生成文件大小: {os.path.getsize(snakemake.output.offshore_province) / (1024*1024):.2f} MB")