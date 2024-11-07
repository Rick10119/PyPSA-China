import geopandas as gpd
import os

# 创建目录
os.makedirs("../data/province_shapes", exist_ok=True)

# 使用新的 GADM 数据 URL
url = "https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_CHN.gpkg"
gdf = gpd.read_file(url, layer=1)  # layer=1 表示省级行政区划

# 打印列名，看看实际的结构
print("数据列名：", gdf.columns.tolist())

# 筛选并重命名省份
gdf = gdf[['NAME_1', 'geometry']]
gdf = gdf.replace({
    'Nei Mongol': 'InnerMongolia',
    'Xinjiang Uygur': 'Xinjiang',
    'Ningxia Hui': 'Ningxia',
    'Xizang': 'Tibet'
})

# 保存为 shapefile（这会自动生成所有必需的文件，包括 .shx）
gdf.to_file("../data/province_shapes/CHN_adm1.shp")