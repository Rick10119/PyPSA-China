from osgeo import gdal
import os
import glob

def merge_tifs(input_dir, output_tif):
    # 获取目录下所有的tif文件
    tif_files = glob.glob(os.path.join(input_dir, "*.tif"))
    print(f"找到以下tif文件：")
    for f in tif_files:
        print(f"- {f}")
    
    if not tif_files:
        raise Exception(f"在目录 {input_dir} 中没有找到tif文件")
    
    # 构建vrt文件（虚拟数据集）
    vrt = gdal.BuildVRT("merged.vrt", tif_files)
    
    # 转换为tif，添加 BIGTIFF=YES
    gdal.Translate(
        output_tif,
        vrt,
        format='GTiff',
        creationOptions=[
            'BIGTIFF=YES',
            'COMPRESS=LZW',
            'TILED=YES',
            'BLOCKXSIZE=256',
            'BLOCKYSIZE=256',
            'NUM_THREADS=ALL_CPUS'
        ]
    )
    
    # 清理临时文件
    vrt = None
    if os.path.exists("merged.vrt"):
        os.remove("merged.vrt")

try:
    input_dir = "data/landuse_availability/GEBCO_2024"
    output_tif = "data/landuse_availability/GEBCO_2024/gebco_2024.tif"
    
    print("开始合并tif文件...")
    merge_tifs(input_dir, output_tif)
    print(f"合并完成，输出文件：{output_tif}")

except Exception as e:
    print(f"处理出错: {str(e)}")
    raise