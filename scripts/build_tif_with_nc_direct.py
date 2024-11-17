from osgeo import gdal
import psutil
import os
import gc

def limit_memory(max_gb=50):
    process = psutil.Process()
    # 设置软限制为50GB
    soft_limit = max_gb * 1024 * 1024 * 1024  # 转换为字节
    
    def check_memory():
        current_memory = process.memory_info().rss
        if current_memory > soft_limit:
            print(f"警告：内存使用超过{max_gb}GB，正在清理...")
            gc.collect()
            if process.memory_info().rss > soft_limit:
                raise MemoryError(f"内存使用超过{max_gb}GB限制")
    
    return check_memory

# 文件路径
input_nc = r'data/landuse_availability/GEBCO_2021/gebco_2021.nc'
output_tif = r'data/landuse_availability/GEBCO_2021/gebco_2021.tif'

try:
    print("开始转换...")
    check_memory = limit_memory(50)  # 设置50GB限制
    
    # 配置GDAL使用更少的内存
    gdal.SetCacheMax(1024 * 1024 * 1024)  # 设置GDAL缓存为1GB
    
    # 使用 GDAL 转换，设置分块处理
    gdal.Translate(
        output_tif,
        input_nc,
        format='GTiff',
        creationOptions=[
            'COMPRESS=LZW',
            'TILED=YES',
            'BLOCKXSIZE=256',  # 设置块大小
            'BLOCKYSIZE=256',
            'NUM_THREADS=ALL_CPUS'  # 使用所有CPU加速处理
        ]
    )
    
    # 检查内存使用
    check_memory()
    print("转换完成")

except Exception as e:
    print(f"处理出错: {str(e)}")
    raise
finally:
    # 清理内存
    gc.collect()