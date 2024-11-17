import psutil
import os
import xarray as xr
import rioxarray
import gc
from snakemake.shell import shell

def monitor_memory():
    process = psutil.Process()
    print(f"当前内存使用: {process.memory_info().rss / 1024 / 1024 / 1024:.2f} GB")

# 从 snakemake 对象获取输入输出文件路径
input_nc = snakemake.input.nc_file
output_tif = snakemake.output.tif_file

try:
    monitor_memory()
    print("开始读取NC文件...")
    
    # 使用chunks参数分块读取
    nc_data = xr.open_dataset(input_nc, chunks={'x': 1000, 'y': 1000})
    
    monitor_memory()
    print("NC文件读取完成，开始转换...")
    
    # 使用优化参数保存tiff
    nc_data.rio.to_raster(
        output_tif,
        tiled=True,
        blockxsize=256,
        blockysize=256,
        compress='LZW'
    )
    
    # 清理内存
    nc_data = None
    gc.collect()
    
    monitor_memory()
    print("转换完成")

except Exception as e:
    print(f"处理出错: {str(e)}")
    raise
finally:
    # 最后清理一次内存
    gc.collect()