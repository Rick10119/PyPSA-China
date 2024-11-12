import xarray as xr
import pandas as pd

# 读取 NetCDF 文件
file_path = r'C:\Users\dell\Documents\GitHub\PyPSA-China\resources\profile_onwind.nc'
data = xr.open_dataset(file_path)

# 打印数据集的变量
print(data.variables)