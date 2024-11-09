import xarray as xr
import os

# 定义文件路径
file1 = "../cutouts/63c61eaf523ec098207308f4544fdc5d.nc"  # 替换为实际文件名
file2 = "../cutouts/China-2020-1.nc"  # 替换为实际文件名
output_file = "../cutouts/China-2020.nc"

# 打开两个文件
ds1 = xr.open_dataset(file1)
ds2 = xr.open_dataset(file2)

# 合并数据集
combined_ds = xr.merge([ds1, ds2])

# 保存合并后的数据集
combined_ds.to_netcdf(output_file)

print(f"合并完成，输出文件为: {output_file}")