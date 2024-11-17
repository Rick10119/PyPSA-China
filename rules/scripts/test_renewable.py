# SPDX-FileCopyrightText: : 2024 The PyPSA-China Authors
#
# SPDX-License-Identifier: MIT

import logging
import progressbar as pgb
import functools
import atlite
import xarray as xr
import geopandas as gpd
from atlite.gis import ExclusionContainer
import numpy as np
import time
import pandas as pd
from functions import pro_names

logger = logging.getLogger(__name__)

# 设置日志
logging.basicConfig(level=logging.INFO)

# 模拟 snakemake 对象
class MockSnakemake:
    class input:
        Build_up_raster = "data/landuse_availability/Build_up.tif"
        Grass_raster = "data/landuse_availability/Grass.tif"
        Bare_raster = "data/landuse_availability/Bare.tif"
        Shrubland_raster = "data/landuse_availability/Shrubland.tif"
        natura1 = 'data/landuse_availability/WDPA_WDOECM_Nov2024_Public_CHN_shp/WDPA_WDOECM_Nov2024_Public_CHN_shp_0/WDPA_WDOECM_Nov2024_Public_CHN_shp-polygons.shp'
        natura2 = 'data/landuse_availability/WDPA_WDOECM_Nov2024_Public_CHN_shp/WDPA_WDOECM_Nov2024_Public_CHN_shp_1/WDPA_WDOECM_Nov2024_Public_CHN_shp-polygons.shp'
        natura3 = 'data/landuse_availability/WDPA_WDOECM_Nov2024_Public_CHN_shp/WDPA_WDOECM_Nov2024_Public_CHN_shp_2/WDPA_WDOECM_Nov2024_Public_CHN_shp-polygons.shp'
        gebco = "data/landuse_availability/GEBCO_2024/gebco_2024.tif"
        provinces_shp = "data/province_shapes/CHN_adm1.shp"
        offshore_province_shapes = "data/resources/regions_offshore_province.geojson"
        offshore_shapes = "data/resources/regions_offshore.geojson"
        cutout = "cutouts/China-2020.nc"

    class output:
        solar_profile = "resources/profile_solar.nc"
        onwind_profile = "resources/profile_onwind.nc"
        offwind_profile = "resources/profile_offwind.nc"

    threads = 4
    config = {
        'atlite': {'show_progress': True},
        'Technique': {
            'solar': True,
            'onwind': True,
            'offwind': True
        },
        'renewable': {
            'solar': {
                'resource': {'method': 'pv'},
                'correction_factor': 1,
                'capacity_per_sqkm': 1.5
            },
            'onwind': {
                'resource': {'method': 'wind'},
                'correction_factor': 1,
                'capacity_per_sqkm': 5
            },
            'offwind': {
                'resource': {'method': 'wind'},
                'correction_factor': 1,
                'capacity_per_sqkm': 8,
                'max_depth': 50,
                'natura': True
            }
        }
    }

snakemake = MockSnakemake()

try:
    # 原始脚本的主要逻辑
    nprocesses = int(snakemake.threads)
    noprogress = not snakemake.config['atlite'].get('show_progress', True)

    cutout = atlite.Cutout(snakemake.input.cutout)
    cutout.prepare()

    # 读取省份数据
    provinces_shp = gpd.read_file(snakemake.input.provinces_shp)[['NAME_1', 'geometry']]

    # 替换名称
    provinces_shp.replace(to_replace={'Nei Mongol': 'InnerMongolia',
                                  'Xinjiang Uygur': 'Xinjiang',
                                  'Ningxia Hui': 'Ningxia',
                                  'Xizang': 'Tibet'}, inplace=True)

    # 设置索引
    provinces_shp.set_index('NAME_1', inplace=True)

    # 检查重复项
    duplicates = provinces_shp.index.duplicated()
    if duplicates.any():
        print("重复的省份名称：", provinces_shp.index[duplicates])
        provinces_shp = provinces_shp[~duplicates]

    # 确保 pro_names 唯一
    pro_names = list(set(pro_names))

    # 重新索引
    provinces_shp = provinces_shp.reindex(pro_names).rename_axis('bus')

    # ... 后续代码与原始脚本相同 ...

except Exception as e:
    logger.error(f"发生错误: {str(e)}", exc_info=True)
    raise
