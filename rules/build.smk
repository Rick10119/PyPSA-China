# added rules

rule build_offshore_province:
    input:
        offshore_shapes="data/resources/regions_offshore.geojson"
    output:
        offshore_province="data/resources/regions_offshore_province.geojson"
    threads: 1
    resources: mem_mb=1000
    script: "scripts/build_offshore_province.py"

rule build_SPH_demand:
    input:
        population="data/population/population.h5"
    output:
        sph_demand="data/heating/SPH_2020.csv"
    threads: 1
    resources: mem_mb=1000
    script: "scripts/build_SPH_demand.py"

rule build_DH_fraction:
    input:
        population="data/population/population.h5",
        sph_demand="data/heating/SPH_2020.csv"
    output:
        dh_fraction="data/heating/DH_fraction_2020.h5"
    threads: 1
    resources: mem_mb=1000
    script: "scripts/DH_fraction_2020.py"

# 添加水电出力特性生成规则

rule build_hydro_profiles:
    input:
        config="config.yaml",  # 配置文件
        dams_data="data/hydro/dams_large.csv",  # 大坝数据
        inflow_data="data/hydro/daily_hydro_inflow_per_dam_1979_2016_m3.pickle"  # 入流数据
    output:
        hydro_p_max_pu="data/p_nom/hydro_p_max_pu.h5"  # 生成的水电出力特性文件
    log:
        "logs/build_hydro_profiles.log"  # 日志文件
    benchmark:
        "benchmarks/build_hydro_profiles"  # 性能基准
    threads: 1
    resources:
        mem_mb=2000  # 资源限制
    script:
        "rules/scripts/build_hydro.py"  # 执行的脚本

# 添加省份形状生成规则

rule build_province_shape:
    input:
        gadm_data="https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_CHN.gpkg"  # GADM 数据 URL
    output:
        province_shapes="data/province_shapes/CHN_adm1.shp"  # 生成的省份形状文件
    log:
        "logs/build_province_shape.log"  # 日志文件
    benchmark:
        "benchmarks/build_province_shape"  # 性能基准
    threads: 1
    resources:
        mem_mb=2000  # 资源限制
    script:
        "rules/scripts/build_provence_shape.py"  # 执行的脚本

# 添加负荷天气年生成规则

rule build_load_weatheryears:
    input:
        load_before="data/load/load_2040_weatheryears_2020_2060_TWh.h5",  # 参考文件（2040年）
        load_after="data/load/load_2055_weatheryears_2020_2060_TWh.h5"  # 参考文件（2055年）
    output:
        load_target="data/load/load_{target_year}_weatheryears_2020_2060_TWh.h5"  # 生成的目标年份负荷文件
    params:
        target_years=[2045, 2050]  # 目标年份列表
    log:
        "logs/build_load_weatheryears.log"  # 日志文件
    benchmark:
        "benchmarks/build_load_weatheryears"  # 性能基准
    threads: 1
    resources:
        mem_mb=2000  # 资源限制
    script:
        "rules/scripts/build_load_weatheryaers.py"  # 执行的脚本



