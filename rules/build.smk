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

rule build_co2_totals:
    output:
        "data/co2_totals.h5"
    script:
        "scripts/build_co2_totals.py"

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



