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
        dh_fraction="data/heating/DH_percent2020.h5"
    threads: 1
    resources: mem_mb=1000
    script: "scripts/DH_percent2020.py"

