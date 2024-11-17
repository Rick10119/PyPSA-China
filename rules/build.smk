# added rules
rule build_tif_with_nc:
    input:
        nc_file="{data_dir}/{file}.nc"
    output:
        tif_file="{data_dir}/{file}.tif"
    threads: 2
    resources: mem_mb=50000
    script: "scripts/build_tif_with_nc.py"

rule build_offshore_province:
    input:
        offshore_shapes="data/resources/regions_offshore.geojson"
    output:
        offshore_province="data/resources/regions_offshore_province.geojson"
    threads: 1
    resources: mem_mb=1000
    script: "scripts/build_offshore_province.py"