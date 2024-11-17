# added rules

rule build_offshore_province:
    input:
        offshore_shapes="data/resources/regions_offshore.geojson"
    output:
        offshore_province="data/resources/regions_offshore_province.geojson"
    threads: 1
    resources: mem_mb=1000
    script: "scripts/build_offshore_province.py"

