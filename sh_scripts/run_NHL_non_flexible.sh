#!/bin/bash
# Non-flexible场景的模拟 (Flexibility: non_constrained, Demand: high, Market: low)

echo "开始运行 NHL_non_flexible 的模拟..."
echo "配置文件: configs/config_NHL_non_flexible.yaml"
echo "配置类型: non-flexible配置"
echo "Scenario: Flexibility=non_constrained, Demand=high, Market=low (NHL)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_NHL_non_flexible.yaml --cores 40

echo "完成 NHL_non_flexible 的模拟"
