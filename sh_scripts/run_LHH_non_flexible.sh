#!/bin/bash
# Non-flexible场景的模拟 (Flexibility: low, Demand: high, Market: high)

echo "开始运行 LHH_non_flexible 的模拟..."
echo "配置文件: configs/config_LHH_non_flexible.yaml"
echo "配置类型: non-flexible配置"
echo "Scenario: Flexibility=low, Demand=high, Market=high (LHH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_LHH_non_flexible.yaml --cores 40

echo "完成 LHH_non_flexible 的模拟"
