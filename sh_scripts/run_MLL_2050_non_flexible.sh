#!/bin/bash
# Non-flexible场景的模拟 (Flexibility: mid, Demand: low, Market: low, Year: 2050)

echo "开始运行 MLL_2050_non_flexible 的模拟..."
echo "配置文件: configs/config_MLL_2050_non_flexible.yaml"
echo "配置类型: non-flexible配置"
echo "Scenario: Flexibility=mid, Demand=low, Market=low, Year=2050 (MLL)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_MLL_2050_non_flexible.yaml --cores 40

echo "完成 MLL_2050_non_flexible 的模拟"
