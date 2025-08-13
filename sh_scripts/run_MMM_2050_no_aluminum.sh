#!/bin/bash
# No aluminum场景的模拟 (Flexibility: mid, Demand: mid, Market: mid, Year: 2050)

echo "开始运行 MMM_2050_no_aluminum 的模拟..."
echo "配置文件: configs/config_MMM_2050_no_aluminum.yaml"
echo "配置类型: no_aluminum配置"
echo "Scenario: Flexibility=mid, Demand=mid, Market=mid, Year=2050 (MMM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_MMM_2050_no_aluminum.yaml --cores 40

echo "完成 MMM_2050_no_aluminum 的模拟"
