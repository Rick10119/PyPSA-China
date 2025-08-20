#!/bin/bash
# No aluminum场景的模拟 (Flexibility: non_constrained, Demand: mid, Market: high, Year: 2050)

echo "开始运行 NMH_2050_no_aluminum 的模拟..."
echo "配置文件: configs/config_NMH_2050_no_aluminum.yaml"
echo "配置类型: no_aluminum配置"
echo "Scenario: Flexibility=non_constrained, Demand=mid, Market=high, Year=2050 (NMH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_NMH_2050_no_aluminum.yaml --cores 40

echo "完成 NMH_2050_no_aluminum 的模拟"
