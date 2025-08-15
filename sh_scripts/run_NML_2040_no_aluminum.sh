#!/bin/bash
# No aluminum场景的模拟 (Flexibility: non_constrained, Demand: mid, Market: low, Year: 2040)

echo "开始运行 NML_2040_no_aluminum 的模拟..."
echo "配置文件: configs/config_NML_2040_no_aluminum.yaml"
echo "配置类型: no_aluminum配置"
echo "Scenario: Flexibility=non_constrained, Demand=mid, Market=low, Year=2040 (NML)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_NML_2040_no_aluminum.yaml --cores 40

echo "完成 NML_2040_no_aluminum 的模拟"
