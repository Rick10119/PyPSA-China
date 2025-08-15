#!/bin/bash
# No aluminum场景的模拟 (Flexibility: low, Demand: mid, Market: low, Year: 2030)

echo "开始运行 LML_2030_no_aluminum 的模拟..."
echo "配置文件: configs/config_LML_2030_no_aluminum.yaml"
echo "配置类型: no_aluminum配置"
echo "Scenario: Flexibility=low, Demand=mid, Market=low, Year=2030 (LML)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_LML_2030_no_aluminum.yaml --cores 40

echo "完成 LML_2030_no_aluminum 的模拟"
