#!/bin/bash
# 60%容量比例的模拟 (Flexibility: low, Demand: mid, Market: low, Year: 2050)

echo "开始运行 LML_2050_60 的模拟..."
echo "配置文件: configs/config_LML_2050_60p.yaml"
echo "配置类型: 60p配置"
echo "Scenario: Flexibility=low, Demand=mid, Market=low, Year=2050 (LML)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_LML_2050_60p.yaml --cores 40

echo "完成 LML_2050_60 的模拟"
