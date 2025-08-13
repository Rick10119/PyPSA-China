#!/bin/bash
# 20%容量比例的模拟 (Flexibility: mid, Demand: mid, Market: mid, Year: 2040)

echo "开始运行 MMM_2040_20 的模拟..."
echo "配置文件: configs/config_MMM_2040_20p.yaml"
echo "配置类型: 20p配置"
echo "Scenario: Flexibility=mid, Demand=mid, Market=mid, Year=2040 (MMM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_MMM_2040_20p.yaml --cores 40

echo "完成 MMM_2040_20 的模拟"
