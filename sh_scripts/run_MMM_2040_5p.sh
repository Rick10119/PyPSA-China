#!/bin/bash
# 5p过剩产能保留比例的模拟 (Flexibility: mid, Demand: mid, Market: mid, Year: 2040)

echo "开始运行 MMM_2040_5p 的模拟..."
echo "配置文件: configs/config_MMM_2040_5p.yaml"
echo "配置类型: 5p配置"
echo "Scenario: Flexibility=mid, Demand=mid, Market=mid, Year=2040 (MMM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_MMM_2040_5p.yaml --cores 40

echo "完成 MMM_2040_5p 的模拟"
