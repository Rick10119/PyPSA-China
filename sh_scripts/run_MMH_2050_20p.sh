#!/bin/bash
# 20p过剩产能保留比例的模拟 (Flexibility: mid, Demand: mid, Market: high, Year: 2050)

echo "开始运行 MMH_2050_20p 的模拟..."
echo "配置文件: configs/config_MMH_2050_20p.yaml"
echo "配置类型: 20p配置"
echo "Scenario: Flexibility=mid, Demand=mid, Market=high, Year=2050 (MMH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_MMH_2050_20p.yaml --cores 40

echo "完成 MMH_2050_20p 的模拟"
