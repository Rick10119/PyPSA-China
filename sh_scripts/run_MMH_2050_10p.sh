#!/bin/bash
# 10%容量比例的模拟 (Flexibility: mid, Demand: mid, Market: high, Year: 2050)

echo "开始运行 MMH_2050_10 的模拟..."
echo "配置文件: configs/config_MMH_2050_10p.yaml"
echo "配置类型: 10p配置"
echo "Scenario: Flexibility=mid, Demand=mid, Market=high, Year=2050 (MMH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_MMH_2050_10p.yaml --cores 40

echo "完成 MMH_2050_10 的模拟"
