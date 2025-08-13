#!/bin/bash
# 90%容量比例的模拟 (Flexibility: mid, Demand: mid, Market: high, Year: 2040)

echo "开始运行 MMH_2040_90 的模拟..."
echo "配置文件: configs/config_MMH_2040_90p.yaml"
echo "配置类型: 90p配置"
echo "Scenario: Flexibility=mid, Demand=mid, Market=high, Year=2040 (MMH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_MMH_2040_90p.yaml --cores 40

echo "完成 MMH_2040_90 的模拟"
