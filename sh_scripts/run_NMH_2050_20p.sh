#!/bin/bash
# 20%容量比例的模拟 (Flexibility: non_constrained, Demand: mid, Market: high, Year: 2050)

echo "开始运行 NMH_2050_20 的模拟..."
echo "配置文件: configs/config_NMH_2050_20p.yaml"
echo "配置类型: 20p配置"
echo "Scenario: Flexibility=non_constrained, Demand=mid, Market=high, Year=2050 (NMH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_NMH_2050_20p.yaml --cores 40

echo "完成 NMH_2050_20 的模拟"
