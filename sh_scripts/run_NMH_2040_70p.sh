#!/bin/bash
# 70p过剩产能保留比例的模拟 (Flexibility: non_constrained, Demand: mid, Market: high, Year: 2040)

echo "开始运行 NMH_2040_70p 的模拟..."
echo "配置文件: configs/config_NMH_2040_70p.yaml"
echo "配置类型: 70p配置"
echo "Scenario: Flexibility=non_constrained, Demand=mid, Market=high, Year=2040 (NMH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_NMH_2040_70p.yaml --cores 40

echo "完成 NMH_2040_70p 的模拟"
