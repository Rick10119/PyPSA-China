#!/bin/bash
# 60p过剩产能保留比例的模拟 (Flexibility: non_constrained, Demand: mid, Market: mid, Year: 2040)

echo "开始运行 NMM_2040_60p 的模拟..."
echo "配置文件: configs/config_NMM_2040_60p.yaml"
echo "配置类型: 60p配置"
echo "Scenario: Flexibility=non_constrained, Demand=mid, Market=mid, Year=2040 (NMM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_NMM_2040_60p.yaml --cores 40

echo "完成 NMM_2040_60p 的模拟"
