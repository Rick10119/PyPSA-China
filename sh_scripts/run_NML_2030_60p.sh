#!/bin/bash
# 60p过剩产能保留比例的模拟 (Flexibility: non_constrained, Demand: mid, Market: low, Year: 2030)

echo "开始运行 NML_2030_60p 的模拟..."
echo "配置文件: configs/config_NML_2030_60p.yaml"
echo "配置类型: 60p配置"
echo "Scenario: Flexibility=non_constrained, Demand=mid, Market=low, Year=2030 (NML)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_NML_2030_60p.yaml --cores 40

echo "完成 NML_2030_60p 的模拟"
