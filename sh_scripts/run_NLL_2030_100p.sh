#!/bin/bash
# 100%容量比例的模拟 (Flexibility: non_constrained, Demand: low, Market: low, Year: 2030)

echo "开始运行 NLL_2030_100p 的模拟..."
echo "配置文件: configs/config_NLL_2030_100p.yaml"
echo "配置类型: 100p配置"
echo "Scenario: Flexibility=non_constrained, Demand=low, Market=low, Year=2030 (NLL)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_NLL_2030_100p.yaml --cores 40

echo "完成 NLL_2030_100p 的模拟"
