#!/bin/bash
# 100%容量比例的模拟 (Flexibility: low, Demand: low, Market: low, Year: 2040)

echo "开始运行 LLL_2040_100p 的模拟..."
echo "配置文件: configs/config_LLL_2040_100p.yaml"
echo "配置类型: 100p配置"
echo "Scenario: Flexibility=low, Demand=low, Market=low, Year=2040 (LLL)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_LLL_2040_100p.yaml --cores 40

echo "完成 LLL_2040_100p 的模拟"
