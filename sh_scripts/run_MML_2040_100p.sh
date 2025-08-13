#!/bin/bash
# 100%容量比例的模拟 (Flexibility: mid, Demand: mid, Market: low, Year: 2040)

echo "开始运行 MML_2040_100 的模拟..."
echo "配置文件: configs/config_MML_2040_100p.yaml"
echo "配置类型: 100p配置"
echo "Scenario: Flexibility=mid, Demand=mid, Market=low, Year=2040 (MML)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_MML_2040_100p.yaml --cores 40

echo "完成 MML_2040_100 的模拟"
