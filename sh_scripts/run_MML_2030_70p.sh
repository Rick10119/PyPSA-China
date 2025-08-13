#!/bin/bash
# 70p过剩产能保留比例的模拟 (Flexibility: mid, Demand: mid, Market: low, Year: 2030)

echo "开始运行 MML_2030_70p 的模拟..."
echo "配置文件: configs/config_MML_2030_70p.yaml"
echo "配置类型: 70p配置"
echo "Scenario: Flexibility=mid, Demand=mid, Market=low, Year=2030 (MML)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_MML_2030_70p.yaml --cores 40

echo "完成 MML_2030_70p 的模拟"
