#!/bin/bash
# 80p过剩产能保留比例的模拟 (Flexibility: low, Demand: mid, Market: high, Year: 2030)

echo "开始运行 LMH_2030_80p 的模拟..."
echo "配置文件: configs/config_LMH_2030_80p.yaml"
echo "配置类型: 80p配置"
echo "Scenario: Flexibility=low, Demand=mid, Market=high, Year=2030 (LMH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_LMH_2030_80p.yaml --cores 40

echo "完成 LMH_2030_80p 的模拟"
