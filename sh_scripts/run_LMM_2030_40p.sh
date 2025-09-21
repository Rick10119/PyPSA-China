#!/bin/bash
# 40p过剩产能保留比例的模拟 (Flexibility: low, Demand: mid, Market: mid, Year: 2030)

echo "开始运行 LMM_2030_40p 的模拟..."
echo "配置文件: configs/config_LMM_2030_40p.yaml"
echo "配置类型: 40p配置"
echo "Scenario: Flexibility=low, Demand=mid, Market=mid, Year=2030 (LMM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_LMM_2030_40p.yaml --cores 40

echo "完成 LMM_2030_40p 的模拟"
