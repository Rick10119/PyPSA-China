#!/bin/bash
# 70p过剩产能保留比例的模拟 (Flexibility: high, Demand: mid, Market: mid, Year: 2050)

echo "开始运行 HMM_2050_70p 的模拟..."
echo "配置文件: configs/config_HMM_2050_70p.yaml"
echo "配置类型: 70p配置"
echo "Scenario: Flexibility=high, Demand=mid, Market=mid, Year=2050 (HMM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_HMM_2050_70p.yaml --cores 40

echo "完成 HMM_2050_70p 的模拟"
