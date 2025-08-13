#!/bin/bash
# 50p过剩产能保留比例的模拟 (Flexibility: high, Demand: mid, Market: mid, Year: 2040)

echo "开始运行 HMM_2040_50p 的模拟..."
echo "配置文件: configs/config_HMM_2040_50p.yaml"
echo "配置类型: 50p配置"
echo "Scenario: Flexibility=high, Demand=mid, Market=mid, Year=2040 (HMM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_HMM_2040_50p.yaml --cores 40

echo "完成 HMM_2040_50p 的模拟"
