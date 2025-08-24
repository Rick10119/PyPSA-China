#!/bin/bash
# No aluminum场景的模拟 (Flexibility: high, Demand: mid, Market: high, Year: 2030)

echo "开始运行 HMH_2030_no_aluminum 的模拟..."
echo "配置文件: configs/config_HMH_2030_no_aluminum.yaml"
echo "配置类型: no_aluminum配置"
echo "Scenario: Flexibility=high, Demand=mid, Market=high, Year=2030 (HMH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_HMH_2030_no_aluminum.yaml --cores 40

echo "完成 HMH_2030_no_aluminum 的模拟"
