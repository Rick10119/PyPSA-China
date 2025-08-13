#!/bin/bash
# 40%容量比例的模拟 (Flexibility: high, Demand: mid, Market: high, Year: 2040)

echo "开始运行 HMH_2040_40 的模拟..."
echo "配置文件: configs/config_HMH_2040_40p.yaml"
echo "配置类型: 40p配置"
echo "Scenario: Flexibility=high, Demand=mid, Market=high, Year=2040 (HMH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_HMH_2040_40p.yaml --cores 40

echo "完成 HMH_2040_40 的模拟"
