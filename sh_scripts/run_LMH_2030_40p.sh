#!/bin/bash
# 40%容量比例的模拟 (Flexibility: low, Demand: mid, Market: high, Year: 2030)

echo "开始运行 LMH_2030_40 的模拟..."
echo "配置文件: configs/config_LMH_2030_40p.yaml"
echo "配置类型: 40p配置"
echo "Scenario: Flexibility=low, Demand=mid, Market=high, Year=2030 (LMH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_LMH_2030_40p.yaml --cores 40

echo "完成 LMH_2030_40 的模拟"
