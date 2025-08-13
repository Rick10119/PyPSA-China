#!/bin/bash
# 100%容量比例的模拟 (Flexibility: low, Demand: mid, Market: high)

echo "开始运行 LMH_100p 的模拟..."
echo "配置文件: configs/config_LMH_100p.yaml"
echo "配置类型: 100p配置"
echo "Scenario: Flexibility=low, Demand=mid, Market=high (LMH)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_LMH_100p.yaml --cores 40

echo "完成 LMH_100p 的模拟"
