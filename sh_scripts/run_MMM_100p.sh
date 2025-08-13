#!/bin/bash
# 100%容量比例的模拟 (Flexibility: mid, Demand: mid, Market: mid)

echo "开始运行 MMM_100p 的模拟..."
echo "配置文件: configs/config_MMM_100p.yaml"
echo "配置类型: 100p配置"
echo "Scenario: Flexibility=mid, Demand=mid, Market=mid (MMM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_MMM_100p.yaml --cores 40

echo "完成 MMM_100p 的模拟"
