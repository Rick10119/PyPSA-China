#!/bin/bash
# Non-flexible场景的模拟 (Flexibility: mid, Demand: mid, Market: low)

echo "开始运行 MML_non_flexible 的模拟..."
echo "配置文件: configs/config_MML_non_flexible.yaml"
echo "配置类型: non-flexible配置"
echo "Scenario: Flexibility=mid, Demand=mid, Market=low (MML)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_MML_non_flexible.yaml --cores 40

echo "完成 MML_non_flexible 的模拟"
