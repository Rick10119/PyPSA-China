#!/bin/bash
# Non-flexible场景的模拟 (Flexibility: high, Demand: high, Market: low)

echo "开始运行 HHL_non_flexible 的模拟..."
echo "配置文件: configs/config_HHL_non_flexible.yaml"
echo "配置类型: non-flexible配置"
echo "Scenario: Flexibility=high, Demand=high, Market=low (HHL)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_HHL_non_flexible.yaml --cores 40

echo "完成 HHL_non_flexible 的模拟"
