#!/bin/bash
# Non-flexible场景的模拟 (Flexibility: non_constrained, Demand: low, Market: mid)

echo "开始运行 NLM_non_flexible 的模拟..."
echo "配置文件: configs/config_NLM_non_flexible.yaml"
echo "配置类型: non-flexible配置"
echo "Scenario: Flexibility=non_constrained, Demand=low, Market=mid (NLM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_NLM_non_flexible.yaml --cores 40

echo "完成 NLM_non_flexible 的模拟"
