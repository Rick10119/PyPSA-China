#!/bin/bash
# 100%容量比例的模拟 (Flexibility: low, Demand: low, Market: mid, Year: 2040)

echo "开始运行 LLM_2040_100p 的模拟..."
echo "配置文件: configs/config_LLM_2040_100p.yaml"
echo "配置类型: 100p配置"
echo "Scenario: Flexibility=low, Demand=low, Market=mid, Year=2040 (LLM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_LLM_2040_100p.yaml --cores 40

echo "完成 LLM_2040_100p 的模拟"
