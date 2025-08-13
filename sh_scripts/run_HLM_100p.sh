#!/bin/bash
# 100%容量比例的模拟 (Flexibility: high, Demand: low, Market: mid)

echo "开始运行 HLM_100p 的模拟..."
echo "配置文件: configs/config_HLM_100p.yaml"
echo "配置类型: 100p配置"
echo "Scenario: Flexibility=high, Demand=low, Market=mid (HLM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_HLM_100p.yaml --cores 40

echo "完成 HLM_100p 的模拟"
