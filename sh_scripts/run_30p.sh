#!/bin/bash
# 运行 30p 过剩产能保留比例的模拟

echo "开始运行 30p 过剩产能保留比例的模拟..."
echo "配置文件: configs/config_30p.yaml"
echo "版本号: 0814.4H.2-MMM-30p"
echo "Scenario: mid-mid-mid (MMM)"
echo "Cap: 30%, Actual: 48.1%"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_30p.yaml --cores 40

echo "完成 30p 过剩产能保留比例的模拟"
