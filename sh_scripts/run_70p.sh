#!/bin/bash
# 运行 70p 过剩产能保留比例的模拟

echo "开始运行 70p 过剩产能保留比例的模拟..."
echo "配置文件: configs/config_70p.yaml"
echo "版本号: 0814.4H.2-MMM-70p"
echo "Scenario: mid-mid-mid (MMM)"
echo "Cap: 70%, Actual: 77.8%"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_70p.yaml --cores 40

echo "完成 70p 过剩产能保留比例的模拟"
