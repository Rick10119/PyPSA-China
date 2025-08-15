#!/bin/bash
# 运行 20p 过剩产能保留比例的模拟

echo "开始运行 20p 过剩产能保留比例的模拟..."
echo "配置文件: configs/config_20p.yaml"
echo "版本号: 0814.4H.2-MMM-20p"
echo "Scenario: mid-mid-mid (MMM)"
echo "Cap: 20%, Actual: 40.7%"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_20p.yaml --cores 40

echo "完成 20p 过剩产能保留比例的模拟"
