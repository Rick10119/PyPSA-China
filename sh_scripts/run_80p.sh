#!/bin/bash
# 运行 80p 过剩产能保留比例的模拟

echo "开始运行 80p 过剩产能保留比例的模拟..."
echo "配置文件: configs/config_80p.yaml"
echo "版本号: 0814.4H.2-MMM-80p"
echo "Scenario: mid-mid-mid (MMM)"
echo "Cap: 80%, Actual: 85.2%"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_80p.yaml --cores 40

echo "完成 80p 过剩产能保留比例的模拟"
