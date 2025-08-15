#!/bin/bash
# 运行 100p 过剩产能保留比例的模拟

echo "开始运行 100p 过剩产能保留比例的模拟..."
echo "配置文件: configs/config_100p.yaml"
echo "版本号: 0814.4H.2-MMM-100p"
echo "Scenario: mid-mid-mid (MMM)"
echo "Cap: 100%, Actual: 100.0%"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_100p.yaml --cores 40

echo "完成 100p 过剩产能保留比例的模拟"
