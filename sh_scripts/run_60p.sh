#!/bin/bash
# 运行 60p 过剩产能保留比例的模拟

echo "开始运行 60p 过剩产能保留比例的模拟..."
echo "配置文件: configs/config_60p.yaml"
echo "版本号: 0814.4H.2-MMM-60p"
echo "Scenario: mid-mid-mid (MMM)"
echo "Cap: 60%, Actual: 70.4%"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_60p.yaml --cores 40

echo "完成 60p 过剩产能保留比例的模拟"
