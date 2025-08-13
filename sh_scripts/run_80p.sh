#!/bin/bash
# 运行 80% 容量比例的模拟

echo "开始运行 80% 容量比例的模拟..."
echo "配置文件: configs/config_80p.yaml"
echo "版本号: 0811.1H.1-MHM-80p"
echo "Scenario: mid-high-mid (MHM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_80p.yaml --cores 40

echo "完成 80% 容量比例的模拟"
