#!/bin/bash
# 运行 90% 容量比例的模拟

echo "开始运行 90% 容量比例的模拟..."
echo "配置文件: configs/config_90p.yaml"
echo "版本号: 0809.1H.1-90p"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_90p.yaml --cores 40

echo "完成 90% 容量比例的模拟"
