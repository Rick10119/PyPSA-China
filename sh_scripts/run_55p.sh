#!/bin/bash
# 运行 55% 容量比例的模拟

echo "开始运行 55% 容量比例的模拟..."
echo "配置文件: configs/config_55p.yaml"
echo "版本号: 0809.1H.1-55p"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_55p.yaml --cores 40

echo "完成 55% 容量比例的模拟"
