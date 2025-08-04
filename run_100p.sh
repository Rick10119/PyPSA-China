#!/bin/bash
# 运行 100% 容量比例的模拟

echo "开始运行 100% 容量比例的模拟..."
echo "配置文件: config_100p.yaml"
echo "版本号: 0723.8H.4-100p"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile config_100p.yaml --cores 40

echo "完成 100% 容量比例的模拟"
