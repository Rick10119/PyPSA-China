#!/bin/bash
# 运行 60% 容量比例的模拟

echo "开始运行 60% 容量比例的模拟..."
echo "配置文件: config_60p.yaml"
echo "版本号: 0723.8H.5-60p"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile config_60p.yaml --cores 40

echo "完成 60% 容量比例的模拟"
