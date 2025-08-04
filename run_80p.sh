#!/bin/bash
# 运行 80% 容量比例的模拟

echo "开始运行 80% 容量比例的模拟..."
echo "配置文件: config_80p.yaml"
echo "版本号: 0723.8H.4-80p"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile config_80p.yaml --cores 40

echo "完成 80% 容量比例的模拟"
