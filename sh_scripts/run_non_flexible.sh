#!/bin/bash
# 运行non-flexible场景的模拟

echo "开始运行non-flexible场景的模拟..."
echo "配置文件: configs/config_non_flexible.yaml"
echo "版本号: 0811.1H.1-non-flexible"
echo "设置: add_aluminum=False, only_other_load=False"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_non_flexible.yaml --cores 40

echo "完成non-flexible场景的模拟"
