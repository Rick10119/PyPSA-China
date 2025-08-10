#!/bin/bash
# 运行不包含电解铝厂的模拟

echo "开始运行不包含电解铝厂的模拟..."
echo "配置文件: configs/config_no_aluminum.yaml"
echo "版本号: 0809.2H.2-no-aluminum"
echo "设置: add_aluminum=False, only_other_load=True"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_no_aluminum.yaml --cores 40

echo "完成不包含电解铝厂的模拟"
