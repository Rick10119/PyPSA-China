#!/bin/bash
# 80%容量比例的模拟 (Flexibility: non_constrained, Demand: mid, Market: mid, Year: 2050)

echo "开始运行 NMM_2050_80 的模拟..."
echo "配置文件: configs/config_NMM_2050_80p.yaml"
echo "配置类型: 80p配置"
echo "Scenario: Flexibility=non_constrained, Demand=mid, Market=mid, Year=2050 (NMM)"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_NMM_2050_80p.yaml --cores 40

echo "完成 NMM_2050_80 的模拟"
