#!/bin/bash
# 批量运行所有容量比例的模拟

echo "开始批量运行所有容量比例的模拟..."
echo

# 运行100%容量比例
echo "=== 运行 100% 容量比例 ==="
./run_100p.sh
echo

# 运行90%容量比例
echo "=== 运行 90% 容量比例 ==="
./run_90p.sh
echo

# 运行80%容量比例
echo "=== 运行 80% 容量比例 ==="
./run_80p.sh
echo

# 运行70%容量比例
echo "=== 运行 70% 容量比例 ==="
./run_70p.sh
echo

# 运行60%容量比例
echo "=== 运行 60% 容量比例 ==="
./run_60p.sh
echo

# 运行55%容量比例
echo "=== 运行 55% 容量比例 ==="
./run_55p.sh
echo

echo "所有容量比例的模拟已完成！"
echo "结果文件位于: results/version-0723.8H.5-*p/"
