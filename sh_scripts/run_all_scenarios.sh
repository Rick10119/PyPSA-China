#!/bin/bash
# 批量运行所有场景的模拟

echo "开始批量运行所有场景的模拟..."
echo "Scenario: mid-high-mid (MHM)"
echo

# 运行不包含电解铝厂的场景
echo "=== 运行不包含电解铝厂的场景 ==="
./sh_scripts/run_no_aluminum.sh
echo

# 运行non-flexible场景
echo "=== 运行non-flexible场景 ==="
./sh_scripts/run_non_flexible.sh
echo

# 运行100%容量比例
echo "=== 运行 100% 容量比例 ==="
./sh_scripts/run_100p.sh
echo

# 运行90%容量比例
echo "=== 运行 90% 容量比例 ==="
./sh_scripts/run_90p.sh
echo

# 运行80%容量比例
echo "=== 运行 80% 容量比例 ==="
./sh_scripts/run_80p.sh
echo

# 运行70%容量比例
echo "=== 运行 70% 容量比例 ==="
./sh_scripts/run_70p.sh
echo

# 运行60%容量比例
echo "=== 运行 60% 容量比例 ==="
./sh_scripts/run_60p.sh
echo

# 运行55%容量比例
echo "=== 运行 55% 容量比例 ==="
./sh_scripts/run_55p.sh
echo

echo "所有场景的模拟已完成！"
echo "结果文件位于: results/version-MHM-*/"
