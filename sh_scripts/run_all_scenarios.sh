#!/bin/bash
# 批量运行所有场景的模拟

echo "开始批量运行所有场景的模拟..."
echo "Scenario: mid-mid-mid (MMM)"
echo

# 运行不包含电解铝厂的场景
echo "=== 运行不包含电解铝厂的场景 ==="
./sh_scripts/run_no_aluminum.sh
echo

# 运行non-flexible场景
echo "=== 运行non-flexible场景 ==="
./sh_scripts/run_non_flexible.sh
echo

# 运行10p过剩产能保留比例
echo "=== 运行 10p 过剩产能保留比例 ==="
./sh_scripts/run_10p.sh
echo

# 运行20p过剩产能保留比例
echo "=== 运行 20p 过剩产能保留比例 ==="
./sh_scripts/run_20p.sh
echo

# 运行30p过剩产能保留比例
echo "=== 运行 30p 过剩产能保留比例 ==="
./sh_scripts/run_30p.sh
echo

# 运行40p过剩产能保留比例
echo "=== 运行 40p 过剩产能保留比例 ==="
./sh_scripts/run_40p.sh
echo

# 运行50p过剩产能保留比例
echo "=== 运行 50p 过剩产能保留比例 ==="
./sh_scripts/run_50p.sh
echo

# 运行60p过剩产能保留比例
echo "=== 运行 60p 过剩产能保留比例 ==="
./sh_scripts/run_60p.sh
echo

# 运行70p过剩产能保留比例
echo "=== 运行 70p 过剩产能保留比例 ==="
./sh_scripts/run_70p.sh
echo

# 运行80p过剩产能保留比例
echo "=== 运行 80p 过剩产能保留比例 ==="
./sh_scripts/run_80p.sh
echo

# 运行90p过剩产能保留比例
echo "=== 运行 90p 过剩产能保留比例 ==="
./sh_scripts/run_90p.sh
echo

# 运行100p过剩产能保留比例
echo "=== 运行 100p 过剩产能保留比例 ==="
./sh_scripts/run_100p.sh
echo

echo "所有场景的模拟已完成！"
echo "结果文件位于: results/version-MMM-*/"
