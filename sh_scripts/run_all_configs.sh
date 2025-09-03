#!/bin/bash
# 批量运行所有配置的模拟

echo "开始批量运行所有配置的模拟..."
echo "总共468个配置: 4种flexibility × 1种demand × 3种market × 3种年份 × 13种容量设置"
echo

config_count=0
total_configs=468


# 运行 LML_2050 场景
echo "=== 运行 LML_2050 场景 ==="
echo "Flexibility: low, Demand: mid, Market: low, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LML_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_LML_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_LML_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_LML_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_LML_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_LML_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_LML_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_LML_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_LML_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_LML_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_LML_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_LML_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LML_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LMM_2050 场景
echo "=== 运行 LMM_2050 场景 ==="
echo "Flexibility: low, Demand: mid, Market: mid, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LMM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_LMM_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_LMM_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_LMM_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_LMM_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_LMM_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_LMM_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_LMM_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_LMM_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_LMM_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_LMM_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_LMM_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LMM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LMH_2050 场景
echo "=== 运行 LMH_2050 场景 ==="
echo "Flexibility: low, Demand: mid, Market: high, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LMH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_LMH_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_LMH_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_LMH_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_LMH_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_LMH_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_LMH_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_LMH_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_LMH_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_LMH_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_LMH_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_LMH_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LMH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MML_2050 场景
echo "=== 运行 MML_2050 场景 ==="
echo "Flexibility: mid, Demand: mid, Market: low, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MML_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_MML_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_MML_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_MML_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_MML_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_MML_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_MML_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_MML_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_MML_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_MML_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_MML_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_MML_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MML_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MMM_2050 场景
echo "=== 运行 MMM_2050 场景 ==="
echo "Flexibility: mid, Demand: mid, Market: mid, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MMM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_MMM_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_MMM_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_MMM_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_MMM_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_MMM_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_MMM_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_MMM_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_MMM_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_MMM_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_MMM_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_MMM_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MMM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MMH_2050 场景
echo "=== 运行 MMH_2050 场景 ==="
echo "Flexibility: mid, Demand: mid, Market: high, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MMH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_MMH_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_MMH_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_MMH_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_MMH_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_MMH_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_MMH_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_MMH_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_MMH_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_MMH_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_MMH_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_MMH_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MMH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HML_2050 场景
echo "=== 运行 HML_2050 场景 ==="
echo "Flexibility: high, Demand: mid, Market: low, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HML_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_HML_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_HML_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_HML_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_HML_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_HML_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_HML_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_HML_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_HML_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_HML_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_HML_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_HML_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HML_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HMM_2050 场景
echo "=== 运行 HMM_2050 场景 ==="
echo "Flexibility: high, Demand: mid, Market: mid, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HMM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_HMM_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_HMM_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_HMM_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_HMM_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_HMM_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_HMM_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_HMM_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_HMM_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_HMM_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_HMM_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_HMM_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HMM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HMH_2050 场景
echo "=== 运行 HMH_2050 场景 ==="
echo "Flexibility: high, Demand: mid, Market: high, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HMH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_HMH_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_HMH_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_HMH_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_HMH_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_HMH_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_HMH_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_HMH_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_HMH_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_HMH_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_HMH_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_HMH_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HMH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NML_2050 场景
echo "=== 运行 NML_2050 场景 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: low, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NML_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_NML_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_NML_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_NML_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_NML_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_NML_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_NML_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_NML_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_NML_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_NML_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_NML_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_NML_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NML_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NMM_2050 场景
echo "=== 运行 NMM_2050 场景 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: mid, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NMM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_NMM_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_NMM_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_NMM_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_NMM_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_NMM_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_NMM_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_NMM_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_NMM_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_NMM_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_NMM_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_NMM_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NMM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NMH_2050 场景
echo "=== 运行 NMH_2050 场景 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: high, Year: 2050"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NMH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_NMH_2050_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 5p 配置
echo "--- 运行 5p 配置 ---"
./sh_scripts/run_NMH_2050_5p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 10p 配置
echo "--- 运行 10p 配置 ---"
./sh_scripts/run_NMH_2050_10p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 20p 配置
echo "--- 运行 20p 配置 ---"
./sh_scripts/run_NMH_2050_20p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 30p 配置
echo "--- 运行 30p 配置 ---"
./sh_scripts/run_NMH_2050_30p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 40p 配置
echo "--- 运行 40p 配置 ---"
./sh_scripts/run_NMH_2050_40p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 50p 配置
echo "--- 运行 50p 配置 ---"
./sh_scripts/run_NMH_2050_50p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 60p 配置
echo "--- 运行 60p 配置 ---"
./sh_scripts/run_NMH_2050_60p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 70p 配置
echo "--- 运行 70p 配置 ---"
./sh_scripts/run_NMH_2050_70p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 80p 配置
echo "--- 运行 80p 配置 ---"
./sh_scripts/run_NMH_2050_80p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 90p 配置
echo "--- 运行 90p 配置 ---"
./sh_scripts/run_NMH_2050_90p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 100p 配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NMH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


echo "所有配置的模拟已完成！"
echo "总共运行了 $total_configs 个配置"
echo "结果文件位于: results/version-*/"
