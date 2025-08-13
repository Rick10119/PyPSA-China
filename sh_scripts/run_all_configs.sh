#!/bin/bash
# 批量运行所有配置的模拟

echo "开始批量运行所有配置的模拟..."
echo "总共72个配置: 4种flexibility × 3种demand × 3种market × 2种容量设置"
echo

config_count=0
total_configs=72


# 运行 LLL 场景
echo "=== 运行 LLL 场景 ==="
echo "Flexibility: low, Demand: low, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LLL_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LLL_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LLM 场景
echo "=== 运行 LLM 场景 ==="
echo "Flexibility: low, Demand: low, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LLM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LLM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LLH 场景
echo "=== 运行 LLH 场景 ==="
echo "Flexibility: low, Demand: low, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LLH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LLH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LML 场景
echo "=== 运行 LML 场景 ==="
echo "Flexibility: low, Demand: mid, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LML_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LML_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LMM 场景
echo "=== 运行 LMM 场景 ==="
echo "Flexibility: low, Demand: mid, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LMM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LMM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LMH 场景
echo "=== 运行 LMH 场景 ==="
echo "Flexibility: low, Demand: mid, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LMH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LMH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LHL 场景
echo "=== 运行 LHL 场景 ==="
echo "Flexibility: low, Demand: high, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LHL_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LHL_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LHM 场景
echo "=== 运行 LHM 场景 ==="
echo "Flexibility: low, Demand: high, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LHM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LHM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LHH 场景
echo "=== 运行 LHH 场景 ==="
echo "Flexibility: low, Demand: high, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LHH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LHH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MLL 场景
echo "=== 运行 MLL 场景 ==="
echo "Flexibility: mid, Demand: low, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MLL_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MLL_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MLM 场景
echo "=== 运行 MLM 场景 ==="
echo "Flexibility: mid, Demand: low, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MLM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MLM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MLH 场景
echo "=== 运行 MLH 场景 ==="
echo "Flexibility: mid, Demand: low, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MLH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MLH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MML 场景
echo "=== 运行 MML 场景 ==="
echo "Flexibility: mid, Demand: mid, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MML_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MML_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MMM 场景
echo "=== 运行 MMM 场景 ==="
echo "Flexibility: mid, Demand: mid, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MMM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MMM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MMH 场景
echo "=== 运行 MMH 场景 ==="
echo "Flexibility: mid, Demand: mid, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MMH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MMH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MHL 场景
echo "=== 运行 MHL 场景 ==="
echo "Flexibility: mid, Demand: high, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MHL_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MHL_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MHM 场景
echo "=== 运行 MHM 场景 ==="
echo "Flexibility: mid, Demand: high, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MHM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MHM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MHH 场景
echo "=== 运行 MHH 场景 ==="
echo "Flexibility: mid, Demand: high, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MHH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MHH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HLL 场景
echo "=== 运行 HLL 场景 ==="
echo "Flexibility: high, Demand: low, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HLL_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HLL_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HLM 场景
echo "=== 运行 HLM 场景 ==="
echo "Flexibility: high, Demand: low, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HLM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HLM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HLH 场景
echo "=== 运行 HLH 场景 ==="
echo "Flexibility: high, Demand: low, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HLH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HLH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HML 场景
echo "=== 运行 HML 场景 ==="
echo "Flexibility: high, Demand: mid, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HML_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HML_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HMM 场景
echo "=== 运行 HMM 场景 ==="
echo "Flexibility: high, Demand: mid, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HMM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HMM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HMH 场景
echo "=== 运行 HMH 场景 ==="
echo "Flexibility: high, Demand: mid, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HMH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HMH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HHL 场景
echo "=== 运行 HHL 场景 ==="
echo "Flexibility: high, Demand: high, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HHL_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HHL_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HHM 场景
echo "=== 运行 HHM 场景 ==="
echo "Flexibility: high, Demand: high, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HHM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HHM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HHH 场景
echo "=== 运行 HHH 场景 ==="
echo "Flexibility: high, Demand: high, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HHH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HHH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NLL 场景
echo "=== 运行 NLL 场景 ==="
echo "Flexibility: non_constrained, Demand: low, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NLL_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NLL_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NLM 场景
echo "=== 运行 NLM 场景 ==="
echo "Flexibility: non_constrained, Demand: low, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NLM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NLM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NLH 场景
echo "=== 运行 NLH 场景 ==="
echo "Flexibility: non_constrained, Demand: low, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NLH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NLH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NML 场景
echo "=== 运行 NML 场景 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NML_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NML_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NMM 场景
echo "=== 运行 NMM 场景 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NMM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NMM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NMH 场景
echo "=== 运行 NMH 场景 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NMH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NMH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NHL 场景
echo "=== 运行 NHL 场景 ==="
echo "Flexibility: non_constrained, Demand: high, Market: low, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NHL_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NHL_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NHM 场景
echo "=== 运行 NHM 场景 ==="
echo "Flexibility: non_constrained, Demand: high, Market: mid, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NHM_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NHM_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NHH 场景
echo "=== 运行 NHH 场景 ==="
echo "Flexibility: non_constrained, Demand: high, Market: high, Year: 2050"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NHH_2050_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NHH_2050_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


echo "所有配置的模拟已完成！"
echo "总共运行了 $total_configs 个配置"
echo "结果文件位于: results/version-*/"
