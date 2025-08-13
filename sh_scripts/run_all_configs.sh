#!/bin/bash
# 批量运行所有配置的模拟

echo "开始批量运行所有配置的模拟..."
echo "总共72个配置: 4种flexibility × 3种demand × 3种market × 2种容量设置"
echo

config_count=0
total_configs=72


# 运行 LLL 场景
echo "=== 运行 LLL 场景 ==="
echo "Flexibility: low, Demand: low, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LLL_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LLL_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LLM 场景
echo "=== 运行 LLM 场景 ==="
echo "Flexibility: low, Demand: low, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LLM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LLM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LLH 场景
echo "=== 运行 LLH 场景 ==="
echo "Flexibility: low, Demand: low, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LLH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LLH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LML 场景
echo "=== 运行 LML 场景 ==="
echo "Flexibility: low, Demand: mid, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LML_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LML_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LMM 场景
echo "=== 运行 LMM 场景 ==="
echo "Flexibility: low, Demand: mid, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LMM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LMM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LMH 场景
echo "=== 运行 LMH 场景 ==="
echo "Flexibility: low, Demand: mid, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LMH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LMH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LHL 场景
echo "=== 运行 LHL 场景 ==="
echo "Flexibility: low, Demand: high, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LHL_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LHL_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LHM 场景
echo "=== 运行 LHM 场景 ==="
echo "Flexibility: low, Demand: high, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LHM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LHM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LHH 场景
echo "=== 运行 LHH 场景 ==="
echo "Flexibility: low, Demand: high, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LHH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_LHH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MLL 场景
echo "=== 运行 MLL 场景 ==="
echo "Flexibility: mid, Demand: low, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MLL_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MLL_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MLM 场景
echo "=== 运行 MLM 场景 ==="
echo "Flexibility: mid, Demand: low, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MLM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MLM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MLH 场景
echo "=== 运行 MLH 场景 ==="
echo "Flexibility: mid, Demand: low, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MLH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MLH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MML 场景
echo "=== 运行 MML 场景 ==="
echo "Flexibility: mid, Demand: mid, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MML_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MML_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MMM 场景
echo "=== 运行 MMM 场景 ==="
echo "Flexibility: mid, Demand: mid, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MMM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MMM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MMH 场景
echo "=== 运行 MMH 场景 ==="
echo "Flexibility: mid, Demand: mid, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MMH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MMH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MHL 场景
echo "=== 运行 MHL 场景 ==="
echo "Flexibility: mid, Demand: high, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MHL_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MHL_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MHM 场景
echo "=== 运行 MHM 场景 ==="
echo "Flexibility: mid, Demand: high, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MHM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MHM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MHH 场景
echo "=== 运行 MHH 场景 ==="
echo "Flexibility: mid, Demand: high, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MHH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MHH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HLL 场景
echo "=== 运行 HLL 场景 ==="
echo "Flexibility: high, Demand: low, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HLL_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HLL_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HLM 场景
echo "=== 运行 HLM 场景 ==="
echo "Flexibility: high, Demand: low, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HLM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HLM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HLH 场景
echo "=== 运行 HLH 场景 ==="
echo "Flexibility: high, Demand: low, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HLH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HLH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HML 场景
echo "=== 运行 HML 场景 ==="
echo "Flexibility: high, Demand: mid, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HML_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HML_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HMM 场景
echo "=== 运行 HMM 场景 ==="
echo "Flexibility: high, Demand: mid, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HMM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HMM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HMH 场景
echo "=== 运行 HMH 场景 ==="
echo "Flexibility: high, Demand: mid, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HMH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HMH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HHL 场景
echo "=== 运行 HHL 场景 ==="
echo "Flexibility: high, Demand: high, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HHL_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HHL_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HHM 场景
echo "=== 运行 HHM 场景 ==="
echo "Flexibility: high, Demand: high, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HHM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HHM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HHH 场景
echo "=== 运行 HHH 场景 ==="
echo "Flexibility: high, Demand: high, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HHH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_HHH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NLL 场景
echo "=== 运行 NLL 场景 ==="
echo "Flexibility: non_constrained, Demand: low, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NLL_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NLL_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NLM 场景
echo "=== 运行 NLM 场景 ==="
echo "Flexibility: non_constrained, Demand: low, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NLM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NLM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NLH 场景
echo "=== 运行 NLH 场景 ==="
echo "Flexibility: non_constrained, Demand: low, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NLH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NLH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NML 场景
echo "=== 运行 NML 场景 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NML_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NML_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NMM 场景
echo "=== 运行 NMM 场景 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NMM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NMM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NMH 场景
echo "=== 运行 NMH 场景 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NMH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NMH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NHL 场景
echo "=== 运行 NHL 场景 ==="
echo "Flexibility: non_constrained, Demand: high, Market: low"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NHL_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NHL_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NHM 场景
echo "=== 运行 NHM 场景 ==="
echo "Flexibility: non_constrained, Demand: high, Market: mid"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NHM_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NHM_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NHH 场景
echo "=== 运行 NHH 场景 ==="
echo "Flexibility: non_constrained, Demand: high, Market: high"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NHH_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_NHH_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


echo "所有配置的模拟已完成！"
echo "总共运行了 $total_configs 个配置"
echo "结果文件位于: results/version-*/"
