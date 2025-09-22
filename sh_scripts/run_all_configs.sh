#!/bin/bash
# 批量运行所有配置的模拟

echo "开始批量运行所有配置的模拟..."
echo "总共39个配置: 4种flexibility × 3种demand × 3种market × 1种容量设置 + 3种market × 1种non-flexible设置"
echo "注意：对于non-flexible情景，由于没有灵活性，不同的flex和demand场景会产生相同的结果"
echo "因此每个market只运行一个non-flexible配置（使用中等水平的flex和demand）"
echo

config_count=0
total_configs=39


# 运行 LLL 场景的100p配置
echo "=== 运行 LLL 场景的100p配置 ==="
echo "Flexibility: low, Demand: low, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LLL_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LLM 场景的100p配置
echo "=== 运行 LLM 场景的100p配置 ==="
echo "Flexibility: low, Demand: low, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LLM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LLH 场景的100p配置
echo "=== 运行 LLH 场景的100p配置 ==="
echo "Flexibility: low, Demand: low, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LLH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LML 场景的100p配置
echo "=== 运行 LML 场景的100p配置 ==="
echo "Flexibility: low, Demand: mid, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LML_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LMM 场景的100p配置
echo "=== 运行 LMM 场景的100p配置 ==="
echo "Flexibility: low, Demand: mid, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LMM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LMH 场景的100p配置
echo "=== 运行 LMH 场景的100p配置 ==="
echo "Flexibility: low, Demand: mid, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LMH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LHL 场景的100p配置
echo "=== 运行 LHL 场景的100p配置 ==="
echo "Flexibility: low, Demand: high, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LHL_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LHM 场景的100p配置
echo "=== 运行 LHM 场景的100p配置 ==="
echo "Flexibility: low, Demand: high, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LHM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 LHH 场景的100p配置
echo "=== 运行 LHH 场景的100p配置 ==="
echo "Flexibility: low, Demand: high, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_LHH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MLL 场景的100p配置
echo "=== 运行 MLL 场景的100p配置 ==="
echo "Flexibility: mid, Demand: low, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MLL_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MLM 场景的100p配置
echo "=== 运行 MLM 场景的100p配置 ==="
echo "Flexibility: mid, Demand: low, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MLM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MLH 场景的100p配置
echo "=== 运行 MLH 场景的100p配置 ==="
echo "Flexibility: mid, Demand: low, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MLH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MML 场景的100p配置
echo "=== 运行 MML 场景的100p配置 ==="
echo "Flexibility: mid, Demand: mid, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MML_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MMM 场景的100p配置
echo "=== 运行 MMM 场景的100p配置 ==="
echo "Flexibility: mid, Demand: mid, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MMM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MMH 场景的100p配置
echo "=== 运行 MMH 场景的100p配置 ==="
echo "Flexibility: mid, Demand: mid, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MMH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MHL 场景的100p配置
echo "=== 运行 MHL 场景的100p配置 ==="
echo "Flexibility: mid, Demand: high, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MHL_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MHM 场景的100p配置
echo "=== 运行 MHM 场景的100p配置 ==="
echo "Flexibility: mid, Demand: high, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MHM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MHH 场景的100p配置
echo "=== 运行 MHH 场景的100p配置 ==="
echo "Flexibility: mid, Demand: high, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_MHH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HLL 场景的100p配置
echo "=== 运行 HLL 场景的100p配置 ==="
echo "Flexibility: high, Demand: low, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HLL_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HLM 场景的100p配置
echo "=== 运行 HLM 场景的100p配置 ==="
echo "Flexibility: high, Demand: low, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HLM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HLH 场景的100p配置
echo "=== 运行 HLH 场景的100p配置 ==="
echo "Flexibility: high, Demand: low, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HLH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HML 场景的100p配置
echo "=== 运行 HML 场景的100p配置 ==="
echo "Flexibility: high, Demand: mid, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HML_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HMM 场景的100p配置
echo "=== 运行 HMM 场景的100p配置 ==="
echo "Flexibility: high, Demand: mid, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HMM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HMH 场景的100p配置
echo "=== 运行 HMH 场景的100p配置 ==="
echo "Flexibility: high, Demand: mid, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HMH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HHL 场景的100p配置
echo "=== 运行 HHL 场景的100p配置 ==="
echo "Flexibility: high, Demand: high, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HHL_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HHM 场景的100p配置
echo "=== 运行 HHM 场景的100p配置 ==="
echo "Flexibility: high, Demand: high, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HHM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 HHH 场景的100p配置
echo "=== 运行 HHH 场景的100p配置 ==="
echo "Flexibility: high, Demand: high, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_HHH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NLL 场景的100p配置
echo "=== 运行 NLL 场景的100p配置 ==="
echo "Flexibility: non_constrained, Demand: low, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NLL_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NLM 场景的100p配置
echo "=== 运行 NLM 场景的100p配置 ==="
echo "Flexibility: non_constrained, Demand: low, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NLM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NLH 场景的100p配置
echo "=== 运行 NLH 场景的100p配置 ==="
echo "Flexibility: non_constrained, Demand: low, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NLH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NML 场景的100p配置
echo "=== 运行 NML 场景的100p配置 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NML_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NMM 场景的100p配置
echo "=== 运行 NMM 场景的100p配置 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NMM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NMH 场景的100p配置
echo "=== 运行 NMH 场景的100p配置 ==="
echo "Flexibility: non_constrained, Demand: mid, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NMH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NHL 场景的100p配置
echo "=== 运行 NHL 场景的100p配置 ==="
echo "Flexibility: non_constrained, Demand: high, Market: low, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NHL_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NHM 场景的100p配置
echo "=== 运行 NHM 场景的100p配置 ==="
echo "Flexibility: non_constrained, Demand: high, Market: mid, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NHM_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 NHH 场景的100p配置
echo "=== 运行 NHH 场景的100p配置 ==="
echo "Flexibility: non_constrained, Demand: high, Market: high, Year: 2040"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_NHH_2040_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# === 运行non-flexible配置 ===
echo "注意：对于non-flexible情景，每个market只运行一个配置（使用中等水平的flex和demand）"
echo "因为不同的flex和demand场景会产生相同的结果"
echo


# 运行 MML 场景的non-flexible配置
echo "=== 运行 MML 场景的non-flexible配置 ==="
echo "Flexibility: mid (中等), Demand: mid (中等), Market: low, Year: 2040"
echo "注意：此配置适用于所有non-flexible情景，因为不同的flex和demand会产生相同的结果"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MML_2040_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MMM 场景的non-flexible配置
echo "=== 运行 MMM 场景的non-flexible配置 ==="
echo "Flexibility: mid (中等), Demand: mid (中等), Market: mid, Year: 2040"
echo "注意：此配置适用于所有non-flexible情景，因为不同的flex和demand会产生相同的结果"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MMM_2040_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


# 运行 MMH 场景的non-flexible配置
echo "=== 运行 MMH 场景的non-flexible配置 ==="
echo "Flexibility: mid (中等), Demand: mid (中等), Market: high, Year: 2040"
echo "注意：此配置适用于所有non-flexible情景，因为不同的flex和demand会产生相同的结果"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_MMH_2040_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo


echo "所有配置的模拟已完成！"
echo "总共运行了 $total_configs 个配置"
echo "结果文件位于: results/version-*/"
echo
echo "配置总结:"
echo "- 100p容量配置: 36个 (4种flex × 3种demand × 3种market)"
echo "- non-flexible配置: 3个 (3种market，每个使用中等flex和demand)"
echo "- 总计: 36 + 3 = 39个配置"
