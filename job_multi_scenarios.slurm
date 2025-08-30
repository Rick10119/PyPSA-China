#!/bin/bash
#SBATCH --job-name=pypsa-china-multi-scenarios  # 作业名称
#SBATCH --nodes=1                # 节点数量
#SBATCH --ntasks=1               # 总任务数
#SBATCH --cpus-per-task=40       # 每个任务的CPU核心数
#SBATCH --mem-per-cpu=15G        # 每个CPU核心的内存
#SBATCH --time=48:00:00          # 总运行时间限制 (48小时)
#SBATCH --mail-type=begin        # 作业开始时发送邮件
#SBATCH --mail-type=end          # 作业结束时发送邮件
#SBATCH --mail-type=fail         # 作业失败时发送邮件
#SBATCH --mail-user=rl8728@princeton.edu
#SBATCH --output=slurm_%j.out    # 标准输出文件
#SBATCH --error=slurm_%j.err     # 标准错误文件

# 设置日志文件
LOG_FILE="multi_scenarios_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== PyPSA-China 多场景批量运行作业开始 ==="
echo "开始时间: $(date)"
echo "作业ID: $SLURM_JOB_ID"
echo "节点: $SLURM_NODELIST"
echo "日志文件: $LOG_FILE"
echo "输出文件: slurm_${SLURM_JOB_ID}.out"
echo "错误文件: slurm_${SLURM_JOB_ID}.err"
echo "注意: 使用 --rerun-incomplete --ignore-incomplete 参数忽略配置文件变化"
echo

# 加载必要的模块
echo "正在加载模块..."
module purge
module load anaconda3/2024.10
conda activate pypsa-plot
module load gurobi/12.0.0

echo "模块加载完成"
echo "Python环境: $(which python)"
echo "Snakemake版本: $(snakemake --version)"
echo

# 定义场景列表和配置
declare -A SCENARIOS=(
    ["no_aluminum"]="configs/config_no_aluminum.yaml"
    ["100p"]="configs/config_100p.yaml"
    ["90p"]="configs/config_90p.yaml"
    ["80p"]="configs/config_80p.yaml"
    ["70p"]="configs/config_70p.yaml"
    ["60p"]="configs/config_60p.yaml"
    ["55p"]="configs/config_55p.yaml"
)

# 记录开始时间
START_TIME=$(date +%s)
SUCCESS_COUNT=0
FAILED_COUNT=0

echo "准备运行以下场景:"
for scenario_name in "${!SCENARIOS[@]}"; do
    echo "  - $scenario_name: ${SCENARIOS[$scenario_name]}"
done
echo

# 循环运行每个场景
for scenario_name in "${!SCENARIOS[@]}"; do
    config_file="${SCENARIOS[$scenario_name]}"
    
    echo "=========================================="
    echo "开始运行场景: $scenario_name"
    echo "配置文件: $config_file"
    echo "开始时间: $(date)"
    echo "=========================================="
    
    # 检查配置文件是否存在
    if [ ! -f "$config_file" ]; then
        echo "错误: 配置文件 $config_file 不存在，跳过此场景"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        continue
    fi
    
    # 检查配置文件语法
    if ! python -c "import yaml; yaml.safe_load(open('$config_file'))" 2>/dev/null; then
        echo "错误: 配置文件 $config_file 语法错误，跳过此场景"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        continue
    fi
    
    # 运行snakemake
    SCENARIO_START=$(date +%s)
    
    echo "正在运行: snakemake --configfile $config_file --cores 40 --rerun-incomplete --ignore-incomplete"
    
    if snakemake --configfile "$config_file" --cores 40 --rerun-incomplete --ignore-incomplete; then
        SCENARIO_END=$(date +%s)
        SCENARIO_DURATION=$((SCENARIO_END - SCENARIO_START))
        echo "场景 $scenario_name 运行成功！"
        echo "运行时间: $((SCENARIO_DURATION / 3600))小时 $((SCENARIO_DURATION % 3600 / 60))分钟 $((SCENARIO_DURATION % 60))秒"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "错误: 场景 $scenario_name 运行失败！"
        echo "继续运行下一个场景..."
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
    
    echo "=========================================="
    echo "场景 $scenario_name 完成"
    echo "完成时间: $(date)"
    echo "=========================================="
    echo
    
    # 在每个场景之间稍作休息，避免资源冲突
    sleep 10
done

# 记录总运行时间
END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))

echo "=== PyPSA-China 多场景批量运行作业完成 ==="
echo "完成时间: $(date)"
echo "总运行时间: $((TOTAL_DURATION / 3600))小时 $((TOTAL_DURATION % 3600 / 60))分钟 $((TOTAL_DURATION % 60))秒"
echo "成功场景数: $SUCCESS_COUNT"
echo "失败场景数: $FAILED_COUNT"
echo "日志文件: $LOG_FILE"
echo "输出文件: slurm_${SLURM_JOB_ID}.out"
echo "错误文件: slurm_${SLURM_JOB_ID}.err"

# 生成运行摘要
echo "=== 运行摘要 ===" > "run_summary_$(date +%Y%m%d_%H%M%S).txt"
echo "作业ID: $SLURM_JOB_ID" >> "run_summary_$(date +%Y%m%d_%H%M%S).txt"
echo "开始时间: $(date -d @$START_TIME)" >> "run_summary_$(date +%Y%m%d_%H%M%S).txt"
echo "完成时间: $(date -d @$END_TIME)" >> "run_summary_$(date +%Y%m%d_%H%M%S).txt"
echo "总运行时间: $((TOTAL_DURATION / 3600))小时 $((TOTAL_DURATION % 3600 / 60))分钟 $((TOTAL_DURATION % 60))秒" >> "run_summary_$(date +%Y%m%d_%H%M%S).txt"
echo "成功场景数: $SUCCESS_COUNT" >> "run_summary_$(date +%Y%m%d_%H%M%S).txt"
echo "失败场景数: $FAILED_COUNT" >> "run_summary_$(date +%Y%m%d_%H%M%S).txt"

echo "运行摘要已保存到: run_summary_$(date +%Y%m%d_%H%M%S).txt"
echo "所有场景运行完成！" 