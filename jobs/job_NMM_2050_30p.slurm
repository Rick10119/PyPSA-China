#!/bin/bash
#SBATCH --job-name=pypsa-china-NMM_2050_30p        # 作业名称
#SBATCH --nodes=1                # 节点数量
#SBATCH --ntasks=1               # 总任务数
#SBATCH --cpus-per-task=40       # 每个任务的CPU核心数
#SBATCH --mem-per-cpu=25G        # 每个CPU核心的内存
#SBATCH --time=24:00:00          # 总运行时间限制
#SBATCH --mail-type=fail         # 作业失败时发送邮件
#SBATCH --mail-user=rl8728@princeton.edu

# 设置日志文件
LOG_FILE="logs/job_NMM_2050_30p_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== PyPSA-China 30%容量比例 (Flexibility: non_constrained, Demand: mid, Market: mid, Year: 2050)作业开始 ==="
echo "开始时间: $(date)"
echo "作业ID: $SLURM_JOB_ID"
echo "节点: $SLURM_NODELIST"
echo "日志文件: $LOG_FILE"
echo "版本: 0815.1H.1-NMM-2050-30p"
echo

# 加载必要的模块
echo "正在加载模块..."
module purge
module load anaconda3/2024.10
conda activate pypsa-plot
module load gurobi/12.0.0

echo "模块加载完成"
echo

# 运行模拟
echo "开始运行30%容量比例 (Flexibility: non_constrained, Demand: mid, Market: mid, Year: 2050)模拟..."
echo "配置文件: configs/config_NMM_2050_30p.yaml"
echo "开始时间: $(date)"
echo

START_TIME=$(date +%s)

if snakemake --configfile configs/config_NMM_2050_30p.yaml --cores 40; then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo "✓ 30%容量比例 (Flexibility: non_constrained, Demand: mid, Market: mid, Year: 2050)模拟运行成功！"
    echo "运行时间: $((DURATION / 3600))小时 $((DURATION % 3600 / 60))分钟 $((DURATION % 60))秒"
else
    echo "✗ 30%容量比例 (Flexibility: non_constrained, Demand: mid, Market: mid, Year: 2050)模拟运行失败！"
    exit 1
fi

echo
echo "=== PyPSA-China 30%容量比例 (Flexibility: non_constrained, Demand: mid, Market: mid, Year: 2050)作业完成 ==="
echo "完成时间: $(date)"
echo "日志文件: $LOG_FILE"
