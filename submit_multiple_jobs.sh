#!/bin/bash
# 批量提交多个独立的SLURM作业

echo "=== 批量提交多个PyPSA-China SLURM作业 ==="
echo "开始时间: $(date)"
echo

# 定义要提交的作业文件（按优先级排序）
JOBS=(
    "jobs/job_no_aluminum.slurm"      # 不包含电解铝厂的基准场景
    "jobs/job_100p.slurm"             # 100%容量比例
    "jobs/job_90p.slurm"              # 90%容量比例
    "jobs/job_80p.slurm"              # 80%容量比例
    "jobs/job_70p.slurm"              # 70%容量比例
    "jobs/job_60p.slurm"              # 60%容量比例
    "jobs/job_55p.slurm"              # 55%容量比例
)

# 检查作业文件是否存在
echo "检查作业文件..."
for job_file in "${JOBS[@]}"; do
    if [ -f "$job_file" ]; then
        echo "  ✓ $job_file"
    else
        echo "  ✗ $job_file (不存在)"
    fi
done
echo

# 创建作业文件（如果不存在）
echo "创建缺失的作业文件..."
for job_file in "${JOBS[@]}"; do
    if [ ! -f "$job_file" ]; then
        echo "创建 $job_file..."
        create_job_file "$job_file"
    fi
done
echo

# 批量提交作业
echo "开始批量提交作业..."
echo

declare -A JOB_IDS
SUBMITTED_COUNT=0
FAILED_COUNT=0

for job_file in "${JOBS[@]}"; do
    echo "提交作业: $job_file"
    
    # 提交作业并获取作业ID
    JOB_ID=$(sbatch "$job_file" | awk '{print $4}')
    
    if [ $? -eq 0 ] && [ -n "$JOB_ID" ]; then
        echo "  ✓ 提交成功，作业ID: $JOB_ID"
        JOB_IDS["$job_file"]="$JOB_ID"
        SUBMITTED_COUNT=$((SUBMITTED_COUNT + 1))
        
        # 等待一小段时间再提交下一个作业，避免系统负载过高
        sleep 2
    else
        echo "  ✗ 提交失败"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
    echo
done

# 显示提交结果
echo "=== 提交结果汇总 ==="
echo "成功提交: $SUBMITTED_COUNT 个作业"
echo "提交失败: $FAILED_COUNT 个作业"
echo

if [ $SUBMITTED_COUNT -gt 0 ]; then
    echo "成功提交的作业:"
    for job_file in "${JOBS[@]}"; do
        if [[ -n "${JOB_IDS[$job_file]}" ]]; then
            echo "  $job_file -> 作业ID: ${JOB_IDS[$job_file]}"
        fi
    done
    echo
    
    # 保存作业信息到文件
    echo "批量作业信息 - $(date)" > "batch_jobs_info.txt"
    echo "========================" >> "batch_jobs_info.txt"
    echo "作业提交顺序:" >> "batch_jobs_info.txt"
    for job_file in "${JOBS[@]}"; do
        if [[ -n "${JOB_IDS[$job_file]}" ]]; then
            echo "$job_file -> ${JOB_IDS[$job_file]}" >> "batch_jobs_info.txt"
        fi
    done
    
    echo "作业信息已保存到: batch_jobs_info.txt"
    echo
    
    # 显示所有作业状态
    echo "当前所有作业状态:"
    squeue -u $USER
    
    echo
    echo "监控命令:"
    echo "  squeue -u $USER                    # 查看用户所有作业"
    echo "  scancel <作业ID>                   # 取消特定作业"
    echo "  tail -f slurm_<作业ID>.out         # 查看特定作业输出"
    echo "  tail -f slurm_<作业ID>.err         # 查看特定作业错误"
    echo
    echo "批量操作命令:"
    echo "  scancel \$(squeue -u $USER -h -o %i)  # 取消所有用户作业"
    echo "  watch -n 10 'squeue -u $USER'         # 每10秒监控作业状态"
fi

echo
echo "=== 批量提交完成 ==="
echo "完成时间: $(date)"

# 函数：创建作业文件
create_job_file() {
    local job_file="$1"
    local scenario_name=$(echo "$job_file" | sed 's/job_\(.*\)\.slurm/\1/')
    
    # 根据场景名称设置不同的参数
    if [ "$scenario_name" = "no_aluminum" ]; then
        local percentage="不包含电解铝厂"
        local config_file="config_no_aluminum.yaml"
        local description="不包含电解铝厂的基准场景"
    else
        local percentage=$(echo "$scenario_name" | sed 's/p$/%/')
        local config_file="config_${scenario_name}.yaml"
        local description="${percentage}容量比例"
    fi
    
    cat > "$job_file" << EOF
#!/bin/bash
#SBATCH --job-name=pypsa-china-${scenario_name}        # 作业名称
#SBATCH --nodes=1                # 节点数量
#SBATCH --ntasks=1               # 总任务数
#SBATCH --cpus-per-task=40       # 每个任务的CPU核心数
#SBATCH --mem-per-cpu=15G        # 每个CPU核心的内存
#SBATCH --time=12:00:00          # 总运行时间限制 (12小时)
#SBATCH --mail-type=begin        # 作业开始时发送邮件
#SBATCH --mail-type=end          # 作业结束时发送邮件
#SBATCH --mail-type=fail         # 作业失败时发送邮件
#SBATCH --mail-user=rl8728@princeton.edu
#SBATCH --output=logs/slurm_${scenario_name}_%j.out    # 标准输出文件
#SBATCH --error=logs/slurm_${scenario_name}_%j.err     # 标准错误文件

# 设置日志文件
LOG_FILE="job_${scenario_name}_\$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "\$LOG_FILE") 2>&1

echo "=== PyPSA-China ${description}作业开始 ==="
echo "开始时间: \$(date)"
echo "作业ID: \$SLURM_JOB_ID"
echo "节点: \$SLURM_NODELIST"
echo "日志文件: \$LOG_FILE"
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
echo "开始运行${description}模拟..."
echo "配置文件: ${config_file}"
echo "开始时间: \$(date)"
echo

START_TIME=\$(date +%s)

if snakemake --configfile ${config_file} --cores 40; then
    END_TIME=\$(date +%s)
    DURATION=\$((END_TIME - START_TIME))
    echo "✓ ${description}模拟运行成功！"
    echo "运行时间: \$((DURATION / 3600))小时 \$((DURATION % 3600 / 60))分钟 \$((DURATION % 60))秒"
else
    echo "✗ ${description}模拟运行失败！"
    exit 1
fi

echo
echo "=== PyPSA-China ${description}作业完成 ==="
echo "完成时间: \$(date)"
echo "日志文件: \$LOG_FILE"
EOF

    # 设置执行权限
    chmod +x "$job_file"
} 