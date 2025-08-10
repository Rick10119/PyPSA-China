#!/bin/bash
# 批量提交容量比例相关的SLURM作业

echo "=== 批量提交PyPSA-China容量比例SLURM作业 ==="
echo "开始时间: $(date)"
echo

# 定义要提交的容量比例作业文件
CAPACITY_JOBS=(
    "jobs/job_100p.slurm"             # 100%容量比例
    "jobs/job_90p.slurm"              # 90%容量比例
    "jobs/job_80p.slurm"              # 80%容量比例
    "jobs/job_70p.slurm"              # 70%容量比例
    "jobs/job_60p.slurm"              # 60%容量比例
    "jobs/job_55p.slurm"              # 55%容量比例
)

# 检查作业文件是否存在
echo "检查容量比例作业文件..."
for job_file in "${CAPACITY_JOBS[@]}"; do
    if [ -f "$job_file" ]; then
        echo "  ✓ $job_file"
    else
        echo "  ✗ $job_file (不存在)"
    fi
done
echo

# 批量提交作业
echo "开始批量提交容量比例作业..."
echo

declare -A JOB_IDS
SUBMITTED_COUNT=0
FAILED_COUNT=0

for job_file in "${CAPACITY_JOBS[@]}"; do
    echo "提交作业: $job_file"
    
    # 提交作业并获取作业ID
    JOB_ID=$(sbatch "$job_file" | awk '{print $4}')
    
    if [ $? -eq 0 ] && [ -n "$JOB_ID" ]; then
        echo "  ✓ 提交成功，作业ID: $JOB_ID"
        JOB_IDS["$job_file"]="$JOB_ID"
        SUBMITTED_COUNT=$((SUBMITTED_COUNT + 1))
        
        # 等待一小段时间再提交下一个作业
        sleep 2
    else
        echo "  ✗ 提交失败"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
    echo
done

# 显示提交结果
echo "=== 提交结果汇总 ==="
echo "成功提交: $SUBMITTED_COUNT 个容量比例作业"
echo "提交失败: $FAILED_COUNT 个作业"
echo

if [ $SUBMITTED_COUNT -gt 0 ]; then
    echo "成功提交的容量比例作业:"
    for job_file in "${CAPACITY_JOBS[@]}"; do
        if [[ -n "${JOB_IDS[$job_file]}" ]]; then
            echo "  $job_file -> 作业ID: ${JOB_IDS[$job_file]}"
        fi
    done
    echo
    
    # 保存作业信息到文件
    echo "容量比例作业信息 - $(date)" > "capacity_jobs_info.txt"
    echo "========================" >> "capacity_jobs_info.txt"
    for job_file in "${CAPACITY_JOBS[@]}"; do
        if [[ -n "${JOB_IDS[$job_file]}" ]]; then
            echo "$job_file -> ${JOB_IDS[$job_file]}" >> "capacity_jobs_info.txt"
        fi
    done
    
    echo "作业信息已保存到: capacity_jobs_info.txt"
    echo
    
    # 显示所有作业状态
    echo "当前所有作业状态:"
    squeue -u $USER
    
    echo
    echo "监控命令:"
    echo "  squeue -u $USER                    # 查看用户所有作业"
    echo "  watch -n 10 'squeue -u $USER'     # 每10秒监控作业状态"
    echo
    echo "查看特定作业输出:"
    for job_file in "${CAPACITY_JOBS[@]}"; do
        if [[ -n "${JOB_IDS[$job_file]}" ]]; then
            echo "  tail -f slurm_${job_file%.*}_${JOB_IDS[$job_file]}.out  # ${job_file}"
        fi
    done
fi

echo
echo "=== 容量比例作业批量提交完成 ==="
echo "完成时间: $(date)" 