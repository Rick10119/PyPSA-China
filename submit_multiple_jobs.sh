#!/bin/bash
# 批量提交多个独立的SLURM作业

echo "=== 批量提交多个PyPSA-China SLURM作业 ==="
echo "开始时间: $(date)"
echo

# 检查jobs文件夹是否存在
if [ ! -d "jobs" ]; then
    echo "✗ jobs文件夹不存在，请先运行 generate_slurm_jobs_advanced.py 生成SLURM作业文件"
    exit 1
fi

# 自动发现所有可用的SLURM作业文件
echo "正在发现可用的SLURM作业文件..."
JOBS=()
for job_file in jobs/job_*.slurm; do
    if [ -f "$job_file" ]; then
        JOBS+=("$job_file")
        echo "  ✓ 发现作业文件: $(basename "$job_file")"
    fi
done

if [ ${#JOBS[@]} -eq 0 ]; then
    echo "✗ 未发现任何SLURM作业文件，请先运行 generate_slurm_jobs_advanced.py 生成作业文件"
    exit 1
fi

echo "共发现 ${#JOBS[@]} 个作业文件"
echo

# 按优先级排序作业（100p优先，然后按字母顺序）
echo "按优先级排序作业文件..."
JOBS_SORTED=()
# 先添加100p作业
for job_file in "${JOBS[@]}"; do
    if [[ "$job_file" == *"100p"* ]]; then
        JOBS_SORTED+=("$job_file")
    fi
done
# 再添加non_flexible作业
for job_file in "${JOBS[@]}"; do
    if [[ "$job_file" == *"non_flexible"* ]]; then
        JOBS_SORTED+=("$job_file")
    fi
done
# 最后添加其他作业
for job_file in "${JOBS[@]}"; do
    if [[ "$job_file" != *"100p"* ]] && [[ "$job_file" != *"non_flexible"* ]]; then
        JOBS_SORTED+=("$job_file")
    fi
done

echo "排序后的作业顺序:"
for i in "${!JOBS_SORTED[@]}"; do
    echo "  $((i+1)). $(basename "${JOBS_SORTED[$i]}")"
done
echo

# 批量提交作业
echo "开始批量提交作业..."
echo

declare -A JOB_IDS
SUBMITTED_COUNT=0
FAILED_COUNT=0

for job_file in "${JOBS_SORTED[@]}"; do
    echo "提交作业: $(basename "$job_file")"
    
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
    for job_file in "${JOBS_SORTED[@]}"; do
        if [[ -n "${JOB_IDS[$job_file]}" ]]; then
            echo "  $(basename "$job_file") -> 作业ID: ${JOB_IDS[$job_file]}"
        fi
    done
    echo
    
    # 保存作业信息到文件
    echo "批量作业信息 - $(date)" > "batch_jobs_info.txt"
    echo "========================" >> "batch_jobs_info.txt"
    echo "作业提交顺序:" >> "batch_jobs_info.txt"
    for job_file in "${JOBS_SORTED[@]}"; do
        if [[ -n "${JOB_IDS[$job_file]}" ]]; then
            echo "$(basename "$job_file") -> ${JOB_IDS[$job_file]}" >> "batch_jobs_info.txt"
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