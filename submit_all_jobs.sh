#!/bin/bash
# 批量提交所有SLURM作业

echo "=== 批量提交PyPSA-China SLURM作业 ==="
echo "开始时间: $(date)"
echo

# 检查必要的文件是否存在
if [ ! -f "job.slurm" ]; then
    echo "错误: job.slurm 文件不存在"
    exit 1
fi

if [ ! -f "run_all_capacities.sh" ]; then
    echo "错误: run_all_capacities.sh 文件不存在"
    exit 1
fi

# 确保脚本有执行权限
chmod +x run_all_capacities.sh

echo "正在提交SLURM作业..."
echo

# 提交作业并记录作业ID
JOB_ID=$(sbatch job.slurm | awk '{print $4}')

if [ $? -eq 0 ]; then
    echo "✓ 作业提交成功！"
    echo "作业ID: $JOB_ID"
    echo "作业名称: pypsa-china-all-scenarios"
    echo "配置文件: job.slurm"
    echo
    
    # 保存作业信息到文件
    echo "作业ID: $JOB_ID" > "job_info.txt"
    echo "提交时间: $(date)" >> "job_info.txt"
    echo "作业名称: pypsa-china-all-scenarios" >> "job_info.txt"
    echo "配置文件: job.slurm" >> "job_info.txt"
    
    echo "作业信息已保存到: job_info.txt"
    echo
    
    # 显示作业状态
    echo "当前作业状态:"
    squeue -j $JOB_ID
    
    echo
    echo "监控作业状态:"
    echo "  squeue -j $JOB_ID          # 查看特定作业状态"
    echo "  squeue -u $USER            # 查看用户所有作业"
    echo "  tail -f slurm_${JOB_ID}.out # 实时查看输出"
    echo "  tail -f slurm_${JOB_ID}.err # 实时查看错误"
    
else
    echo "✗ 作业提交失败！"
    exit 1
fi

echo
echo "=== 批量提交完成 ==="
echo "完成时间: $(date)" 