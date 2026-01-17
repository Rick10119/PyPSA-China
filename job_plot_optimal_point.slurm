#!/bin/bash
#SBATCH --job-name=optimal-point    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=40        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=15G         # memory per cpu-core (4G is default)
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
#SBATCH --qos=high               # high priority queue
#SBATCH --qos=urgent             # urgent priority queue (if available)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=rl8728@princeton.edu


module purge
module load anaconda3/2024.10
conda activate pypsa-china

# 创建日志目录
mkdir -p logs

echo "=========================================="
echo "开始运行最优点分析作业"
echo "作业ID: $SLURM_JOB_ID"
echo "开始时间: $(date)"
echo "=========================================="


echo "=========================================="
echo "运行最优点分析脚本..."
echo "=========================================="

# 运行最优点分析脚本
python scripts/plot_capacity_optimal_point.py

if [ $? -eq 0 ]; then
    echo "✅ 最优点分析脚本运行成功"
else
    echo "❌ 最优点分析脚本运行失败"
    exit 1
fi

echo "=========================================="
echo "最优点分析脚本运行完成！"
echo "结束时间: $(date)"
echo "输出目录: results/optimal_point_analysis"
echo "=========================================="

# 显示生成的文件
echo "生成的文件列表："
ls -la results/optimal_point_analysis/plots/ 2>/dev/null || echo "plots目录不存在"
ls -la results/optimal_point_analysis/ 2>/dev/null || echo "optimal_point_analysis目录不存在"

echo "=========================================="
echo "作业完成！"
echo "=========================================="
