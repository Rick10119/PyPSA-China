#!/bin/bash
#SBATCH --job-name=scenario-costs    # create a short name for your job
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
module load gurobi/13.0.0

# 创建日志目录
mkdir -p logs

echo "=========================================="
echo "开始运行多灵活性多市场比较分析作业"
echo "作业ID: $SLURM_JOB_ID"
echo "开始时间: $(date)"
echo "=========================================="


echo "=========================================="
echo "运行多灵活性多市场比较分析脚本..."
echo "=========================================="

# 运行多灵活性多市场比较分析脚本
python scripts/plot_capacity_MMM_2050.py

if [ $? -eq 0 ]; then
    echo "✅ 多灵活性多市场比较分析脚本运行成功"
else
    echo "❌ 多灵活性多市场比较分析脚本运行失败"
    exit 1
fi

echo "=========================================="
echo "多灵活性多市场比较分析脚本运行完成！"
echo "结束时间: $(date)"
echo "输出目录: results/flexibility_market_comparison"
echo "=========================================="

# 显示生成的文件
echo "生成的文件列表："
ls -la results/flexibility_market_analysis/plots/ 2>/dev/null || echo "plots目录不存在"
ls -la results/flexibility_market_analysis/ 2>/dev/null || echo "flexibility_market_comparison目录不存在"

echo "=========================================="
echo "作业完成！"
echo "=========================================="
