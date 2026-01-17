#!/bin/bash
#SBATCH --job-name=scenario-analysis    # create a short name for your job
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
conda activate pypsa-plot
module load gurobi/12.0.0

# 创建日志目录
mkdir -p logs

echo "=========================================="
echo "开始运行场景分析作业"
echo "作业ID: $SLURM_JOB_ID"
echo "开始时间: $(date)"
echo "=========================================="

# 创建输出目录
mkdir -p results/scenario_analysis

echo "=========================================="
echo "运行场景成本变化分析脚本..."
echo "=========================================="

# 运行场景成本变化分析脚本
python scripts/plot_value_scenario_comparison.py \
    --config config.yaml \
    --output results/scenario_analysis \
    --file-type costs

if [ $? -eq 0 ]; then
    echo "✅ 场景成本变化分析脚本运行成功"
else
    echo "❌ 场景成本变化分析脚本运行失败"
    exit 1
fi

echo "=========================================="
echo "运行场景容量变化分析脚本..."
echo "=========================================="

# 运行场景容量变化分析脚本
python scripts/plot_value_scenario_comparison.py \
    --config config.yaml \
    --output results/scenario_analysis \
    --file-type capacities

if [ $? -eq 0 ]; then
    echo "✅ 场景容量变化分析脚本运行成功"
else
    echo "❌ 场景容量变化分析脚本运行失败"
    exit 1
fi

echo "=========================================="
echo "所有场景分析脚本运行完成！"
echo "结束时间: $(date)"
echo "输出目录: results/scenario_analysis"
echo "=========================================="

# 显示生成的文件
echo "生成的文件列表："
echo "场景对比图表："
ls -la results/scenario_analysis/scenario_plots/ 2>/dev/null || echo "scenario_plots目录不存在"
echo ""
echo "摘要表格："
ls -la results/scenario_analysis/summary_tables/ 2>/dev/null || echo "summary_tables目录不存在"
echo ""
echo "完整目录结构："
ls -la results/scenario_analysis/ 2>/dev/null || echo "scenario_analysis目录不存在"

echo "=========================================="
echo "作业完成！"
echo "=========================================="
