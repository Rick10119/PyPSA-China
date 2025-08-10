#!/usr/bin/env python3
"""
自动生成PyPSA-China的SLURM作业文件
支持不同容量比例和不包含电解铝厂的场景
"""

import os
import re
from pathlib import Path

def generate_slurm_job(scenario_name, config_file, description, output_dir="."):
    """
    生成单个SLURM作业文件
    
    Args:
        scenario_name (str): 场景名称（如 '100p', 'no_aluminum'）
        config_file (str): 配置文件名称
        description (str): 场景描述
        output_dir (str): 输出目录
    """
    
    # SLURM作业模板
    slurm_template = f"""#!/bin/bash
#SBATCH --job-name=pypsa-china-{scenario_name}        # 作业名称
#SBATCH --nodes=1                # 节点数量
#SBATCH --ntasks=1               # 总任务数
#SBATCH --cpus-per-task=40       # 每个任务的CPU核心数
#SBATCH --mem-per-cpu=15G        # 每个CPU核心的内存
#SBATCH --time=12:00:00          # 总运行时间限制 (12小时)
#SBATCH --mail-type=begin        # 作业开始时发送邮件
#SBATCH --mail-type=end          # 作业结束时发送邮件
#SBATCH --mail-type=fail         # 作业失败时发送邮件
#SBATCH --mail-user=rl8728@princeton.edu
#SBATCH --output=slurm_{scenario_name}_%j.out    # 标准输出文件
#SBATCH --error=slurm_{scenario_name}_%j.err     # 标准错误文件

# 设置日志文件
LOG_FILE="job_{scenario_name}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== PyPSA-China {description}作业开始 ==="
echo "开始时间: $(date)"
echo "作业ID: $SLURM_JOB_ID"
echo "节点: $SLURM_NODELIST"
echo "日志文件: $LOG_FILE"
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
echo "开始运行{description}模拟..."
echo "配置文件: {config_file}"
echo "开始时间: $(date)"
echo

START_TIME=$(date +%s)

if snakemake --configfile {config_file} --cores 40; then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo "✓ {description}模拟运行成功！"
    echo "运行时间: $((DURATION / 3600))小时 $((DURATION % 3600 / 60))分钟 $((DURATION % 60))秒"
else
    echo "✗ {description}模拟运行失败！"
    exit 1
fi

echo
echo "=== PyPSA-China {description}作业完成 ==="
echo "完成时间: $(date)"
echo "日志文件: $LOG_FILE"
"""
    
    # 生成文件名
    job_filename = f"job_{scenario_name}.slurm"
    job_path = Path(output_dir) / job_filename
    
    # 写入文件
    with open(job_path, 'w', encoding='utf-8') as f:
        f.write(slurm_template)
    
    # 设置执行权限
    os.chmod(job_path, 0o755)
    
    print(f"✓ 已生成: {job_filename}")
    return job_filename

def generate_all_slurm_jobs():
    """
    生成所有场景的SLURM作业文件
    """
    print("开始生成PyPSA-China的SLURM作业文件...")
    print()
    
    # 定义所有场景
    scenarios = [
        {
            "name": "no_aluminum",
            "config": "config_no_aluminum.yaml",
            "description": "不包含电解铝厂的基准场景"
        },
        {
            "name": "100p",
            "config": "config_100p.yaml",
            "description": "100%容量比例"
        },
        {
            "name": "90p",
            "config": "config_90p.yaml",
            "description": "90%容量比例"
        },
        {
            "name": "80p",
            "config": "config_80p.yaml",
            "description": "80%容量比例"
        },
        {
            "name": "70p",
            "config": "config_70p.yaml",
            "description": "70%容量比例"
        },
        {
            "name": "60p",
            "config": "config_60p.yaml",
            "description": "60%容量比例"
        },
        {
            "name": "55p",
            "config": "config_55p.yaml",
            "description": "55%容量比例"
        }
    ]
    
    generated_files = []
    
    # 生成每个场景的作业文件
    for scenario in scenarios:
        filename = generate_slurm_job(
            scenario["name"],
            scenario["config"],
            scenario["description"]
        )
        generated_files.append(filename)
    
    print()
    print("=== 生成完成 ===")
    print(f"共生成 {len(generated_files)} 个SLURM作业文件:")
    for filename in generated_files:
        print(f"  - {filename}")
    
    print()
    print("使用方法:")
    print("1. 提交单个作业: sbatch job_100p.slurm")
    print("2. 批量提交所有作业: ./submit_multiple_jobs.sh")
    print("3. 只提交容量比例作业: ./submit_capacity_jobs.sh")
    
    return generated_files

def generate_custom_slurm_job(scenario_name, config_file, description, 
                             nodes=1, cpus_per_task=40, mem_per_cpu="15G", 
                             time_limit="12:00:00", mail_user="rl8728@princeton.edu"):
    """
    生成自定义参数的SLURM作业文件
    
    Args:
        scenario_name (str): 场景名称
        config_file (str): 配置文件名称
        description (str): 场景描述
        nodes (int): 节点数量
        cpus_per_task (int): 每个任务的CPU核心数
        mem_per_cpu (str): 每个CPU核心的内存
        time_limit (str): 时间限制
        mail_user (str): 邮件用户
    """
    
    slurm_template = f"""#!/bin/bash
#SBATCH --job-name=pypsa-china-{scenario_name}        # 作业名称
#SBATCH --nodes={nodes}                # 节点数量
#SBATCH --ntasks=1               # 总任务数
#SBATCH --cpus-per-task={cpus_per_task}       # 每个任务的CPU核心数
#SBATCH --mem-per-cpu={mem_per_cpu}        # 每个CPU核心的内存
#SBATCH --time={time_limit}          # 总运行时间限制
#SBATCH --mail-type=begin        # 作业开始时发送邮件
#SBATCH --mail-type=end          # 作业结束时发送邮件
#SBATCH --mail-type=fail         # 作业失败时发送邮件
#SBATCH --mail-user={mail_user}
#SBATCH --output=slurm_{scenario_name}_%j.out    # 标准输出文件
#SBATCH --error=slurm_{scenario_name}_%j.err     # 标准错误文件

# 设置日志文件
LOG_FILE="job_{scenario_name}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== PyPSA-China {description}作业开始 ==="
echo "开始时间: $(date)"
echo "作业ID: $SLURM_JOB_ID"
echo "节点: $SLURM_NODELIST"
echo "日志文件: $LOG_FILE"
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
echo "开始运行{description}模拟..."
echo "配置文件: {config_file}"
echo "开始时间: $(date)"
echo

START_TIME=$(date +%s)

if snakemake --configfile {config_file} --cores {cpus_per_task}; then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo "✓ {description}模拟运行成功！"
    echo "运行时间: $((DURATION / 3600))小时 $((DURATION % 3600 / 60))分钟 $((DURATION % 60))秒"
else
    echo "✗ {description}模拟运行失败！"
    exit 1
fi

echo
echo "=== PyPSA-China {description}作业完成 ==="
echo "完成时间: $(date)"
echo "日志文件: $LOG_FILE"
"""
    
    # 生成文件名
    job_filename = f"job_{scenario_name}.slurm"
    
    # 写入文件
    with open(job_filename, 'w', encoding='utf-8') as f:
        f.write(slurm_template)
    
    # 设置执行权限
    os.chmod(job_filename, 0o755)
    
    print(f"✓ 已生成自定义作业文件: {job_filename}")
    return job_filename

def main():
    """
    主函数
    """
    print("PyPSA-China SLURM作业文件生成器")
    print("=" * 50)
    print()
    
    # 检查是否在正确的目录
    if not os.path.exists("config.yaml"):
        print("警告: 未找到config.yaml文件，请确保在PyPSA-China项目根目录下运行此脚本")
        print()
    
    # 生成所有标准作业文件
    generated_files = generate_all_slurm_jobs()
    
    print()
    print("可选功能:")
    print("如需生成自定义参数的作业文件，可以调用:")
    print("  generate_custom_slurm_job('custom_name', 'config.yaml', '描述', nodes=2, cpus_per_task=32)")
    
    return generated_files

if __name__ == "__main__":
    main() 