#!/usr/bin/env python3
"""
生成不同电解铝厂容量比例的配置文件
"""

import yaml
import os
import shutil
from pathlib import Path

def generate_capacity_configs():
    """
    生成不同容量比例的配置文件
    """
    # 容量比例设置
    capacity_ratios = [1.0, 0.9, 0.8, 0.7, 0.6, 0.55]  # 100%, 90%, 80%, 70%, 60%
    
    # 读取原始配置文件
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    base_version = config['version']
    
    for ratio in capacity_ratios:
        # 创建新的配置副本
        new_config = config.copy()
        
        # 更新版本号，添加容量比例标识
        percentage = int(ratio * 100)
        new_config['version'] = f"{base_version}-{percentage}p"
        
        # 更新容量比例设置
        new_config['aluminum_capacity_ratio'] = ratio
        if 'aluminum' in new_config:
            new_config['aluminum']['capacity_ratio'] = ratio
        
        # 生成新的配置文件名
        config_filename = f"config_{percentage}p.yaml"
        
        # 保存新的配置文件
        with open(config_filename, 'w', encoding='utf-8') as f:
            yaml.dump(new_config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
        
        print(f"已生成配置文件: {config_filename}")
        print(f"  版本号: {new_config['version']}")
        print(f"  容量比例: {percentage}%")
        print()

def create_run_scripts():
    """
    为每个容量比例创建运行脚本
    """
    capacity_ratios = [1.0, 0.9, 0.8, 0.7, 0.6, 0.55]
    
    for ratio in capacity_ratios:
        percentage = int(ratio * 100)
        config_file = f"config_{percentage}p.yaml"
        
        # 创建运行脚本
        script_content = f"""#!/bin/bash
# 运行 {percentage}% 容量比例的模拟

echo "开始运行 {percentage}% 容量比例的模拟..."
echo "配置文件: {config_file}"
echo "版本号: 0723.8H.4-{percentage}p"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile {config_file} --cores 40

echo "完成 {percentage}% 容量比例的模拟"
"""
        
        script_filename = f"run_{percentage}p.sh"
        with open(script_filename, 'w') as f:
            f.write(script_content)
        
        # 设置执行权限
        os.chmod(script_filename, 0o755)
        
        print(f"已生成运行脚本: {script_filename}")

def create_batch_run_script():
    """
    创建批量运行所有容量比例的脚本
    """
    script_content = """#!/bin/bash
# 批量运行所有容量比例的模拟

echo "开始批量运行所有容量比例的模拟..."
echo

# 运行100%容量比例
echo "=== 运行 100% 容量比例 ==="
./run_100p.sh
echo

# 运行90%容量比例
echo "=== 运行 90% 容量比例 ==="
./run_90p.sh
echo

# 运行80%容量比例
echo "=== 运行 80% 容量比例 ==="
./run_80p.sh
echo

# 运行70%容量比例
echo "=== 运行 70% 容量比例 ==="
./run_70p.sh
echo

# 运行60%容量比例
echo "=== 运行 60% 容量比例 ==="
./run_60p.sh
echo

echo "所有容量比例的模拟已完成！"
echo "结果文件位于: results/version-0723.8H.4-*p/"
"""
    
    with open("run_all_capacities.sh", 'w') as f:
        f.write(script_content)
    
    os.chmod("run_all_capacities.sh", 0o755)
    print("已生成批量运行脚本: run_all_capacities.sh")

def main():
    """
    主函数
    """
    print("开始生成不同容量比例的配置文件...")
    print()
    
    # 生成配置文件
    generate_capacity_configs()
    
    # 生成运行脚本
    create_run_scripts()
    
    # 生成批量运行脚本
    create_batch_run_script()
    
    print("所有文件生成完成！")
    print()
    print("使用方法:")
    print("1. 运行单个容量比例: ./run_100p.sh")
    print("2. 运行所有容量比例: ./run_all_capacities.sh")
    print("3. 手动运行: snakemake --configfile config_100p.yaml --cores 40")

if __name__ == "__main__":
    main() 