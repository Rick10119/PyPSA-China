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
    
    # 创建configs文件夹（如果不存在）
    configs_dir = Path('configs')
    configs_dir.mkdir(exist_ok=True)
    
    # 读取原始配置文件
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    base_version = config['version']
    
    # 获取当前scenario设置
    current_scenario = config.get('current_scenario', {})
    smelter_flex = current_scenario.get('smelter_flexibility', 'mid')
    primary_demand = current_scenario.get('primary_demand', 'high')
    market_opp = current_scenario.get('market_opportunity', 'mid')
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[smelter_flex]}{demand_map[primary_demand]}{market_map[market_opp]}"
    
    for ratio in capacity_ratios:
        # 创建新的配置副本
        new_config = config.copy()
        
        # 更新版本号，添加scenario和容量比例标识
        percentage = int(ratio * 100)
        new_config['version'] = f"{base_version}-{scenario_suffix}-{percentage}p"
        
        # 更新容量比例设置
        new_config['aluminum_capacity_ratio'] = ratio
        if 'aluminum' in new_config:
            new_config['aluminum']['capacity_ratio'] = ratio
        
        # 生成新的配置文件名
        config_filename = configs_dir / f"config_{percentage}p.yaml"
        
        # 保存新的配置文件
        with open(config_filename, 'w', encoding='utf-8') as f:
            yaml.dump(new_config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
        
        print(f"已生成配置文件: {config_filename}")
        print(f"  版本号: {new_config['version']}")
        print(f"  Scenario: {smelter_flex}-{primary_demand}-{market_opp} ({scenario_suffix})")
        print(f"  容量比例: {percentage}%")
        print()

def generate_no_aluminum_config():
    """
    生成不包含电解铝厂的配置文件
    """
    # 读取原始配置文件
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    base_version = config['version']
    
    # 获取当前scenario设置
    current_scenario = config.get('current_scenario', {})
    smelter_flex = current_scenario.get('smelter_flexibility', 'mid')
    primary_demand = current_scenario.get('primary_demand', 'high')
    market_opp = current_scenario.get('market_opportunity', 'mid')
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[smelter_flex]}{demand_map[primary_demand]}{market_map[market_opp]}"
    
    # 创建新的配置副本
    new_config = config.copy()
    
    # 更新版本号，添加scenario和no-aluminum标识
    new_config['version'] = f"{base_version}-{scenario_suffix}-no-aluminum"
    
    # 设置不包含电解铝厂
    new_config['add_aluminum'] = False
    new_config['only_other_load'] = True
    
    # 移除电解铝相关的配置
    if 'aluminum' in new_config:
        del new_config['aluminum']
    
    # 生成新的配置文件名
    config_filename = Path('configs') / "config_no_aluminum.yaml"
    
    # 保存新的配置文件
    with open(config_filename, 'w', encoding='utf-8') as f:
        yaml.dump(new_config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
    
    print(f"已生成配置文件: {config_filename}")
    print(f"  版本号: {new_config['version']}")
    print(f"  设置: add_aluminum=False, only_other_load=True")
    print()

def generate_non_flexible_config():
    """
    生成non-flexible场景的配置文件（不包含电解铝厂，也不包含其他负荷灵活性）
    """
    # 读取原始配置文件
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    base_version = config['version']
    
    # 获取当前scenario设置
    current_scenario = config.get('current_scenario', {})
    smelter_flex = current_scenario.get('smelter_flexibility', 'mid')
    primary_demand = current_scenario.get('primary_demand', 'high')
    market_opp = current_scenario.get('market_opportunity', 'mid')
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[smelter_flex]}{demand_map[primary_demand]}{market_map[market_opp]}"
    
    # 创建新的配置副本
    new_config = config.copy()
    
    # 更新版本号，添加scenario和non-flexible标识
    new_config['version'] = f"{base_version}-{scenario_suffix}-non-flexible"
    
    # 设置non-flexible场景
    new_config['add_aluminum'] = False
    new_config['only_other_load'] = False
    
    # 移除电解铝相关的配置
    if 'aluminum' in new_config:
        del new_config['aluminum']
    
    # 生成新的配置文件名
    config_filename = Path('configs') / "config_non_flexible.yaml"
    
    # 保存新的配置文件
    with open(config_filename, 'w', encoding='utf-8') as f:
        yaml.dump(new_config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
    
    print(f"已生成配置文件: {config_filename}")
    print(f"  版本号: {new_config['version']}")
    print(f"  设置: add_aluminum=False, only_other_load=False")
    print()

def create_run_sh_scripts(base_version):
    """
    为每个容量比例创建运行脚本
    """
    # 创建sh_scripts文件夹（如果不存在）
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 读取原始配置文件获取scenario信息
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    current_scenario = config.get('current_scenario', {})
    smelter_flex = current_scenario.get('smelter_flexibility', 'mid')
    primary_demand = current_scenario.get('primary_demand', 'high')
    market_opp = current_scenario.get('market_opportunity', 'mid')
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[smelter_flex]}{demand_map[primary_demand]}{market_map[market_opp]}"
    
    capacity_ratios = [1.0, 0.9, 0.8, 0.7, 0.6, 0.55]
    
    for ratio in capacity_ratios:
        percentage = int(ratio * 100)
        config_file = f"configs/config_{percentage}p.yaml"
        
        # 创建运行脚本
        script_content = f"""#!/bin/bash
# 运行 {percentage}% 容量比例的模拟

echo "开始运行 {percentage}% 容量比例的模拟..."
echo "配置文件: {config_file}"
echo "版本号: {base_version}-{scenario_suffix}-{percentage}p"
echo "Scenario: {smelter_flex}-{primary_demand}-{market_opp} ({scenario_suffix})"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile {config_file} --cores 40

echo "完成 {percentage}% 容量比例的模拟"
"""
        
        # 生成脚本文件名
        script_filename = sh_scripts_dir / f"run_{percentage}p.sh"
        
        # 保存脚本文件
        with open(script_filename, 'w') as f:
            f.write(script_content)
        
        # 设置执行权限
        os.chmod(script_filename, 0o755)
        
        print(f"已生成运行脚本: {script_filename}")

def create_no_aluminum_run_script(base_version):
    """
    为不包含电解铝厂的场景创建运行脚本
    """
    # 创建sh_scripts文件夹（如果不存在）
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 读取原始配置文件获取scenario信息
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    current_scenario = config.get('current_scenario', {})
    smelter_flex = current_scenario.get('smelter_flexibility', 'mid')
    primary_demand = current_scenario.get('primary_demand', 'high')
    market_opp = current_scenario.get('market_opportunity', 'mid')
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[smelter_flex]}{demand_map[primary_demand]}{market_map[market_opp]}"
    
    script_content = f"""#!/bin/bash
# 运行不包含电解铝厂的模拟

echo "开始运行不包含电解铝厂的模拟..."
echo "配置文件: configs/config_no_aluminum.yaml"
echo "版本号: {base_version}-{scenario_suffix}-no-aluminum"
echo "Scenario: {smelter_flex}-{primary_demand}-{market_opp} ({scenario_suffix})"
echo "设置: add_aluminum=False, only_other_load=True"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_no_aluminum.yaml --cores 40

echo "完成不包含电解铝厂的模拟"
"""
    
    script_filename = sh_scripts_dir / "run_no_aluminum.sh"
    with open(script_filename, 'w') as f:
        f.write(script_content)
    
    # 设置执行权限
    os.chmod(script_filename, 0o755)
    
    print(f"已生成运行脚本: {script_filename}")

def create_non_flexible_run_script(base_version):
    """
    为non-flexible场景创建运行脚本
    """
    # 创建sh_scripts文件夹（如果不存在）
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 读取原始配置文件获取scenario信息
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    current_scenario = config.get('current_scenario', {})
    smelter_flex = current_scenario.get('smelter_flexibility', 'mid')
    primary_demand = current_scenario.get('primary_demand', 'high')
    market_opp = current_scenario.get('market_opportunity', 'mid')
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[smelter_flex]}{demand_map[primary_demand]}{market_map[market_opp]}"
    
    script_content = f"""#!/bin/bash
# 运行non-flexible场景的模拟

echo "开始运行non-flexible场景的模拟..."
echo "配置文件: configs/config_non_flexible.yaml"
echo "版本号: {base_version}-{scenario_suffix}-non-flexible"
echo "Scenario: {smelter_flex}-{primary_demand}-{market_opp} ({scenario_suffix})"
echo "设置: add_aluminum=False, only_other_load=False"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile configs/config_non_flexible.yaml --cores 40

echo "完成non-flexible场景的模拟"
"""
    
    script_filename = sh_scripts_dir / "run_non_flexible.sh"
    with open(script_filename, 'w') as f:
        f.write(script_content)
    
    # 设置执行权限
    os.chmod(script_filename, 0o755)
    
    print(f"已生成运行脚本: {script_filename}")

def create_batch_run_script():
    """
    创建批量运行所有容量比例的脚本
    """
    # 创建sh_scripts文件夹（如果不存在）
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 读取原始配置文件获取scenario信息
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    current_scenario = config.get('current_scenario', {})
    smelter_flex = current_scenario.get('smelter_flexibility', 'mid')
    primary_demand = current_scenario.get('primary_demand', 'high')
    market_opp = current_scenario.get('market_opportunity', 'mid')
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[smelter_flex]}{demand_map[primary_demand]}{market_map[market_opp]}"
    
    script_content = f"""#!/bin/bash
# 批量运行所有容量比例的模拟

echo "开始批量运行所有容量比例的模拟..."
echo "Scenario: {smelter_flex}-{primary_demand}-{market_opp} ({scenario_suffix})"
echo

# 运行100%容量比例
echo "=== 运行 100% 容量比例 ==="
./sh_scripts/run_100p.sh
echo

# 运行90%容量比例
echo "=== 运行 90% 容量比例 ==="
./sh_scripts/run_90p.sh
echo

# 运行80%容量比例
echo "=== 运行 80% 容量比例 ==="
./sh_scripts/run_80p.sh
echo

# 运行70%容量比例
echo "=== 运行 70% 容量比例 ==="
./sh_scripts/run_70p.sh
echo

# 运行60%容量比例
echo "=== 运行 60% 容量比例 ==="
./sh_scripts/run_60p.sh
echo

# 运行55%容量比例
echo "=== 运行 55% 容量比例 ==="
./sh_scripts/run_55p.sh
echo

echo "所有容量比例的模拟已完成！"
echo "结果文件位于: results/version-{scenario_suffix}-*p/"
"""
    
    script_filename = sh_scripts_dir / "run_all_capacities.sh"
    with open(script_filename, 'w') as f:
        f.write(script_content)
    
    os.chmod(script_filename, 0o755)
    print("已生成批量运行脚本: sh_scripts/run_all_capacities.sh")

def create_all_scenarios_run_script():
    """
    创建运行所有场景（包括不包含电解铝厂的场景）的脚本
    """
    # 创建sh_scripts文件夹（如果不存在）
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 读取原始配置文件获取scenario信息
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    current_scenario = config.get('current_scenario', {})
    smelter_flex = current_scenario.get('smelter_flexibility', 'mid')
    primary_demand = current_scenario.get('primary_demand', 'high')
    market_opp = current_scenario.get('market_opportunity', 'mid')
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[smelter_flex]}{demand_map[primary_demand]}{market_map[market_opp]}"
    
    script_content = f"""#!/bin/bash
# 批量运行所有场景的模拟

echo "开始批量运行所有场景的模拟..."
echo "Scenario: {smelter_flex}-{primary_demand}-{market_opp} ({scenario_suffix})"
echo

# 运行不包含电解铝厂的场景
echo "=== 运行不包含电解铝厂的场景 ==="
./sh_scripts/run_no_aluminum.sh
echo

# 运行non-flexible场景
echo "=== 运行non-flexible场景 ==="
./sh_scripts/run_non_flexible.sh
echo

# 运行100%容量比例
echo "=== 运行 100% 容量比例 ==="
./sh_scripts/run_100p.sh
echo

# 运行90%容量比例
echo "=== 运行 90% 容量比例 ==="
./sh_scripts/run_90p.sh
echo

# 运行80%容量比例
echo "=== 运行 80% 容量比例 ==="
./sh_scripts/run_80p.sh
echo

# 运行70%容量比例
echo "=== 运行 70% 容量比例 ==="
./sh_scripts/run_70p.sh
echo

# 运行60%容量比例
echo "=== 运行 60% 容量比例 ==="
./sh_scripts/run_60p.sh
echo

# 运行55%容量比例
echo "=== 运行 55% 容量比例 ==="
./sh_scripts/run_55p.sh
echo

echo "所有场景的模拟已完成！"
echo "结果文件位于: results/version-{scenario_suffix}-*/"
"""
    
    script_filename = sh_scripts_dir / "run_all_scenarios.sh"
    with open(script_filename, 'w') as f:
        f.write(script_content)
    
    os.chmod(script_filename, 0o755)
    print("已生成所有场景运行脚本: sh_scripts/run_all_scenarios.sh")

def main():
    """
    主函数
    """
    print("开始生成不同场景的配置文件...")
    print()
    
    # 读取原始配置文件获取版本号
    with open('config.yaml', 'r', encoding='utf-8') as f:
        base_config = yaml.safe_load(f)
    base_version = base_config['version']
    
    # 生成容量比例配置文件
    print("=== 生成容量比例配置文件 ===")
    generate_capacity_configs()
    
    # 生成不包含电解铝厂的配置文件
    print("=== 生成不包含电解铝厂的配置文件 ===")
    generate_no_aluminum_config()

    # 生成non-flexible场景的配置文件
    print("=== 生成non-flexible场景的配置文件 ===")
    generate_non_flexible_config()
    
    # 生成运行脚本
    print("=== 生成运行脚本 ===")
    create_run_sh_scripts(base_version)
    create_no_aluminum_run_script(base_version)
    create_non_flexible_run_script(base_version)
    
    # 生成批量运行脚本
    print("=== 生成批量运行脚本 ===")
    create_batch_run_script()
    create_all_scenarios_run_script()
    
    print("所有文件生成完成！")
    print()
    print("使用方法:")
    print("1. 运行单个容量比例: ./sh_scripts/run_100p.sh")
    print("2. 运行不包含电解铝厂的场景: ./sh_scripts/run_no_aluminum.sh")
    print("3. 运行non-flexible场景: ./sh_scripts/run_non_flexible.sh")
    print("4. 运行所有容量比例: ./sh_scripts/run_all_capacities.sh")
    print("5. 运行所有场景: ./sh_scripts/run_all_scenarios.sh")
    print("6. 手动运行: snakemake --configfile configs/config_100p.yaml --cores 40")

if __name__ == "__main__":
    main()
