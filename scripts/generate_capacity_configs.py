#!/usr/bin/env python3
"""
生成不同电解铝厂容量比例的配置文件
固定为MMMM场景，使用过剩产能保留比例计算方法
"""

import yaml
import os
import shutil
from pathlib import Path

def calculate_actual_capacity_ratio(year: int, cap_ratio: float, demand_level: str) -> float:
    """
    计算实际容量比例
    
    Args:
        year: 年份
        cap_ratio: 过剩产能保留比例 (如0.1表示10%)
        demand_level: 需求级别 ('mid')
        
    Returns:
        实际容量比例
    """
    # 总产能
    total_capacity = 4500
    
    # 根据年份设置需求 (使用实际数据)
    demand_by_year = {
        "2030": 2902.417177819193,
        "2040": 1508.1703393209764,
        "2050": 1166.6836345743664,
    }
    
    demand = demand_by_year.get(str(year), 0)
    
    # 计算实际容量比例: demand/capacity × (1-cap) + cap
    actual_ratio = (demand / total_capacity) * (1 - cap_ratio) + cap_ratio
    
    return actual_ratio

def generate_capacity_configs():
    """
    生成不同容量比例的配置文件
    固定为MMMM场景，使用过剩产能保留比例
    """
    # 过剩产能保留比例 (Cap值)
    cap_ratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # 10%, 20%, ..., 100%
    
    # 固定为MMMM场景
    flex = 'mid'
    demand_level = 'mid'
    market = 'mid'
    year = 2050  # 固定年份为2050
    
    # 创建configs文件夹（如果不存在）
    configs_dir = Path('configs')
    configs_dir.mkdir(exist_ok=True)
    
    # 读取原始配置文件
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    base_version = config['version']
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[flex]}{demand_map[demand_level]}{market_map[market]}"
    
    for cap_ratio in cap_ratios:
        # 计算实际容量比例
        actual_capacity_ratio = calculate_actual_capacity_ratio(year, cap_ratio, demand_level)
        
        # 创建新的配置副本
        new_config = config.copy()
        
        # 更新版本号，添加scenario和容量比例标识
        cap_percentage = int(cap_ratio * 100)
        new_config['version'] = f"{base_version}-{scenario_suffix}-{cap_percentage}p"
        
        # 更新容量比例设置
        new_config['aluminum_capacity_ratio'] = actual_capacity_ratio
        if 'aluminum' in new_config:
            new_config['aluminum']['capacity_ratio'] = actual_capacity_ratio
            
            # 更新当前scenario设置
            if 'current_scenario' not in new_config['aluminum']:
                new_config['aluminum']['current_scenario'] = {}
            
            new_config['aluminum']['current_scenario']['smelter_flexibility'] = flex
            new_config['aluminum']['current_scenario']['primary_demand'] = demand_level
            new_config['aluminum']['current_scenario']['market_opportunity'] = market
        
        # 生成新的配置文件名
        config_filename = configs_dir / f"config_{cap_percentage}p.yaml"
        
        # 保存新的配置文件
        with open(config_filename, 'w', encoding='utf-8') as f:
            yaml.dump(new_config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
        
        print(f"已生成配置文件: {config_filename}")
        print(f"  版本号: {new_config['version']}")
        print(f"  Scenario: {flex}-{demand_level}-{market} ({scenario_suffix})")
        print(f"  Cap: {cap_percentage}%, Actual: {actual_capacity_ratio:.1%}")
        print()

def generate_no_aluminum_config():
    """
    生成不包含电解铝厂的配置文件
    固定为MMMM场景
    """
    # 固定为MMMM场景
    flex = 'mid'
    demand_level = 'mid'
    market = 'mid'
    year = 2050
    
    # 读取原始配置文件
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    base_version = config['version']
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[flex]}{demand_map[demand_level]}{market_map[market]}"
    
    # 创建新的配置副本
    new_config = config.copy()
    
    # 更新版本号，添加scenario和no-aluminum标识
    new_config['version'] = f"{base_version}-{scenario_suffix}-no-aluminum"
    
    # 设置不包含电解铝厂
    new_config['add_aluminum'] = False
    new_config['only_other_load'] = True
    
    # 保留电解铝配置但不启用，以便记录情景参数
    if 'aluminum' not in new_config:
        new_config['aluminum'] = {}
    
    # 设置情景参数（即使不启用电解铝功能）
    if 'current_scenario' not in new_config['aluminum']:
        new_config['aluminum']['current_scenario'] = {}
    
    new_config['aluminum']['current_scenario']['smelter_flexibility'] = flex
    new_config['aluminum']['current_scenario']['primary_demand'] = demand_level
    new_config['aluminum']['current_scenario']['market_opportunity'] = market
    
    # 设置容量比例为0，表示不启用电解铝功能
    new_config['aluminum']['capacity_ratio'] = 0.0
    
    # 生成新的配置文件名
    config_filename = Path('configs') / "config_no_aluminum.yaml"
    
    # 保存新的配置文件
    with open(config_filename, 'w', encoding='utf-8') as f:
        yaml.dump(new_config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
    
    print(f"已生成配置文件: {config_filename}")
    print(f"  版本号: {new_config['version']}")
    print(f"  Scenario: {flex}-{demand_level}-{market} ({scenario_suffix})")
    print(f"  设置: add_aluminum=False, only_other_load=True")
    print()

def generate_non_flexible_config():
    """
    生成non-flexible场景的配置文件（不包含电解铝厂，也不包含其他负荷灵活性）
    固定为MMMM场景
    """
    # 固定为MMMM场景
    flex = 'mid'
    demand_level = 'mid'
    market = 'mid'
    year = 2050
    
    # 读取原始配置文件
    with open('config.yaml', 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    base_version = config['version']
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[flex]}{demand_map[demand_level]}{market_map[market]}"
    
    # 创建新的配置副本
    new_config = config.copy()
    
    # 更新版本号，添加scenario和non-flexible标识
    new_config['version'] = f"{base_version}-{scenario_suffix}-non-flexible"
    
    # 设置non-flexible场景
    new_config['add_aluminum'] = False
    new_config['only_other_load'] = False
    
    # 保留电解铝配置但不启用，以便记录情景参数
    if 'aluminum' not in new_config:
        new_config['aluminum'] = {}
    
    # 设置情景参数（即使不启用电解铝功能）
    if 'current_scenario' not in new_config['aluminum']:
        new_config['aluminum']['current_scenario'] = {}
    
    new_config['aluminum']['current_scenario']['smelter_flexibility'] = flex
    new_config['aluminum']['current_scenario']['primary_demand'] = demand_level
    new_config['aluminum']['current_scenario']['market_opportunity'] = market
    
    # 设置容量比例为0，表示不启用电解铝功能
    new_config['aluminum']['capacity_ratio'] = 0.0
    
    # 生成新的配置文件名
    config_filename = Path('configs') / "config_non_flexible.yaml"
    
    # 保存新的配置文件
    with open(config_filename, 'w', encoding='utf-8') as f:
        yaml.dump(new_config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
    
    print(f"已生成配置文件: {config_filename}")
    print(f"  版本号: {new_config['version']}")
    print(f"  Scenario: {flex}-{demand_level}-{market} ({scenario_suffix})")
    print(f"  设置: add_aluminum=False, only_other_load=False")
    print()

def create_run_sh_scripts(base_version):
    """
    为每个容量比例创建运行脚本
    固定为MMMM场景
    """
    # 创建sh_scripts文件夹（如果不存在）
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 固定为MMMM场景
    flex = 'mid'
    demand_level = 'mid'
    market = 'mid'
    year = 2050
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[flex]}{demand_map[demand_level]}{market_map[market]}"
    
    # 过剩产能保留比例 (Cap值)
    cap_ratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    for cap_ratio in cap_ratios:
        cap_percentage = int(cap_ratio * 100)
        config_file = f"configs/config_{cap_percentage}p.yaml"
        
        # 计算实际容量比例
        actual_capacity_ratio = calculate_actual_capacity_ratio(year, cap_ratio, demand_level)
        
        # 创建运行脚本
        script_content = f"""#!/bin/bash
# 运行 {cap_percentage}p 过剩产能保留比例的模拟

echo "开始运行 {cap_percentage}p 过剩产能保留比例的模拟..."
echo "配置文件: {config_file}"
echo "版本号: {base_version}-{scenario_suffix}-{cap_percentage}p"
echo "Scenario: {flex}-{demand_level}-{market} ({scenario_suffix})"
echo "Cap: {cap_percentage}%, Actual: {actual_capacity_ratio:.1%}"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile {config_file} --cores 40

echo "完成 {cap_percentage}p 过剩产能保留比例的模拟"
"""
        
        # 生成脚本文件名
        script_filename = sh_scripts_dir / f"run_{cap_percentage}p.sh"
        
        # 保存脚本文件
        with open(script_filename, 'w') as f:
            f.write(script_content)
        
        # 设置执行权限
        os.chmod(script_filename, 0o755)
        
        print(f"已生成运行脚本: {script_filename}")

def create_no_aluminum_run_script(base_version):
    """
    为不包含电解铝厂的场景创建运行脚本
    固定为MMMM场景
    """
    # 创建sh_scripts文件夹（如果不存在）
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 固定为MMMM场景
    flex = 'mid'
    demand_level = 'mid'
    market = 'mid'
    year = 2050
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[flex]}{demand_map[demand_level]}{market_map[market]}"
    
    script_content = f"""#!/bin/bash
# 运行不包含电解铝厂的模拟

echo "开始运行不包含电解铝厂的模拟..."
echo "配置文件: configs/config_no_aluminum.yaml"
echo "版本号: {base_version}-{scenario_suffix}-no-aluminum"
echo "Scenario: {flex}-{demand_level}-{market} ({scenario_suffix})"
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
    固定为MMMM场景
    """
    # 创建sh_scripts文件夹（如果不存在）
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 固定为MMMM场景
    flex = 'mid'
    demand_level = 'mid'
    market = 'mid'
    year = 2050
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[flex]}{demand_map[demand_level]}{market_map[market]}"
    
    script_content = f"""#!/bin/bash
# 运行non-flexible场景的模拟

echo "开始运行non-flexible场景的模拟..."
echo "配置文件: configs/config_non_flexible.yaml"
echo "版本号: {base_version}-{scenario_suffix}-non-flexible"
echo "Scenario: {flex}-{demand_level}-{market} ({scenario_suffix})"
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
    固定为MMMM场景
    """
    # 创建sh_scripts文件夹（如果不存在）
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 固定为MMMM场景
    flex = 'mid'
    demand_level = 'mid'
    market = 'mid'
    year = 2050
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[flex]}{demand_map[demand_level]}{market_map[market]}"
    
    script_content = f"""#!/bin/bash
# 批量运行所有容量比例的模拟

echo "开始批量运行所有容量比例的模拟..."
echo "Scenario: {flex}-{demand_level}-{market} ({scenario_suffix})"
echo

# 运行10p过剩产能保留比例
echo "=== 运行 10p 过剩产能保留比例 ==="
./sh_scripts/run_10p.sh
echo

# 运行20p过剩产能保留比例
echo "=== 运行 20p 过剩产能保留比例 ==="
./sh_scripts/run_20p.sh
echo

# 运行30p过剩产能保留比例
echo "=== 运行 30p 过剩产能保留比例 ==="
./sh_scripts/run_30p.sh
echo

# 运行40p过剩产能保留比例
echo "=== 运行 40p 过剩产能保留比例 ==="
./sh_scripts/run_40p.sh
echo

# 运行50p过剩产能保留比例
echo "=== 运行 50p 过剩产能保留比例 ==="
./sh_scripts/run_50p.sh
echo

# 运行60p过剩产能保留比例
echo "=== 运行 60p 过剩产能保留比例 ==="
./sh_scripts/run_60p.sh
echo

# 运行70p过剩产能保留比例
echo "=== 运行 70p 过剩产能保留比例 ==="
./sh_scripts/run_70p.sh
echo

# 运行80p过剩产能保留比例
echo "=== 运行 80p 过剩产能保留比例 ==="
./sh_scripts/run_80p.sh
echo

# 运行90p过剩产能保留比例
echo "=== 运行 90p 过剩产能保留比例 ==="
./sh_scripts/run_90p.sh
echo

# 运行100p过剩产能保留比例
echo "=== 运行 100p 过剩产能保留比例 ==="
./sh_scripts/run_100p.sh
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
    固定为MMMM场景
    """
    # 创建sh_scripts文件夹（如果不存在）
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 固定为MMMM场景
    flex = 'mid'
    demand_level = 'mid'
    market = 'mid'
    year = 2050
    
    # 将scenario映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[flex]}{demand_map[demand_level]}{market_map[market]}"
    
    script_content = f"""#!/bin/bash
# 批量运行所有场景的模拟

echo "开始批量运行所有场景的模拟..."
echo "Scenario: {flex}-{demand_level}-{market} ({scenario_suffix})"
echo

# 运行不包含电解铝厂的场景
echo "=== 运行不包含电解铝厂的场景 ==="
./sh_scripts/run_no_aluminum.sh
echo

# 运行non-flexible场景
echo "=== 运行non-flexible场景 ==="
./sh_scripts/run_non_flexible.sh
echo

# 运行10p过剩产能保留比例
echo "=== 运行 10p 过剩产能保留比例 ==="
./sh_scripts/run_10p.sh
echo

# 运行20p过剩产能保留比例
echo "=== 运行 20p 过剩产能保留比例 ==="
./sh_scripts/run_20p.sh
echo

# 运行30p过剩产能保留比例
echo "=== 运行 30p 过剩产能保留比例 ==="
./sh_scripts/run_30p.sh
echo

# 运行40p过剩产能保留比例
echo "=== 运行 40p 过剩产能保留比例 ==="
./sh_scripts/run_40p.sh
echo

# 运行50p过剩产能保留比例
echo "=== 运行 50p 过剩产能保留比例 ==="
./sh_scripts/run_50p.sh
echo

# 运行60p过剩产能保留比例
echo "=== 运行 60p 过剩产能保留比例 ==="
./sh_scripts/run_60p.sh
echo

# 运行70p过剩产能保留比例
echo "=== 运行 70p 过剩产能保留比例 ==="
./sh_scripts/run_70p.sh
echo

# 运行80p过剩产能保留比例
echo "=== 运行 80p 过剩产能保留比例 ==="
./sh_scripts/run_80p.sh
echo

# 运行90p过剩产能保留比例
echo "=== 运行 90p 过剩产能保留比例 ==="
./sh_scripts/run_90p.sh
echo

# 运行100p过剩产能保留比例
echo "=== 运行 100p 过剩产能保留比例 ==="
./sh_scripts/run_100p.sh
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
    
    print()
    print("使用方法:")
    print("1. 运行单个过剩产能保留比例: ./sh_scripts/run_100p.sh")
    print("2. 运行不包含电解铝厂的场景: ./sh_scripts/run_no_aluminum.sh")
    print("3. 运行non-flexible场景: ./sh_scripts/run_non_flexible.sh")
    print("4. 运行所有过剩产能保留比例: ./sh_scripts/run_all_capacities.sh")
    print("5. 运行所有场景: ./sh_scripts/run_all_scenarios.sh")
    print("6. 手动运行: snakemake --configfile configs/config_100p.yaml --cores 40")
    print()
    print("配置说明:")
    print("- 固定场景: MMMM (mid-mid-mid-mid)")
    print("- 年份: 2050")
    print("- 过剩产能保留比例: 10p, 20p, 30p, 40p, 50p, 60p, 70p, 80p, 90p, 100p")
    print("- Cap含义: 过剩产能保留比例，如10p=保留10%过剩产能")
    print("- 实际容量比例: demand/capacity × (1-cap) + cap")
    print("- 总配置数: 13个 (10个容量比例 + 2个特殊场景 + 1个non-flexible)")
    # print("7. 提交SLURM作业: sbatch jobs/job_100p.slurm")
    # print("8. 批量提交所有作业: ./submit_multiple_jobs.sh")

if __name__ == "__main__":
    main()
