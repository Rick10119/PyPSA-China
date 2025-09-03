#!/usr/bin/env python3
"""
生成不同电解铝厂灵活性、市场机会、年份和容量比例的配置文件
"""

import yaml
import os
import shutil
from pathlib import Path

def clean_directories():
    """
    清空configs和sh_scripts文件夹中的所有文件
    """
    # 清理configs文件夹
    configs_dir = Path('configs')
    if configs_dir.exists():
        # 删除所有文件
        for file_path in configs_dir.glob('*'):
            if file_path.is_file():
                file_path.unlink()
                print(f"已删除: {file_path}")
        print("configs文件夹已清空")
    else:
        configs_dir.mkdir(exist_ok=True)
        print("创建configs文件夹")
    
    # 清理sh_scripts文件夹
    sh_scripts_dir = Path('sh_scripts')
    if sh_scripts_dir.exists():
        # 删除所有.sh文件
        for file_path in sh_scripts_dir.glob('*.sh'):
            if file_path.is_file():
                file_path.unlink()
                print(f"已删除: {file_path}")
        print("sh_scripts文件夹中的.sh文件已清空")
    else:
        sh_scripts_dir.mkdir(exist_ok=True)
        print("创建sh_scripts文件夹")

def generate_capacity_test_configs():
    """
    生成468个配置文件：4种flexibility × 1种demand × 3种market × 3种年份 × 13种容量设置
    """
    # 创建configs文件夹（如果不存在）
    configs_dir = Path('configs')
    configs_dir.mkdir(exist_ok=True)
    
    # 读取原始配置文件
    with open('config.yaml', 'r', encoding='utf-8') as f:
        base_config = yaml.safe_load(f)
    
    base_version = base_config['version']
    
    # 定义所有可能的组合
    flexibility_levels = ['low', 'mid', 'high', 'non_constrained']  # L, M, H, N
    demand_level = 'mid'  # 固定为M
    market_levels = ['low', 'mid', 'high']  # L, M, H
    years = [2050]  # 年份
    
    # 定义过剩产能保留比例 (Cap值)
    # Cap=10%意味着保留过剩产能的10%，对应实际容量比例约为70%
    # Cap=20%意味着保留过剩产能的20%，对应实际容量比例约为80%
    # 等等
    cap_ratios = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # 10%, 20%, ..., 100%
    
    # 将级别映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'mid': 'M'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    config_count = 0
    
    for flex in flexibility_levels:
        for market in market_levels:
            for year in years:
                scenario_suffix = f"{flex_map[flex]}{demand_map[demand_level]}{market_map[market]}"
                
                # 生成non-flexible配置
                config_non_flex = generate_single_config(
                    base_config, base_version, flex, demand_level, market, 
                    cap_ratio=None, year=year, config_type='non_flexible'
                )
                
                # 保存non-flexible配置文件
                config_filename_non_flex = configs_dir / f"config_{scenario_suffix}_{year}_non_flexible.yaml"
                
                with open(config_filename_non_flex, 'w', encoding='utf-8') as f:
                    yaml.dump(config_non_flex, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
                
                config_count += 1
                print(f"已生成配置文件 {config_count}/468:")
                print(f"  {scenario_suffix}_{year}_non_flexible.yaml - 版本: {config_non_flex['version']}")
                print(f"  Flexibility: {flex}, Demand: {demand_level}, Market: {market}, Year: {year}, Capacity: non_flexible")
                print()
                
                # 生成no aluminum配置
                config_no_al = generate_single_config(
                    base_config, base_version, flex, demand_level, market, 
                    cap_ratio=None, year=year, config_type='no_aluminum'
                )
                
                # 保存no aluminum配置文件
                config_filename_no_al = configs_dir / f"config_{scenario_suffix}_{year}_no_aluminum.yaml"
                
                with open(config_filename_no_al, 'w', encoding='utf-8') as f:
                    yaml.dump(config_no_al, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
                
                config_count += 1
                print(f"已生成配置文件 {config_count}/468:")
                print(f"  {scenario_suffix}_{year}_no_aluminum.yaml - 版本: {config_no_al['version']}")
                print(f"  Flexibility: {flex}, Demand: {demand_level}, Market: {market}, Year: {year}, Capacity: no_aluminum")
                print()
                
                # 生成各种过剩产能保留比例的配置
                for cap_ratio in cap_ratios:
                    # 计算实际容量比例
                    actual_capacity_ratio = calculate_actual_capacity_ratio(year, cap_ratio, demand_level)
                    
                    # 生成配置文件
                    config = generate_single_config(
                        base_config, base_version, flex, demand_level, market, 
                        cap_ratio, year, config_type='capacity', actual_capacity_ratio=actual_capacity_ratio
                    )
                    
                    # 保存配置文件 - 保持原有命名模式
                    cap_percentage = int(cap_ratio * 100)
                    config_filename = configs_dir / f"config_{scenario_suffix}_{year}_{cap_percentage}p.yaml"
                    
                    with open(config_filename, 'w', encoding='utf-8') as f:
                        yaml.dump(config, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
                    
                    config_count += 1
                    print(f"已生成配置文件 {config_count}/468:")
                    print(f"  {scenario_suffix}_{year}_{cap_percentage}p.yaml - 版本: {config['version']}")
                    print(f"  Flexibility: {flex}, Demand: {demand_level}, Market: {market}, Year: {year}, Cap: {cap_percentage}%, Actual: {actual_capacity_ratio:.1%}")
                    print()
    
    print(f"总共生成了 {config_count} 个配置文件")

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
    
    demand = demand_by_year.get(str(year), 2902.417177819193)
    
    # 计算实际容量比例: demand/capacity × (1-cap) + cap
    actual_ratio = (demand / total_capacity) * (1 - cap_ratio) + cap_ratio
    
    return actual_ratio

def generate_single_config(base_config, base_version, flex, demand, market, cap_ratio, year, config_type, actual_capacity_ratio=None):
    """
    生成单个配置文件
    
    Args:
        base_config: 基础配置
        base_version: 基础版本号
        flex: 灵活性级别 (low/mid/high/non_constrained)
        demand: 需求级别 (固定为mid)
        market: 市场机会级别 (low/mid/high)
        cap_ratio: 过剩产能保留比例 (0.1-1.0) 或 None
        year: 年份 (2050)
        config_type: 配置类型 ('non_flexible', 'no_aluminum', 'capacity')
        actual_capacity_ratio: 实际容量比例 (如果config_type为'capacity')
    """
    # 创建新的配置副本
    new_config = base_config.copy()
    
    # 将级别映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'mid': 'M'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[flex]}{demand_map[demand]}{market_map[market]}"
    
    if config_type == 'non_flexible':
        # non-flexible配置
        new_config['version'] = f"{base_version}-{scenario_suffix}-{year}-non_flexible"
        new_config['add_aluminum'] = False
        new_config['only_other_load'] = False
        
        # 保留电解铝配置但不启用，以便记录情景参数
        if 'aluminum' not in new_config:
            new_config['aluminum'] = {}
        
        # 设置情景参数（即使不启用电解铝功能）
        if 'current_scenario' not in new_config['aluminum']:
            new_config['aluminum']['current_scenario'] = {}
        
        new_config['aluminum']['current_scenario']['smelter_flexibility'] = flex
        new_config['aluminum']['current_scenario']['primary_demand'] = demand
        new_config['aluminum']['current_scenario']['market_opportunity'] = market
        
        # 设置容量比例为0，表示不启用电解铝功能
        new_config['aluminum']['capacity_ratio'] = 0.0
    elif config_type == 'no_aluminum':
        # no aluminum配置
        new_config['version'] = f"{base_version}-{scenario_suffix}-{year}-no_aluminum"
        new_config['add_aluminum'] = False
        new_config['only_other_load'] = True
        
        # 保留电解铝配置但不启用，以便记录情景参数
        if 'aluminum' not in new_config:
            new_config['aluminum'] = {}
        
        # 设置情景参数（即使不启用电解铝功能）
        if 'current_scenario' not in new_config['aluminum']:
            new_config['aluminum']['current_scenario'] = {}
        
        new_config['aluminum']['current_scenario']['smelter_flexibility'] = flex
        new_config['aluminum']['current_scenario']['primary_demand'] = demand
        new_config['aluminum']['current_scenario']['market_opportunity'] = market
        
        # 设置容量比例为0，表示不启用电解铝功能
        new_config['aluminum']['capacity_ratio'] = 0.0
    else:
        # 容量配置 - 保持原有命名模式
        cap_percentage = int(cap_ratio * 100)
        new_config['version'] = f"{base_version}-{scenario_suffix}-{year}-{cap_percentage}p"
        
        # 设置电解铝厂参数
        new_config['add_aluminum'] = True
        new_config['only_other_load'] = True
        new_config['aluminum_capacity_ratio'] = actual_capacity_ratio
        
        # 如果是non_constrained配置，设置相应的参数
        if flex == 'non_constrained':
            new_config['iterative_optimization'] = False
            new_config['aluminum_commitment'] = False
        
        # 更新电解铝配置
        if 'aluminum' in new_config:
            new_config['aluminum']['capacity_ratio'] = actual_capacity_ratio
            
            # 更新当前scenario设置
            if 'current_scenario' not in new_config['aluminum']:
                new_config['aluminum']['current_scenario'] = {}
            
            new_config['aluminum']['current_scenario']['smelter_flexibility'] = flex
            new_config['aluminum']['current_scenario']['primary_demand'] = demand
            new_config['aluminum']['current_scenario']['market_opportunity'] = market
    
    # 更新年份相关设置
    new_config['costs']['year'] = year
    
    # 根据年份设置对应的planning_horizon
    if 'scenario' in new_config:
        new_config['scenario']['planning_horizons'] = [year]
    
    return new_config

def create_run_scripts():
    """
    为每个配置创建运行脚本
    """
    # 创建sh_scripts文件夹（如果不存在）
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 定义所有可能的组合
    flexibility_levels = ['low', 'mid', 'high', 'non_constrained']
    demand_level = 'mid'
    market_levels = ['low', 'mid', 'high']
    years = [2050]
    cap_ratios = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    # 将级别映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'mid': 'M'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    script_count = 0
    
    for flex in flexibility_levels:
        for market in market_levels:
            for year in years:
                scenario_suffix = f"{flex_map[flex]}{demand_map[demand_level]}{market_map[market]}"
                
                # 创建non-flexible运行脚本
                script_non_flex = create_single_run_script(scenario_suffix, year, "non_flexible", flex, demand_level, market)
                script_filename_non_flex = sh_scripts_dir / f"run_{scenario_suffix}_{year}_non_flexible.sh"
                
                with open(script_filename_non_flex, 'w') as f:
                    f.write(script_non_flex)
                os.chmod(script_filename_non_flex, 0o755)
                
                script_count += 1
                print(f"已生成运行脚本 {script_count}/468:")
                print(f"  run_{scenario_suffix}_{year}_non_flexible.sh")
                
                # 创建no aluminum运行脚本
                script_no_al = create_single_run_script(scenario_suffix, year, "no_aluminum", flex, demand_level, market)
                script_filename_no_al = sh_scripts_dir / f"run_{scenario_suffix}_{year}_no_aluminum.sh"
                
                with open(script_filename_no_al, 'w') as f:
                    f.write(script_no_al)
                os.chmod(script_filename_no_al, 0o755)
                
                script_count += 1
                print(f"已生成运行脚本 {script_count}/468:")
                print(f"  run_{scenario_suffix}_{year}_no_aluminum.sh")
                
                # 创建各种过剩产能保留比例的运行脚本
                for cap_ratio in cap_ratios:
                    cap_percentage = int(cap_ratio * 100)
                    
                    # 创建运行脚本 - 保持原有命名模式
                    script_content = create_single_run_script(scenario_suffix, year, f"{cap_percentage}p", flex, demand_level, market)
                    script_filename = sh_scripts_dir / f"run_{scenario_suffix}_{year}_{cap_percentage}p.sh"
                    
                    with open(script_filename, 'w') as f:
                        f.write(script_content)
                    os.chmod(script_filename, 0o755)
                    
                    script_count += 1
                    print(f"已生成运行脚本 {script_count}/468:")
                    print(f"  run_{scenario_suffix}_{year}_{cap_percentage}p.sh")
    
    print(f"总共生成了 {script_count} 个运行脚本")

def create_single_run_script(scenario_suffix, year, capacity_type, flex, demand, market):
    """
    创建单个运行脚本
    
    Args:
        scenario_suffix: 场景标识符 (如 LMM)
        year: 年份
        capacity_type: 容量类型 (non_flexible, no_aluminum, 10p, 20p等)
        flex: 灵活性级别
        demand: 需求级别
        market: 市场机会级别
    """
    if capacity_type == "non_flexible":
        config_file = f"configs/config_{scenario_suffix}_{year}_non_flexible.yaml"
        description = f"Non-flexible场景的模拟 (Flexibility: {flex}, Demand: {demand}, Market: {market}, Year: {year})"
        config_desc = "non_flexible配置"
    elif capacity_type == "no_aluminum":
        config_file = f"configs/config_{scenario_suffix}_{year}_no_aluminum.yaml"
        description = f"No aluminum场景的模拟 (Flexibility: {flex}, Demand: {demand}, Market: {market}, Year: {year})"
        config_desc = "no_aluminum配置"
    else:
        config_file = f"configs/config_{scenario_suffix}_{year}_{capacity_type}.yaml"
        description = f"{capacity_type}过剩产能保留比例的模拟 (Flexibility: {flex}, Demand: {demand}, Market: {market}, Year: {year})"
        config_desc = f"{capacity_type}配置"
    
    script_content = f"""#!/bin/bash
# {description}

echo "开始运行 {scenario_suffix}_{year}_{capacity_type} 的模拟..."
echo "配置文件: {config_file}"
echo "配置类型: {config_desc}"
echo "Scenario: Flexibility={flex}, Demand={demand}, Market={market}, Year={year} ({scenario_suffix})"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile {config_file} --cores 40

echo "完成 {scenario_suffix}_{year}_{capacity_type} 的模拟"
"""
    
    return script_content

def create_batch_run_script():
    """
    创建批量运行所有配置的脚本
    """
    sh_scripts_dir = Path('sh_scripts')
    sh_scripts_dir.mkdir(exist_ok=True)
    
    # 定义所有可能的组合
    flexibility_levels = ['low', 'mid', 'high', 'non_constrained']
    demand_level = 'mid'
    market_levels = ['low', 'mid', 'high']
    years = [2050]
    cap_ratios = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    # 将级别映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'mid': 'M'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    script_content = """#!/bin/bash
# 批量运行所有配置的模拟

echo "开始批量运行所有配置的模拟..."
echo "总共468个配置: 4种flexibility × 1种demand × 3种market × 3种年份 × 13种容量设置"
echo

config_count=0
total_configs=468

"""
    
    for flex in flexibility_levels:
        for market in market_levels:
            for year in years:
                scenario_suffix = f"{flex_map[flex]}{demand_map[demand_level]}{market_map[market]}"
                
                script_content += f"""
# 运行 {scenario_suffix}_{year} 场景
echo "=== 运行 {scenario_suffix}_{year} 场景 ==="
echo "Flexibility: {flex}, Demand: {demand_level}, Market: {market}, Year: {year}"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_{scenario_suffix}_{year}_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

# 运行no aluminum配置
echo "--- 运行 no aluminum 配置 ---"
./sh_scripts/run_{scenario_suffix}_{year}_no_aluminum.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

"""
                
                for cap_ratio in cap_ratios:
                    cap_percentage = int(cap_ratio * 100)
                    
                    script_content += f"""
# 运行 {cap_percentage}p 配置
echo "--- 运行 {cap_percentage}p 配置 ---"
./sh_scripts/run_{scenario_suffix}_{year}_{cap_percentage}p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

"""
    
    script_content += """
echo "所有配置的模拟已完成！"
echo "总共运行了 $total_configs 个配置"
echo "结果文件位于: results/version-*/"
"""
    
    script_filename = sh_scripts_dir / "run_all_configs.sh"
    with open(script_filename, 'w') as f:
        f.write(script_content)
    
    os.chmod(script_filename, 0o755)
    print("已生成批量运行脚本: sh_scripts/run_all_configs.sh")

def main():
    """
    主函数
    """
    print("开始生成不同场景的配置文件...")
    print()
    
    # 清空configs文件夹
    print("=== 清空configs文件夹 ===")
    clean_directories()
    print()
    
    # 生成468个配置文件
    print("=== 生成468个配置文件 ===")
    generate_capacity_test_configs()
    print()
    
    # 生成运行脚本
    print("=== 生成运行脚本 ===")
    create_run_scripts()
    print()
    
    # 生成批量运行脚本
    print("=== 生成批量运行脚本 ===")
    create_batch_run_script()
    print()
    
    print("所有文件生成完成！")
    print()
    
    # 调用高级SLURM作业文件生成器
    print("=== 生成SLURM作业文件 ===")
    try:
        import subprocess
        import sys
        
        # 直接运行脚本而不是导入
        result = subprocess.run([sys.executable, "scripts/generate_slurm_jobs_advanced.py"], 
                              capture_output=True, text=True, cwd=os.getcwd())
        
        if result.returncode == 0:
            print("✓ 成功生成SLURM作业文件")
            # 提取生成的文件数量
            if "共生成" in result.stdout:
                for line in result.stdout.split('\n'):
                    if "共生成" in line and "个SLURM作业文件" in line:
                        print(line.strip())
                        break
        else:
            print(f"✗ 生成SLURM作业文件时出错: {result.stderr}")
            print("请手动运行: python scripts/generate_slurm_jobs_advanced.py")
    except Exception as e:
        print(f"✗ 生成SLURM作业文件时出错: {e}")
        print("请手动运行: python scripts/generate_slurm_jobs_advanced.py")
    
    print()
    print("使用方法:")
    print("1. 运行单个配置: ./sh_scripts/run_LMM_2030_10p.sh")
    print("2. 运行所有配置: ./sh_scripts/run_all_configs.sh")
    print("3. 手动运行: snakemake --configfile configs/config_LMM_2030_10p.yaml --cores 40")
    print("4. 提交SLURM作业: sbatch jobs/job_LMM_2030_10p.slurm")
    print("5. 批量提交所有作业: ./submit_multiple_jobs.sh")
    print()
    print("配置说明:")
    print("- 灵活性(Flexibility): L=low, M=mid, H=high, N=non_constrained")
    print("- 需求(Demand): 固定为M=mid")
    print("- 市场机会(Market): L=low, M=mid, H=high")
    print("- 年份(Year): 2050")
    print("- 容量设置: non_flexible, no_aluminum, 10p, 20p, ..., 100p (13种)")
    print("- Cap含义: 过剩产能保留比例，如10p=保留10%过剩产能")
    print("- 实际容量比例: demand/capacity × (1-cap) + cap")
    print("- 总配置数: 4×1×3×3×13 = 468个")

if __name__ == "__main__":
    main()
