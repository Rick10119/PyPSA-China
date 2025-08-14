#!/usr/bin/env python3
"""
生成不同电解铝厂灵活性、需求和市场机会组合的配置文件
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

def generate_value_test_configs():
    """
    生成48个配置文件：4种flexibility × 3种demand × 3种market × 1种容量设置 + 3种market × 1种non-flexible设置
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
    demand_levels = ['low', 'mid', 'high']  # L, M, H
    market_levels = ['low', 'mid', 'high']  # L, M, H
    year = 2050  # 固定为2050年
    
    # 将级别映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    config_count = 0
    
    # 生成100p容量配置（4×3×3 = 36个）
    for flex in flexibility_levels:
        for demand in demand_levels:
            for market in market_levels:
                # 生成100p容量配置
                config_100p = generate_single_config(
                    base_config, base_version, flex, demand, market, 
                    capacity_ratio=1.0, is_non_flexible=False, year=year
                )
                
                # 保存100p配置文件
                scenario_suffix = f"{flex_map[flex]}{demand_map[demand]}{market_map[market]}"
                config_filename_100p = configs_dir / f"config_{scenario_suffix}_{year}_100p.yaml"
                with open(config_filename_100p, 'w', encoding='utf-8') as f:
                    yaml.dump(config_100p, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
                
                config_count += 1
                print(f"已生成配置文件 {config_count}/48:")
                print(f"  {scenario_suffix}_{year}_100p.yaml - 版本: {config_100p['version']}")
                print()
    
    # 为non-flexible情景，每个market只生成一个配置文件（使用中等水平的flex和demand）
    print("=== 生成non-flexible配置文件 ===")
    print("注意：对于non-flexible情景，由于没有灵活性，不同的flex和demand场景会产生相同的结果")
    print("因此每个market只生成一个配置文件（使用中等水平的flex和demand）")
    print()
    
    for market in market_levels:
        # 对于non-flexible情景，使用中等水平的flex和demand
        flex = 'mid'  # 使用中等灵活性
        demand = 'mid'  # 使用中等需求
        
        # 生成non-flexible配置
        config_non_flex = generate_single_config(
            base_config, base_version, flex, demand, market, 
            capacity_ratio=1.0, is_non_flexible=True, year=year
        )
        
        # 保存non-flexible配置文件
        scenario_suffix = f"{flex_map[flex]}{demand_map[demand]}{market_map[market]}"
        config_filename_non_flex = configs_dir / f"config_{scenario_suffix}_{year}_non_flexible.yaml"
        with open(config_filename_non_flex, 'w', encoding='utf-8') as f:
            yaml.dump(config_non_flex, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
        
        config_count += 1
        print(f"已生成non-flexible配置文件 {config_count}/48:")
        print(f"  {scenario_suffix}_{year}_non_flexible.yaml - 版本: {config_non_flex['version']}")
        print(f"  (使用中等flex和demand，适用于所有non-flexible情景)")
        print()
    
    print(f"总共生成了 {config_count} 个配置文件")
    print("配置说明:")
    print("- 100p容量配置: 36个 (4种flex × 3种demand × 3种market)")
    print("- non-flexible配置: 3个 (3种market，每个使用中等flex和demand)")
    print("- 总计: 36 + 3 = 39个配置文件")

def generate_single_config(base_config, base_version, flex, demand, market, capacity_ratio, is_non_flexible, year):
    """
    生成单个配置文件
    
    Args:
        base_config: 基础配置
        base_version: 基础版本号
        flex: 灵活性级别 (low/mid/high/non_constrained)
        demand: 需求级别 (low/mid/high)
        market: 市场机会级别 (low/mid/high)
        capacity_ratio: 容量比例
        is_non_flexible: 是否为non-flexible配置
        year: 年份 (固定为2050)
    """
    # 创建新的配置副本
    new_config = base_config.copy()
    
    # 将级别映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    scenario_suffix = f"{flex_map[flex]}{demand_map[demand]}{market_map[market]}"
    
    if is_non_flexible:
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
        
        # 设置容量比例为1.0
        new_config['aluminum']['capacity_ratio'] = 1.0
    else:
        # 100p容量配置
        new_config['version'] = f"{base_version}-{scenario_suffix}-{year}-100p"
        new_config['add_aluminum'] = True
        new_config['only_other_load'] = True
        new_config['aluminum_capacity_ratio'] = capacity_ratio
        
        # 如果是non_constrained配置，设置相应的参数
        if flex == 'non_constrained':
            new_config['iterative_optimization'] = False
            new_config['aluminum_commitment'] = False
        
        # 确保电解铝配置存在
        if 'aluminum' not in new_config:
            new_config['aluminum'] = {}
        
        # 设置容量比例
        new_config['aluminum']['capacity_ratio'] = capacity_ratio
        
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
    demand_levels = ['low', 'mid', 'high']
    market_levels = ['low', 'mid', 'high']
    year = 2050  # 固定为2050年
    
    # 将级别映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    script_count = 0
    
    # 为100p配置生成运行脚本（4×3×3 = 36个）
    for flex in flexibility_levels:
        for demand in demand_levels:
            for market in market_levels:
                scenario_suffix = f"{flex_map[flex]}{demand_map[demand]}{market_map[market]}"
                
                # 创建100p运行脚本
                script_100p = create_single_run_script(scenario_suffix, year, "100p", flex, demand, market)
                script_filename_100p = sh_scripts_dir / f"run_{scenario_suffix}_{year}_100p.sh"
                with open(script_filename_100p, 'w') as f:
                    f.write(script_100p)
                os.chmod(script_filename_100p, 0o755)
                
                script_count += 1
                print(f"已生成100p运行脚本 {script_count}/36:")
                print(f"  run_{scenario_suffix}_{year}_100p.sh")
    
    print()
    print("=== 生成non-flexible运行脚本 ===")
    print("注意：对于non-flexible情景，每个market只生成一个运行脚本")
    print()
    
    # 为non-flexible配置生成运行脚本（3个，每个market一个）
    for market in market_levels:
        # 对于non-flexible情景，使用中等水平的flex和demand
        flex = 'mid'  # 使用中等灵活性
        demand = 'mid'  # 使用中等需求
        
        scenario_suffix = f"{flex_map[flex]}{demand_map[demand]}{market_map[market]}"
        
        # 创建non-flexible运行脚本
        script_non_flex = create_single_run_script(scenario_suffix, year, "non_flexible", flex, demand, market)
        script_filename_non_flex = sh_scripts_dir / f"run_{scenario_suffix}_{year}_non_flexible.sh"
        with open(script_filename_non_flex, 'w') as f:
            f.write(script_non_flex)
        os.chmod(script_filename_non_flex, 0o755)
        
        script_count += 1
        print(f"已生成non-flexible运行脚本 {script_count}/39:")
        print(f"  run_{scenario_suffix}_{year}_non_flexible.sh")
        print(f"  (适用于所有non-flexible情景，使用中等flex和demand)")
    
    print(f"总共生成了 {script_count} 个运行脚本")
    print("脚本说明:")
    print("- 100p运行脚本: 36个 (4种flex × 3种demand × 3种market)")
    print("- non-flexible运行脚本: 3个 (3种market，每个使用中等flex和demand)")
    print("- 总计: 36 + 3 = 39个运行脚本")

def create_single_run_script(scenario_suffix, year, config_type, flex, demand, market):
    """
    创建单个运行脚本
    
    Args:
        scenario_suffix: 场景标识符 (如 LMM)
        year: 年份
        config_type: 配置类型 (100p 或 non_flexible)
        flex: 灵活性级别
        demand: 需求级别
        market: 市场机会级别
    """
    config_file = f"configs/config_{scenario_suffix}_{year}_{config_type}.yaml"
    
    if config_type == "100p":
        description = f"100%容量比例的模拟 (Flexibility: {flex}, Demand: {demand}, Market: {market}, Year: {year})"
        config_desc = "100p配置"
    else:
        description = f"Non-flexible场景的模拟 (Flexibility: {flex}, Demand: {demand}, Market: {market}, Year: {year})"
        config_desc = "non-flexible配置"
    
    script_content = f"""#!/bin/bash
# {description}

echo "开始运行 {scenario_suffix}_{year}_{config_type} 的模拟..."
echo "配置文件: {config_file}"
echo "配置类型: {config_desc}"
echo "Scenario: Flexibility={flex}, Demand={demand}, Market={market}, Year={year} ({scenario_suffix})"
echo

# 使用指定的配置文件运行snakemake
snakemake --configfile {config_file} --cores 40

echo "完成 {scenario_suffix}_{year}_{config_type} 的模拟"
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
    demand_levels = ['low', 'mid', 'high']
    market_levels = ['low', 'mid', 'high']
    year = 2050  # 固定为2050年
    
    # 将级别映射为简短的标识符
    flex_map = {'low': 'L', 'mid': 'M', 'high': 'H', 'non_constrained': 'N'}
    demand_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    market_map = {'low': 'L', 'mid': 'M', 'high': 'H'}
    
    script_content = """#!/bin/bash
# 批量运行所有配置的模拟

echo "开始批量运行所有配置的模拟..."
echo "总共39个配置: 4种flexibility × 3种demand × 3种market × 1种容量设置 + 3种market × 1种non-flexible设置"
echo "注意：对于non-flexible情景，由于没有灵活性，不同的flex和demand场景会产生相同的结果"
echo "因此每个market只运行一个non-flexible配置（使用中等水平的flex和demand）"
echo

config_count=0
total_configs=39

"""
    
    # 运行100p配置（36个）
    for flex in flexibility_levels:
        for demand in demand_levels:
            for market in market_levels:
                scenario_suffix = f"{flex_map[flex]}{demand_map[demand]}{market_map[market]}"
                
                script_content += f"""
# 运行 {scenario_suffix} 场景的100p配置
echo "=== 运行 {scenario_suffix} 场景的100p配置 ==="
echo "Flexibility: {flex}, Demand: {demand}, Market: {market}, Year: {year}"
echo

# 运行100p配置
echo "--- 运行 100p 配置 ---"
./sh_scripts/run_{scenario_suffix}_{year}_100p.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

"""
    
    # 运行non-flexible配置（3个）
    script_content += """
# === 运行non-flexible配置 ===
echo "注意：对于non-flexible情景，每个market只运行一个配置（使用中等水平的flex和demand）"
echo "因为不同的flex和demand场景会产生相同的结果"
echo

"""
    
    for market in market_levels:
        # 对于non-flexible情景，使用中等水平的flex和demand
        flex = 'mid'  # 使用中等灵活性
        demand = 'mid'  # 使用中等需求
        
        scenario_suffix = f"{flex_map[flex]}{demand_map[demand]}{market_map[market]}"
        
        script_content += f"""
# 运行 {scenario_suffix} 场景的non-flexible配置
echo "=== 运行 {scenario_suffix} 场景的non-flexible配置 ==="
echo "Flexibility: {flex} (中等), Demand: {demand} (中等), Market: {market}, Year: {year}"
echo "注意：此配置适用于所有non-flexible情景，因为不同的flex和demand会产生相同的结果"
echo

# 运行non-flexible配置
echo "--- 运行 non-flexible 配置 ---"
./sh_scripts/run_{scenario_suffix}_{year}_non_flexible.sh
config_count=$((config_count + 1))
echo "进度: $config_count/$total_configs"
echo

"""
    
    script_content += """
echo "所有配置的模拟已完成！"
echo "总共运行了 $total_configs 个配置"
echo "结果文件位于: results/version-*/"
echo
echo "配置总结:"
echo "- 100p容量配置: 36个 (4种flex × 3种demand × 3种market)"
echo "- non-flexible配置: 3个 (3种market，每个使用中等flex和demand)"
echo "- 总计: 36 + 3 = 39个配置"
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
    
    # 生成39个配置文件
    print("=== 生成39个配置文件 ===")
    generate_value_test_configs()
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
    print("1. 运行单个配置: ./sh_scripts/run_LMM_2050_100p.sh")
    print("2. 运行所有配置: ./sh_scripts/run_all_configs.sh")
    print("3. 手动运行: snakemake --configfile configs/config_LMM_2050_100p.yaml --cores 40")
    print("4. 提交SLURM作业: sbatch jobs/job_LMM_2050_100p.slurm")
    print("5. 批量提交所有作业: ./submit_multiple_jobs.sh")
    print()
    print("配置说明:")
    print("- 灵活性(Flexibility): L=low, M=mid, H=high, N=non_constrained")
    print("- 需求(Demand): L=low, M=mid, H=high")
    print("- 市场机会(Market): L=low, M=mid, H=high")
    print("- 年份(Year): 固定为2050")
    print("- 容量设置: 100p(100%容量), non_flexible(基准组)")
    print("- 总配置数: 4×3×3×1 + 3×1 = 39个")
    print()
    print("优化说明:")
    print("- 对于non-flexible情景，由于没有灵活性，不同的flex和demand场景会产生相同的结果")
    print("- 因此每个market只生成一个non-flexible配置文件（使用中等水平的flex和demand）")
    print("- 这样可以减少重复计算，提高效率")

if __name__ == "__main__":
    main()
