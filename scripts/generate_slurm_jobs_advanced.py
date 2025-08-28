#!/usr/bin/env python3
"""
高级SLURM作业文件生成器
支持从配置文件自动读取场景信息，并生成对应的SLURM作业文件
"""

import os
import yaml
import glob
from pathlib import Path
from typing import List, Dict, Any

class SlurmJobGenerator:
    """SLURM作业文件生成器"""
    
    def __init__(self, base_config="config.yaml", output_dir="jobs"):
        """
        初始化生成器
        
        Args:
            base_config (str): 基础配置文件路径
            output_dir (str): 输出目录
        """
        self.base_config = base_config
        self.output_dir = Path(output_dir)
        
        # 清理并重新创建jobs文件夹
        if self.output_dir.exists():
            # 删除所有.slurm文件
            for file_path in self.output_dir.glob('*.slurm'):
                if file_path.is_file():
                    file_path.unlink()
                    print(f"已删除旧的SLURM文件: {file_path}")
            print("jobs文件夹中的旧SLURM文件已清空")
        else:
            self.output_dir.mkdir(exist_ok=True)
            print("创建jobs文件夹")
        
        self.scenarios = []
        
    def discover_scenarios(self):
        """自动发现所有可用的场景配置文件"""
        print("正在发现场景配置文件...")
        
        # 查找configs文件夹中的所有config_*.yaml文件
        configs_dir = Path("configs")
        if not configs_dir.exists():
            print("  ✗ configs文件夹不存在，请先运行 generate_capacity_configs.py 生成配置文件")
            return []
        
        config_files = list(configs_dir.glob("config_*.yaml"))
        config_files.sort()
        
        discovered_scenarios = []
        
        for config_file in config_files:
            # 从文件名提取场景名称
            if config_file.name == "config.yaml":
                continue
                
            scenario_name = config_file.stem.replace("config_", "")
            
            # 读取配置文件获取更多信息
            try:
                with open(config_file, 'r', encoding='utf-8') as f:
                    config_data = yaml.safe_load(f)
                
                # 提取场景信息
                version = config_data.get('version', 'unknown')
                description = self._generate_description(scenario_name, config_data)
                
                scenario_info = {
                    "name": scenario_name,
                    "config_file": str(config_file),  # 使用完整路径
                    "description": description,
                    "version": version,
                    "config_data": config_data
                }
                
                discovered_scenarios.append(scenario_info)
                print(f"  ✓ 发现场景: {scenario_name} -> {description}")
                
            except Exception as e:
                print(f"  ✗ 读取配置文件 {config_file} 失败: {e}")
        
        self.scenarios = discovered_scenarios
        print(f"共发现 {len(discovered_scenarios)} 个场景")
        print()
        
        return discovered_scenarios
    
    def _generate_description(self, scenario_name: str, config_data: Dict[str, Any]) -> str:
        """根据场景名称和配置数据生成描述"""
        
        # 新的命名规则: config_LMM_2030_100p.yaml 或 config_LMM_2030_no_aluminum.yaml 或 config_LMM_2050_non_flexible.yaml
        if '_' in scenario_name:
            parts = scenario_name.split('_')
            if len(parts) >= 3:
                # 第一部分是flexibility+demand+market组合 (如 LMM)
                scenario_code = parts[0]
                # 第二部分是年份 (如 2030, 2050)
                year = parts[1]
                # 第三部分及之后是容量类型 (如 100p, no_aluminum, non_flexible)
                capacity_part = '_'.join(parts[2:])
                
                # 解析scenario代码
                flex_map = {'L': 'low', 'M': 'mid', 'H': 'high', 'N': 'non_constrained'}
                if len(scenario_code) == 3:
                    flex = flex_map.get(scenario_code[0], 'unknown')
                    demand = flex_map.get(scenario_code[1], 'unknown')
                    market = flex_map.get(scenario_code[2], 'unknown')
                    
                    # 解析容量类型
                    if capacity_part == 'no_aluminum':
                        return f"No aluminum场景 (Flexibility: {flex}, Demand: {demand}, Market: {market}, Year: {year})"
                    elif capacity_part == 'non_flexible':
                        return f"Non-flexible基准组 (Flexibility: {flex}, Demand: {demand}, Market: {market}, Year: {year})"
                    elif capacity_part.endswith('p'):
                        try:
                            percentage = int(capacity_part.replace('p', ''))
                            return f"{percentage}%容量比例 (Flexibility: {flex}, Demand: {demand}, Market: {market}, Year: {year})"
                        except ValueError:
                            pass
                    else:
                        return f"{capacity_part}配置 (Flexibility: {flex}, Demand: {demand}, Market: {market}, Year: {year})"
        
        # 兼容旧的命名规则
        if scenario_name == "no_aluminum":
            return "不包含电解铝厂的基准场景"
        
        # 检查是否是容量比例场景
        if scenario_name.endswith('p'):
            try:
                percentage = int(scenario_name.replace('p', ''))
                return f"{percentage}%容量比例"
            except ValueError:
                pass
        
        # 检查配置数据中的特殊设置
        if 'aluminum_capacity_ratio' in config_data:
            ratio = config_data['aluminum_capacity_ratio']
            if isinstance(ratio, (int, float)):
                percentage = int(ratio * 100)
                return f"{percentage}%容量比例"
        
        if config_data.get('add_aluminum') == False:
            return "不包含电解铝厂的场景"
        
        # 默认描述
        return f"场景: {scenario_name}"
    
    def generate_slurm_job(self, scenario: Dict[str, Any], 
                          template_params: Dict[str, Any] = None) -> str:
        """
        生成单个SLURM作业文件
        
        Args:
            scenario (Dict): 场景信息
            template_params (Dict): 模板参数覆盖
            
        Returns:
            str: 生成的文件名
        """
        
        # 默认模板参数
        default_params = {
            "nodes": 1,
            "ntasks": 1,
            "cpus_per_task": 40,
            "mem_per_cpu": "25G",
            "time_limit": "12:00:00",
            "mail_user": "rl8728@princeton.edu",
            "modules": [
                "module purge",
                "module load anaconda3/2024.10",
                "conda activate pypsa-plot",
                "module load gurobi/12.0.0"
            ]
        }
        
        # 合并用户参数
        if template_params:
            default_params.update(template_params)
        
        # 生成SLURM作业内容
        slurm_content = self._generate_slurm_content(scenario, default_params)
        
        # 写入文件
        job_filename = f"job_{scenario['name']}.slurm"
        job_path = self.output_dir / job_filename
        
        with open(job_path, 'w', encoding='utf-8') as f:
            f.write(slurm_content)
        
        # 设置执行权限
        os.chmod(job_path, 0o755)
        
        print(f"✓ 已生成: {self.output_dir}/{job_filename}")
        return job_filename
    
    def _generate_slurm_content(self, scenario: Dict[str, Any], 
                               params: Dict[str, Any]) -> str:
        """生成SLURM作业文件内容"""
        
        # 构建模块加载命令
        module_commands = "\n".join(params["modules"])
        
        slurm_content = f"""#!/bin/bash
#SBATCH --job-name=pypsa-china-{scenario['name']}        # 作业名称
#SBATCH --nodes={params['nodes']}                # 节点数量
#SBATCH --ntasks={params['ntasks']}               # 总任务数
#SBATCH --cpus-per-task={params['cpus_per_task']}       # 每个任务的CPU核心数
#SBATCH --mem-per-cpu={params['mem_per_cpu']}        # 每个CPU核心的内存
#SBATCH --time={params['time_limit']}          # 总运行时间限制
#SBATCH --mail-type=fail         # 作业失败时发送邮件
#SBATCH --mail-user={params['mail_user']}

# 设置日志文件
LOG_FILE="logs/job_{scenario['name']}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== PyPSA-China {scenario['description']}作业开始 ==="
echo "开始时间: $(date)"
echo "作业ID: $SLURM_JOB_ID"
echo "节点: $SLURM_NODELIST"
echo "日志文件: $LOG_FILE"
echo "版本: {scenario['version']}"
echo

# 加载必要的模块
echo "正在加载模块..."
{module_commands}

echo "模块加载完成"
echo

# 运行模拟
echo "开始运行{scenario['description']}模拟..."
echo "配置文件: {scenario['config_file']}"
echo "开始时间: $(date)"
echo

START_TIME=$(date +%s)

if snakemake --configfile {scenario['config_file']} -np; then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo "✓ {scenario['description']}模拟运行成功！"
    echo "运行时间: $((DURATION / 3600))小时 $((DURATION % 3600 / 60))分钟 $((DURATION % 60))秒"
else
    echo "✗ {scenario['description']}模拟运行失败！"
    exit 1
fi

echo
echo "=== PyPSA-China {scenario['description']}作业完成 ==="
echo "完成时间: $(date)"
echo "日志文件: $LOG_FILE"
"""
        
        return slurm_content
    
    def generate_all_jobs(self, template_params: Dict[str, Any] = None) -> List[str]:
        """
        生成所有场景的SLURM作业文件
        
        Args:
            template_params (Dict): 模板参数覆盖
            
        Returns:
            List[str]: 生成的文件名列表
        """
        
        if not self.scenarios:
            self.discover_scenarios()
        
        print("开始生成所有场景的SLURM作业文件...")
        print()
        
        generated_files = []
        
        for scenario in self.scenarios:
            filename = self.generate_slurm_job(scenario, template_params)
            generated_files.append(filename)
        
        print()
        print("=== 生成完成 ===")
        print(f"共生成 {len(generated_files)} 个SLURM作业文件:")
        for filename in generated_files:
            print(f"  - {filename}")
        
        return generated_files
    
    def generate_custom_job(self, scenario_name: str, config_file: str, 
                           description: str, **kwargs) -> str:
        """
        生成自定义场景的SLURM作业文件
        
        Args:
            scenario_name (str): 场景名称
            config_file (str): 配置文件名称
            description (str): 场景描述
            **kwargs: 其他参数
            
        Returns:
            str: 生成的文件名
        """
        
        custom_scenario = {
            "name": scenario_name,
            "config_file": config_file,
            "description": description,
            "version": "custom",
            "config_data": {}
        }
        
        return self.generate_slurm_job(custom_scenario, kwargs)

def main():
    """主函数"""
    print("PyPSA-China 高级SLURM作业文件生成器")
    print("=" * 60)
    print()
    
    # 检查是否在正确的目录
    if not os.path.exists("config.yaml"):
        print("警告: 未找到config.yaml文件，请确保在PyPSA-China项目根目录下运行此脚本")
        print()
    
    # 创建生成器实例
    generator = SlurmJobGenerator()
    
    # 自动发现场景
    scenarios = generator.discover_scenarios()
    
    if not scenarios:
        print("未发现任何场景配置文件，请先运行 generate_capacity_configs.py 生成配置文件")
        return
    
    # 生成所有作业文件
    generated_files = generator.generate_all_jobs()
    
    print()
    print("使用方法:")
    print("1. 提交单个作业: sbatch jobs/job_100p.slurm")
    print("2. 批量提交所有作业: ./submit_multiple_jobs.sh")
    print("3. 只提交容量比例作业: ./submit_capacity_jobs.sh")
    
    print()
    print("高级用法示例:")
    print("# 生成自定义参数的作业文件")
    print("generator = SlurmJobGenerator()")
    print("generator.generate_custom_job('test', 'config_test.yaml', '测试场景', cpus_per_task=32, time_limit='6:00:00')")
    
    return generated_files

if __name__ == "__main__":
    main() 