"""
情景维度工具函数
用于处理config.yaml中定义的情景维度参数
"""

import yaml
from typing import Dict, Any, Optional


def load_config(config_path: str = "config.yaml") -> Dict[str, Any]:
    """
    加载配置文件
    
    Args:
        config_path: 配置文件路径
        
    Returns:
        配置字典
    """
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def get_scenario_params(config: Dict[str, Any], 
                       smelter_flexibility: str = None,
                       primary_demand: str = None, 
                       grid_interaction: str = None) -> Dict[str, Any]:
    """
    获取指定情景的参数
    
    Args:
        config: 配置字典
        smelter_flexibility: 电解铝厂灵活性情景 ('low', 'mid', 'high')
        primary_demand: 原铝需求情景 ('low', 'mid', 'high')
        grid_interaction: 电网交互情景 ('low', 'mid', 'high')
        
    Returns:
        情景参数字典
    """
    # 如果没有指定，使用当前配置中的默认值
    current_scenario = config['aluminum']['current_scenario']
    
    if smelter_flexibility is None:
        smelter_flexibility = current_scenario['smelter_flexibility']
    if primary_demand is None:
        primary_demand = current_scenario['primary_demand']
    if grid_interaction is None:
        grid_interaction = current_scenario['grid_interaction']
    
    scenario_dims = config['aluminum']['scenario_dimensions']
    
    return {
        'smelter_flexibility': scenario_dims['smelter_flexibility'][smelter_flexibility],
        'primary_demand': scenario_dims['primary_demand'][primary_demand],
        'grid_interaction': scenario_dims['grid_interaction'][grid_interaction],
        'scenario_names': {
            'smelter_flexibility': smelter_flexibility,
            'primary_demand': primary_demand,
            'grid_interaction': grid_interaction
        }
    }


def get_smelter_params(config: Dict[str, Any], scenario: str = None) -> Dict[str, Any]:
    """
    获取电解铝厂参数
    
    Args:
        config: 配置字典
        scenario: 情景 ('low', 'mid', 'high')
        
    Returns:
        电解铝厂参数字典
    """
    if scenario is None:
        scenario = config['aluminum']['current_scenario']['smelter_flexibility']
    
    return config['aluminum']['scenario_dimensions']['smelter_flexibility'][scenario]


def get_demand_params(config: Dict[str, Any], scenario: str = None) -> Dict[str, Any]:
    """
    获取原铝需求参数
    
    Args:
        config: 配置字典
        scenario: 情景 ('low', 'mid', 'high')
        
    Returns:
        需求参数字典
    """
    if scenario is None:
        scenario = config['aluminum']['current_scenario']['primary_demand']
    
    return config['aluminum']['scenario_dimensions']['primary_demand'][scenario]


def get_grid_interaction_params(config: Dict[str, Any], scenario: str = None) -> Dict[str, Any]:
    """
    获取电网交互参数
    
    Args:
        config: 配置字典
        scenario: 情景 ('low', 'mid', 'high')
        
    Returns:
        电网交互参数字典
    """
    if scenario is None:
        scenario = config['aluminum']['current_scenario']['grid_interaction']
    
    return config['aluminum']['scenario_dimensions']['grid_interaction'][scenario]


def list_all_scenarios() -> Dict[str, list]:
    """
    列出所有可用的情景
    
    Returns:
        包含所有情景的字典
    """
    return {
        'smelter_flexibility': ['low', 'mid', 'high'],
        'primary_demand': ['low', 'mid', 'high'],
        'grid_interaction': ['low', 'mid', 'high']
    }


def generate_scenario_combinations() -> list:
    """
    生成所有可能的情景组合
    
    Returns:
        情景组合列表
    """
    scenarios = list_all_scenarios()
    combinations = []
    
    for sf in scenarios['smelter_flexibility']:
        for pd in scenarios['primary_demand']:
            for gi in scenarios['grid_interaction']:
                combinations.append({
                    'smelter_flexibility': sf,
                    'primary_demand': pd,
                    'grid_interaction': gi,
                    'name': f"{sf}-{pd}-{gi}"
                })
    
    return combinations


def print_scenario_summary(config: Dict[str, Any], scenario_params: Dict[str, Any]):
    """
    打印情景摘要
    
    Args:
        config: 配置字典
        scenario_params: 情景参数字典
    """
    print("=== 情景参数摘要 ===")
    print(f"电解铝厂灵活性: {scenario_params['scenario_names']['smelter_flexibility']}")
    print(f"原铝需求: {scenario_params['scenario_names']['primary_demand']}")
    print(f"电网交互: {scenario_params['scenario_names']['grid_interaction']}")
    print()
    
    print("电解铝厂参数:")
    sf_params = scenario_params['smelter_flexibility']
    print(f"  过载率: {sf_params['overcapacity_rate']} p.u.")
    print(f"  最小功率: {sf_params['p_min']} p.u.")
    print(f"  允许重启: {sf_params['allow_restart']}")
    print(f"  重启成本: ${sf_params['restart_costs']:,}/MW")
    print()
    
    print("需求参数:")
    pd_params = scenario_params['primary_demand']
    print(f"  国内需求比例: {pd_params['domestic_demand_ratio']*100}%")
    print(f"  出口率: {pd_params['export_rate']*100}%")
    print(f"  回收率: {pd_params['recycling_rate']*100}%")
    print(f"  产品寿命: {pd_params['product_lifetime']} 年")
    print()
    
    print("电网交互参数:")
    gi_params = scenario_params['grid_interaction']
    print(f"  VRE成本降低: {gi_params['vre_cost_reduction']*100}%")
    print(f"  电池成本降低: {gi_params['battery_cost_reduction']*100}%")
    print(f"  H2存储成本降低: {gi_params['h2_storage_cost_reduction']*100}%")
    print(f"  其他灵活需求: {gi_params['other_flexible_demand']*100}%")


if __name__ == "__main__":
    # 测试函数
    config = load_config()
    
    # 获取默认情景参数
    default_params = get_scenario_params(config)
    print_scenario_summary(config, default_params)
    
    # 生成所有组合
    combinations = generate_scenario_combinations()
    print(f"\n总共有 {len(combinations)} 种情景组合:")
    for combo in combinations[:5]:  # 只显示前5个
        print(f"  {combo['name']}")
    if len(combinations) > 5:
        print(f"  ... 还有 {len(combinations)-5} 种组合")