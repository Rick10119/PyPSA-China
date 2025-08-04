"""
情景维度工具函数
用于处理config.yaml中定义的情景维度参数
"""

import yaml
from typing import Dict, Any, Optional
import json


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


def get_aluminum_smelter_operational_params(config: Dict[str, Any], 
                                          smelter_flexibility: str = None,
                                          al_smelter_p_nom: float = None) -> Dict[str, Any]:
    """
    获取电解铝厂运行参数，包括启动成本、待机成本、停机成本等
    
    Args:
        config: 配置字典
        smelter_flexibility: 电解铝厂灵活性情景 ('low', 'mid', 'high')
        al_smelter_p_nom: 电解铝厂容量 (MW)，用于计算成本
        
    Returns:
        电解铝厂运行参数字典
    """
    # 获取电解铝厂灵活性参数
    smelter_params = get_smelter_params(config, smelter_flexibility)
    
    # 基础参数
    operational_params = {
        'p_min_pu': smelter_params['p_min_pu'],
        'capital_cost': 163432.8,
        'stand_by_cost': 1.5,#$/MW/h
        'marginal_cost': 1,
    }
    
    # 如果提供了容量，计算成本参数
    if al_smelter_p_nom is not None:
        operational_params.update({
            'start_up_cost': 0.5 * smelter_params['restart_cost'] * al_smelter_p_nom,
            'shut_down_cost': 0.5 * smelter_params['restart_cost'] * al_smelter_p_nom,
            'stand_by_cost': operational_params['stand_by_cost'] * al_smelter_p_nom,
        })
    else:
        # 只返回单位成本
        operational_params.update({
            'start_up_cost_per_mw': 0.5 * smelter_params['restart_cost'],
            'shut_down_cost_per_mw': 0.5 * smelter_params['restart_cost'],
        })
    
    return operational_params


def get_aluminum_demand_for_year(config: Dict[str, Any], 
                                year: str,
                                primary_demand_scenario: str = None,
                                aluminum_demand_json_path: str = "data/aluminum_demand/aluminum_primary_demand_all_scenarios.json") -> float:
    """
    获取指定年份和情景的原铝需求
    
    Args:
        config: 配置字典
        year: 年份 (如 '2030', '2050')
        primary_demand_scenario: 原铝需求情景 ('low', 'mid', 'high')
        aluminum_demand_json_path: 铝需求JSON文件路径
        
    Returns:
        原铝需求 (吨)
    """
    if primary_demand_scenario is None:
        primary_demand_scenario = config['aluminum']['current_scenario']['primary_demand']
    
    # 加载铝需求数据
    with open(aluminum_demand_json_path, 'r', encoding='utf-8') as f:
        aluminum_demand_data = json.load(f)
    
    # 获取需求数据 (10kt)
    primary_demand_10kt = aluminum_demand_data['primary_aluminum_demand'][primary_demand_scenario][year]
    
    # 转换为吨
    primary_demand_tons = primary_demand_10kt * 10000
    
    return primary_demand_tons


def get_aluminum_load_for_network(config: Dict[str, Any],
                                 year: str,
                                 network_snapshots,
                                 nodes,
                                 production_ratio,
                                 primary_demand_scenario: str = None,
                                 aluminum_demand_json_path: str = "data/aluminum_demand/aluminum_primary_demand_all_scenarios.json") -> Dict[str, Any]:
    """
    获取电解铝负荷数据，用于网络建模
    
    Args:
        config: 配置字典
        year: 年份
        network_snapshots: 网络时间序列
        nodes: 节点列表
        production_ratio: 各省份生产比例
        primary_demand_scenario: 原铝需求情景
        aluminum_demand_json_path: 铝需求JSON文件路径
        
    Returns:
        包含铝负荷数据的字典
    """
    # 获取原铝需求
    primary_demand_tons = get_aluminum_demand_for_year(
        config, year, primary_demand_scenario, aluminum_demand_json_path
    )
    
    # 计算全国铝负荷 (MW)
    hours_per_year = 8760
    national_al_load = primary_demand_tons / hours_per_year
    
    # 创建铝负荷时间序列
    al_load_values = np.tile(
        national_al_load * production_ratio.values,
        (len(network_snapshots), 1)
    )
    
    aluminum_load = pd.DataFrame(
        data=al_load_values,
        index=network_snapshots,
        columns=production_ratio.index
    )
    
    return {
        'aluminum_load': aluminum_load,
        'national_al_load': national_al_load,
        'primary_demand_tons': primary_demand_tons
    }


def get_current_scenario_name(config: Dict[str, Any]) -> str:
    """
    获取当前情景组合的名称
    
    Args:
        config: 配置字典
        
    Returns:
        情景组合名称 (如 'mid-mid-mid')
    """
    current = config['aluminum']['current_scenario']
    return f"{current['smelter_flexibility']}-{current['primary_demand']}-{current['grid_interaction']}"


def validate_scenario_params(config: Dict[str, Any], 
                           smelter_flexibility: str = None,
                           primary_demand: str = None,
                           grid_interaction: str = None) -> bool:
    """
    验证情景参数是否有效
    
    Args:
        config: 配置字典
        smelter_flexibility: 电解铝厂灵活性情景
        primary_demand: 原铝需求情景
        grid_interaction: 电网交互情景
        
    Returns:
        参数是否有效
    """
    valid_scenarios = ['low', 'mid', 'high']
    
    if smelter_flexibility and smelter_flexibility not in valid_scenarios:
        return False
    if primary_demand and primary_demand not in valid_scenarios:
        return False
    if grid_interaction and grid_interaction not in valid_scenarios:
        return False
    
    return True


def print_scenario_summary(config: Dict[str, Any], scenario_params: Dict[str, Any]):
    """
    打印情景参数摘要
    
    Args:
        config: 配置字典
        scenario_params: 情景参数字典
    """
    print("=== 情景参数摘要 ===")
    print(f"情景组合: {get_current_scenario_name(config)}")
    
    # 电解铝厂参数
    smelter = scenario_params['smelter_flexibility']
    print(f"\n电解铝厂运行灵活性:")
    print(f"  最小功率: {smelter['p_min']}")
    print(f"  启动成本: {smelter['start_up_cost']} $/MW")
    print(f"  停机成本: {smelter['shut_down_cost']} $/MW")
    
    # 需求参数
    demand = scenario_params['primary_demand']
    print(f"\n原铝需求:")
    print(f"  国内需求比例: {demand['domestic_demand_ratio']}")
    print(f"  出口率: {demand['export_rate']}")
    print(f"  回收率: {demand['recycling_rate']}")
    print(f"  产品寿命: {demand['product_lifetime']} 年")
    
    # 电网交互参数
    grid = scenario_params['grid_interaction']
    print(f"\n电网交互市场机会:")
    print(f"  VRE成本降低: {grid['vre_cost_reduction']*100}%")
    print(f"  电池成本降低: {grid['battery_cost_reduction']*100}%")
    print(f"  H2存储成本降低: {grid['h2_storage_cost_reduction']*100}%")
    print(f"  其他灵活需求: {grid['other_flexible_demand']*100}%")


# 添加必要的导入
import numpy as np
import pandas as pd