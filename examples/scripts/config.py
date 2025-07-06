# config.py

# 配置参数
CONFIG = {
    "years": [2025],
    "resolution": 1,  # 8小时
    "original_cost": 64.09, # 原始成本
    "al_co2_limit": 10 * 1e6 * 1.0, # kgCO2/year, 200百万吨（2023） * 碳排放目标
    "al_demand": 5.5,  # 铝需求，10%负荷
    # "al_excess_rate": [0, 0.05, 0.1, 0.15, 0.2],  # 铝过剩率
    "al_excess_rate": 0.3,  # 铝过剩率
    
    # 切负荷设置
    "enable_load_shedding": True,  # 是否启用切负荷
    "load_shedding_cost": 1e2,  # 普通负荷切负荷成本（欧元/kWh）
    "al_load_shedding_cost": 1e3,  # 铝负荷切负荷成本（欧元/kWh）
    
    # 电解槽成本
    "al_capital_cost": 0 * 6210,  # 电解槽年准化资本成本（$/MW/年）
    
    # 铝存储成本($/MWh), 计算方式：
    # 1/13.4 是每MWh生产的铝(吨)
    # 0.74 是一吨铝需要的空间（立方米）
    # 0.8 是存储价格（0.8元/平方米/天）
    # 1e-6 是转换为百万
    # 7.55 是欧元兑人民币汇率
    # 24 是换算成小时
    "al_marginal_cost_storage": 3.91479E-05 * 7.3 / 7.55,
    "al_committable": True,
    # storage limit
    "al_storage_limit": 24 * 30, # hours
    
    # 电解槽启动成本, 计算方式：
    # $65M ~ 450,000 mt per year
    "al_start_up_cost": 16000 * 7.3 / 7.55, # euro/MW
    "al_p_min_pu": [0.9],  # 改回列表以支持pmin脚本
    # 启动时间
    # "al_start_up_time": 24, # 小时
    
    # 迭代优化参数
    "max_iterations": 10,
    "convergence_tolerance": 0.01,  # 目标函数相对变化收敛阈值（1%）
    "verbose": True,
    
    # 默认成本
    "default_costs": {
        "FOM": 0,
        "VOM": 0,
        "efficiency": 1,
        "fuel": 0,
        "investment": 0,
        "lifetime": 25,
        "CO2 intensity": 0,
        "discount rate": 0.07,
    }
} 