# config.py

# 配置参数
CONFIG = {
    "years": [2025],
    "resolution": 4,  # 4小时
    "al_demand": 5.5,  # 铝需求，10%负荷
    "al_excess_rate": [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],  # 铝过剩率
    
    # 电解槽成本
    "al_capital_cost": 0,  # 电解槽年准化资本成本（百万欧元/MW/年）
    
    # 铝存储成本, 计算方式：
    # 1/13.4 是每MWh生产的铝
    # 单位转换为百万欧元/吨
    # 0.74 是一吨铝需要的空间（立方米）
    # 0.8 是存储价格（0.8元/平方米/天）
    # 1e-6 是转换为百万
    # 7.55 是欧元兑人民币汇率
    # 24 是换算成小时
    "al_marginal_cost_storage": 0,
    
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