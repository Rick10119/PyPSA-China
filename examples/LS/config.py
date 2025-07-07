# config.py
"""
 使用Pypsa数据时，因为单位为欧元，因此按照2020年人民币对欧元平均汇率(约7.85)进行转换。

 新能源机组参数采用中国2020年数据，投资成本使用中国2020数据，边际成本采用PyPSA数据（边际成本相对很小）。

 火电机组参数采用中国2020年数据，投资成本暂定为3000万元/MW，边际成本简化为分段成本的均值。

 水电机组参数采用中国2020年数据，投资成本暂时使用Pypsa数据（中国数据待搜索），边际成本简化为分段成本的均值。

 储能单元参数采用中国2020年数据，投资成本和边际成本需要搜索。

"""
# 配置参数
CONFIG = {
    ## 经济参数

    # "snapshots": 24 * 7,  # 模拟时长：7天，每小时一个数据点
    "country_code": "CN", # 中国
    "technologies": ["solar", "onwind", "thermal", "hydro"], # 候选发电技术
    "coal_types": ["coal_320", "coal_350", "coal_660", "coal_1000"], # 火电机组类型
    "hydro_types": ["hydro_300", "hydro_400"], # 水电机组类型
    "storage_types": ["store_100","store_200","store_300","store_400","store_500","store_600_1","store_600_2"], # 储能类型
    "costs": {
        "solar": {
            "capital_cost": 450,   # 投资成本 [万元/MW]  中国2020数据
            "marginal_cost": 0.001,  # 边际成本 [万元/MWh]  pypsa数据
            "FOM": 1.58  # 固定运营成本 [%/year]  pypsa数据
        },
        "onwind": {
            "capital_cost": 800, # 投资成本 [万元/MW]  中国2020数据
            "marginal_cost": 0.002, # 边际成本 [万元/MWh]  pypsa数据
            "FOM": 1.25  # 固定运营成本 [%/year]  pypsa数据
        },
        # 火电成本投资不是特别可靠，暂定均为3000万元/MW
        "thermal": {
            "startup_cost": 380,  # 启动成本 [万元/次]
            "shutdown_cost": 240, # 关机成本 [万元/次]
            "coal_320": {
                "capital_cost": 3000, # 投资成本 [万元/MW]
                "marginal_cost": 0.447, # 边际成本 [万元/MWh]  
            },
            "coal_350": {
                "capital_cost": 3000, # 投资成本 [万元/MW]
                "marginal_cost": 0.4145, # 边际成本 [万元/MWh]  
            },
            "coal_660": {
                "capital_cost": 3000, # 投资成本 [万元/MW]
                "marginal_cost": 0.429, # 边际成本 [万元/MWh]  
            },
            "coal_1000": {
                "capital_cost": 3000, # 投资成本 [万元/MW]
                "marginal_cost": 0.401, # 边际成本 [万元/MWh]  
            }
        },
        # 水电投资成本需要查找，边际成本和投资成本都比较有问题
        "hydro": {
            "startup_cost": 120, # 启动成本 [万元/次]
            "shutdown_cost": 20, # 关机成本 [万元/次]
            "FOM": 1.0, # 固定运营成本 [%/year]  pypsa数据
            "hydro_300": {
                "capital_cost": 1733, # 投资成本 [万元/MW] Pypsa数据
                "marginal_cost": 0.341, # 边际成本 [万元/MWh]   
            },
            "hydro_400": {
                "capital_cost": 1733, # 投资成本 [万元/MW] Pypsa数据
                "marginal_cost": 0.341, # 边际成本 [万元/MWh]  
            }
        }
    },

    ## 机组参数

    # 新能源机组相关参数
    # 安徽省新能源机组参数配置，粗略认为整个安徽省只有一台光伏发电机组和一台风力发电机组
    # 出力为安徽省2020年光伏和风力发电的平均出力
    "p_nom_solar": 13700,  # 光伏发电机组的容量 [MW]
    "p_nom_onwind": 4120,   # 风力发电机组的容量 [MW]

    # 热电机组相关参数
    # 安徽省火电机组参数配置
    "thermal_units": {
        "min_up_time": 2,      # 最小运行时间 [h]
        "min_down_time": 1,    # 最小停机时间 [h]
        "coal_320": {
            "capacity": 32000,     # 机组容量 [MW]
            "p_min": 12500/32000,  # 最小出力 p.u
            "ramp_rate": 19500,    # 爬坡率 [MW/h]
            "count": 18            # 已建数量
        },
        "coal_350": {
            "capacity": 35000,     # 机组容量 [MW]
            "p_min": 16500/35000,  # 最小出力 p.u
            "ramp_rate": 18500,    # 爬坡率 [MW/h]
            "count": 12            # 已建数量
        },
        "coal_660": {
            "capacity": 66000,     # 机组容量 [MW]
            "p_min": 18000/66000,  # 最小出力 p.u
            "ramp_rate": 48000,    # 爬坡率 [MW/h]
            "count": 60            # 已建数量
        },
        "coal_1000": {
            "capacity": 100000,    # 机组容量 [MW]
            "p_min": 45000/100000, # 最小出力 p.u
            "ramp_rate": 55000,    # 爬坡率 [MW/h]
            "count": 12            # 已建数量
        }
    },

    ## 水电参数
    # 水电机组相关参数
    "hydro_units": {
        "hydro_300": {
            "capacity": 30000,      # 水电机组容量 [MW]
            "p_min": 0.5,          # 最小出力 p.u
            "ramp_rate": 15000,      # 爬坡率 [MW/h]
            "count": 3             # 已建数量
        },
        "hydro_400": {
            "capacity": 40000,      # 水电机组容量 [MW]
            "p_min": 0.5,          # 最小出力 p.u
            "ramp_rate": 20000,      # 爬坡率 [MW/h]
            "count": 7             # 已建数量
        }
    },

    ## 储能参数 需要确认单位
    # 储能单元相关参数 - 使用矩阵格式定义
    "storage_params": {
        # 全局参数
        "efficiency_dispatch": 0.95,  # 放电效率
        "efficiency_store": 0.95,     # 充电效率
        "soc_max": 0.9,               # SOC上限
        "soc_min": 0.1,               # SOC下限
        
        # 矩阵定义: 列为储能类型，行为参数
        "matrix": {
            # 储能类型           store_100  store_200  store_300  store_400  store_500  store_600_1  store_600_2
            "capacity":         [100,      200,       300,       400,       500,       600,         600],      # 容量 [MW]
            "p_charge_max":     [50,       100,       150,       200,       200,       300,         200],      # 最大充电功率 [MW]
            "p_discharge_max":  [50,       100,       150,       200,       200,       300,         200],      # 最大放电功率 [MW]
            "degradation_cost": [66.79,    66.79,     66.79,     66.79,     25,        25,         66.79],    # 退化成本 [万元/循环]
            "count":            [2,        7,         1,         4,         1,         1,           1],        # 已建数量
            "eff_override":     [False,    False,     False,     False,     True,     True,        False]     # 是否覆盖全局效率
        },
        
        # 特殊效率覆盖
        "special_efficiency": {
                "efficiency_dispatch": 0.85, 
                "efficiency_store": 0.85
        }
    },

    ## LS配置参数（AI拟定）
    # Define considered products
    "materials": [
        # Upstream raw materials
        "bauxite", "iron_ore", "coke", "limestone", "oil", "clay",
        "alumina", "silica_sand",
        
        # Primary products
        "aluminum", "steel", "basic_chemicals", "cement", "glass",
        
        # Downstream raw materials 
        "copper",
        
        # Downstream products
        "vehicles", "machine_tools", "heavy_machinery", "buildings",
        "aircraft_parts", "consumer_electronics"
    ],


    # Define initial inventory for each product at start of optimization
    "initial_inventories": { 
        # Upstream raw materials
        "bauxite": 10e8,
        "iron_ore": 10e8,
        "coke": 10e8,
        "limestone": 10e8,
        "oil": 10e8,
        "clay": 10e8,
        "alumina": 10e8,
        "silica_sand": 10e8,
        # Primary products
        "aluminum": 1,
        "steel": 1,
        "basic_chemicals": 1,
        "cement": 1,
        "glass": 1,
        # Downstream raw materials
        "copper": 10e8,
        # Downstream products 
        "vehicles": 0,
        "machine_tools": 0,
        "heavy_machinery": 0,
        "buildings": 0,
        "aircraft_parts": 0,
        "consumer_electronics": 0
    },


    # 定义每种产品的需求量，先认为原料的需求量为0，
    # 2级产品的需求量为一个较大的数值，表示下游产业的需求
    "demand_matrix": [
        # 原料
        0,      # 铝土矿
        0,      # 铁矿石
        0,      # 焦炭
        0,      # 石灰石
        0,      # 石油
        0,      # 黏土
        0,      # 氧化铝 
        0,      # 石英砂
        # 1级产品 （假设需求较小）
        300,   # 铝材
        1200,  # 钢铁
        75,    # 基础化学品
        500,   # 水泥
        250,   # 玻璃
        # ##
        0,      # 铜材
        # 2级产品
        800,   # 整车
        700,   # 数控机床
        700,   # 重型机械
        1200,  # 建筑物
        1000,  # 飞机部件
        3000   # 消费电子产品
    ],



    # 定义一个1×n矩阵，代表每个产品的均价，由于原料假定为无需求，因此不设置均价
    "price_matrix": [
        # 单位：万元（AI拟定价格）
        # 原料
        0,      # 铝土矿
        0,      # 铁矿石
        0,      # 焦炭
        0,      # 石灰石
        0,      # 石油
        0,      # 黏土
        0,      # 氧化铝
        0,      # 石英砂
        # 1级产品
        1.9,    # 铝材
        0.4,    # 钢铁
        0.7,    # 基础化学品
        0.04,   # 水泥
        0.2,    # 玻璃
        0,      # 铜材
        # 2级产品
        12,   # 整车
        35,   # 数控机床
        20,   # 重型机械
        15,   # 建筑物 （这个ai一开始定的1m^2 = 1.5万，但我觉得可以高一些）
        80,   # 飞机部件
        100   # 消费电子产品 （这个很高，而且低能耗用的原料也少，可能会成为决定性因素）
    ],

    # 定义各产品生产的单位能耗（MWh/单位产品）
    "energy_consumption": [
        # 上游产业原料
        1,      # 铝土矿
        1,      # 铁矿石 
        1,      # 焦炭
        1,      # 石灰石
        1,      # 石油
        1,      # 黏土
        1,      # 氧化铝
        1,      # 石英砂
        # 上游产业产品
        13.6,   # 铝材
        4.6,    # 钢铁
        5.3,    # 基础化学品
        0.8,    # 水泥
        1.8,    # 玻璃
        # 下游产业原料
        1,      # 铜材
        # 下游产业产品
        0.14,   # 整车
        2.5,    # 数控机床
        3.0,    # 重型机械
        4.0,    # 建筑物
        7.5,    # 飞机部件
        0.1     # 消费电子产品
    ],

    # 定义一个20x20的关系矩阵
    "relationship_matrix": [
        # 铝土矿  铁矿石  焦炭    石灰石  石油    黏土    氧化铝  石英砂  铝材    钢铁    基础化学品 水泥   玻璃   铜材   整车   数控机床 重型机械 建筑物 飞机部件 消费电子产品
        [0,      0,      0,      0,      0,      0,      0,      0,      4.5,    0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0],      # 铝土矿
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      1.6,    0,      0,      0,      0,      0,      0,      0,      0,      0,      0],      # 铁矿石
        [0,      0,      0,      0,      0,      0,      0,      0,      0.01,   0.5,    0.01,   0,      0.01,   0,      0,      0,      0,      0,      0,      0],      # 焦炭
        [0,      0,      0,      0,      0,      0,      0,      0,      0.01,   0.25,   0,      1.4,    0.2,    0,      0,      0,      0,      0,      0,      0],      # 石灰石
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0.001,  1.5,    0.02,   0,      0,      0,      0,      0,      0,      0,      0],      # 石油
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0.01,   0.05,   0,      0,      0,      0,      0,      0,      0,      0],      # 黏土
        [0,      0,      0,      0,      0,      0,      0,      0,      2,      0,      0.01,   0.01,   0.05,   0,      0,      0,      0,      0,      0,      0],      # 氧化铝
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0.7,    0,      0,      0,      0,      0,      0,      0],      # 石英砂
        [0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0.18,   0.05,   0.05,   0,      0.6,    0.04],   # 铝材
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0.85,   4.75,   4.75,   0.07,   0.15,   0],      # 钢铁
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0.25,   0.02,   0.02,   0.005,  0.1,    0.05],   # 基础化学品
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0,      0.4,    0,      0],      # 水泥
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0.04,   0,      0,      0.015,  0,      0.07],   # 玻璃
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0.1,    0.1,    0,      0.02,   0.015],  # 铜材
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0,      0],      # 整车
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0,      0],      # 数控机床
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0,      0],      # 重型机械
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0,      0],      # 建筑物
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1,      0],      # 飞机部件
        [0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      1]       # 消费电子产品
    ]

} 

