# 情景维度配置指南

本文档介绍如何在PyPSA-China项目中配置和使用三个情景维度：电解铝厂运行灵活性、原铝需求和电网交互市场机会。

## 概述

项目现在支持三个独立的情景维度，每个维度都有三个级别（low、mid、high），总共可以产生27种不同的情景组合。

### 三个情景维度

1. **电解铝厂运行灵活性 (Smelter operational flexibility)**
   - 控制电解铝厂的运行参数
   - 包括过载率、最小功率、重启允许性和重启成本

2. **原铝需求 (Primary aluminum demand)**
   - 控制铝需求的结构
   - 包括国内需求比例、出口率、回收率和产品寿命

3. **电网交互市场机会 (Grid-interaction market opportunity)**
   - 控制电网交互相关的成本参数
   - 包括VRE、电池、H2存储成本降低和其他灵活需求

## 配置文件设置

### config.yaml中的配置结构

```yaml
aluminum:
  # 情景维度设置
  scenario_dimensions:
    # 1. 电解铝厂运行灵活性
    smelter_flexibility:
      low:
        overcapacity_rate: 0.9  # p.u.
        p_min: 0.9  # p.u.
        allow_restart: false
        restart_costs: 110000  # $/MW
      mid:
        overcapacity_rate: 0.9  # p.u.
        p_min: 0.9  # p.u.
        allow_restart: true
        restart_costs: 16000  # $/MW
      high:
        overcapacity_rate: 0.7  # p.u.
        p_min: 0.7  # p.u.
        allow_restart: true
        restart_costs: 3200  # $/MW
    
    # 2. 原铝需求
    primary_demand:
      low:
        domestic_demand_ratio: 0.8  # 80%
        export_rate: 0.2  # 20%
        recycling_rate: 0.2  # 20%
        product_lifetime: 20  # years
      mid:
        domestic_demand_ratio: 0.7  # 70%
        export_rate: 0.3  # 30%
        recycling_rate: 0.16  # 16%
        product_lifetime: 16  # years
      high:
        domestic_demand_ratio: 0.6  # 60%
        export_rate: 0.4  # 40%
        recycling_rate: 0.12  # 12%
        product_lifetime: 12  # years
    
    # 3. 电网交互市场机会
    grid_interaction:
      low:
        vre_cost_reduction: 0.0  # 成本降低比例
        battery_cost_reduction: 0.0  # 成本降低比例
        h2_storage_cost_reduction: 0.0  # 成本降低比例
        other_flexible_demand: 0.0  # 其他灵活需求比例
      mid:
        vre_cost_reduction: 0.1  # 10%成本降低
        battery_cost_reduction: 0.1  # 10%成本降低
        h2_storage_cost_reduction: 0.1  # 10%成本降低
        other_flexible_demand: 0.05  # 5%其他灵活需求
      high:
        vre_cost_reduction: 0.2  # 20%成本降低
        battery_cost_reduction: 0.2  # 20%成本降低
        h2_storage_cost_reduction: 0.2  # 20%成本降低
        other_flexible_demand: 0.1  # 10%其他灵活需求

  # 当前选择的情景组合 (默认设置为mid-mid-mid)
  current_scenario:
    smelter_flexibility: "mid"
    primary_demand: "mid"
    grid_interaction: "mid"
```

## 使用方法

### 1. 基本使用

```python
from scripts.scenario_utils import load_config, get_scenario_params

# 加载配置
config = load_config()

# 获取默认情景参数
default_params = get_scenario_params(config)

# 获取特定情景参数
high_flex_params = get_scenario_params(
    config,
    smelter_flexibility="high",
    primary_demand="mid",
    grid_interaction="low"
)
```

### 2. 获取单个维度的参数

```python
from scripts.scenario_utils import (
    get_smelter_params,
    get_demand_params,
    get_grid_interaction_params
)

# 获取电解铝厂参数
smelter_params = get_smelter_params(config, "high")

# 获取需求参数
demand_params = get_demand_params(config, "low")

# 获取电网交互参数
grid_params = get_grid_interaction_params(config, "mid")
```

### 3. 生成所有情景组合

```python
from scripts.scenario_utils import generate_scenario_combinations

# 生成所有27种情景组合
combinations = generate_scenario_combinations()

for combo in combinations:
    print(f"情景: {combo['name']}")
    print(f"  电解铝厂灵活性: {combo['smelter_flexibility']}")
    print(f"  原铝需求: {combo['primary_demand']}")
    print(f"  电网交互: {combo['grid_interaction']}")
```

### 4. 在PyPSA网络中使用

```python
# 获取情景参数
params = get_scenario_params(config, "high", "mid", "low")

# 在PyPSA网络中使用电解铝厂参数
smelter = params['smelter_flexibility']
network.generators.loc[aluminum_generators, 'p_min_pu'] = smelter['p_min']
network.generators.loc[aluminum_generators, 'start_up_cost'] = smelter['restart_costs']

# 使用需求参数
demand = params['primary_demand']
aluminum_demand = base_demand * demand['domestic_demand_ratio']

# 使用电网交互参数
grid = params['grid_interaction']
vre_cost_adjusted = base_vre_cost * (1 - grid['vre_cost_reduction'])
```

## 情景组合示例

### 代表性情景组合

1. **保守情景 (low-low-low)**
   - 低电解铝厂灵活性
   - 低原铝需求
   - 低电网交互机会

2. **基准情景 (mid-mid-mid)**
   - 中等电解铝厂灵活性
   - 中等原铝需求
   - 中等电网交互机会

3. **激进情景 (high-high-high)**
   - 高电解铝厂灵活性
   - 高原铝需求
   - 高电网交互机会

4. **高灵活性-低需求-高电网交互 (high-low-high)**
   - 适合研究电解铝厂作为电网灵活性资源的情景

5. **低灵活性-高需求-低电网交互 (low-high-low)**
   - 适合研究传统电解铝厂运行模式的情景

## 参数说明

### 电解铝厂运行灵活性参数

- **overcapacity_rate**: 过载率，表示电解铝厂可以超过额定功率运行的比例
- **p_min**: 最小功率比例，电解铝厂运行时的最小功率要求
- **allow_restart**: 是否允许电解铝厂重启
- **restart_costs**: 重启成本，单位为$/MW

### 原铝需求参数

- **domestic_demand_ratio**: 国内需求比例，相对于总需求的比例
- **export_rate**: 出口率，出口需求占总需求的比例
- **recycling_rate**: 回收率，回收铝占总需求的比例
- **product_lifetime**: 产品寿命，铝产品的平均使用寿命（年）

### 电网交互参数

- **vre_cost_reduction**: VRE（可变可再生能源）成本降低比例
- **battery_cost_reduction**: 电池成本降低比例
- **h2_storage_cost_reduction**: H2存储成本降低比例
- **other_flexible_demand**: 其他灵活需求占总需求的比例

## 运行示例

运行示例脚本来查看所有功能：

```bash
python examples/scenario_example.py
```

运行工具函数测试：

```bash
python scripts/scenario_utils.py
```

## 注意事项

1. **默认设置**: 当前默认情景设置为 `mid-mid-mid`
2. **参数验证**: 所有参数都在合理范围内，但建议在使用前验证
3. **向后兼容**: 新的情景维度配置不会影响现有的功能
4. **扩展性**: 可以轻松添加新的情景级别或维度

## 故障排除

如果遇到问题，请检查：

1. 配置文件格式是否正确
2. 情景名称是否拼写正确（low、mid、high）
3. 参数值是否在合理范围内
4. 依赖包是否正确安装（yaml等） 