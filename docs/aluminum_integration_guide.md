# 电解铝需求与装机容量融合到模型说明文档

## 1. 概述

本文档详细说明电解铝（Aluminum）的需求数据和装机容量如何融合到PyPSA-China模型中。电解铝作为高耗能工业负荷，在模型中通过以下方式集成：

- **需求数据**：从情景数据文件中读取各年份的原铝需求
- **装机容量**：从现有电解铝厂数据中读取各省份的装机容量
- **模型组件**：通过Link、Store、Load等组件在PyPSA网络中表示
- **优化算法**：使用迭代优化算法处理电解铝的启停约束和灵活性

## 2. 数据来源

### 2.1 需求数据

**文件路径**：`data/aluminum_demand/aluminum_demand_all_scenarios.json`

**数据结构**：
```json
{
  "primary_aluminum_demand": {
    "low": {
      "2030": 1909.875997700051,  // 单位：10kt (万吨)
      "2050": 628.2308472342488,
      ...
    },
    "mid": {...},
    "high": {...}
  }
}
```

**数据说明**：
- 包含三个需求情景：`low`、`mid`、`high`
- 覆盖年份：2030-2060
- 单位：10kt（万吨）
- 数据来源：EP模型计算结果

### 2.2 装机容量数据

**文件路径**：`data/p_nom/al_smelter_p_max.csv`

**数据结构**：
```csv
Province,p_nom
Anhui,161.492652
Gansu,300.120354
Shandong,505.749267
...
```

**数据说明**：
- 单位：10kt/年（万吨/年）
- 包含各省份的电解铝年产量数据
- 用于计算各省份的装机容量和 production ratio

## 3. 需求数据处理流程

### 3.1 需求数据读取

需求数据通过 `scripts/scenario_utils.py` 中的函数读取：

```python
def get_aluminum_demand_for_year(config, year, primary_demand_scenario=None, 
                                  aluminum_demand_json_path=...):
    """
    获取指定年份和情景的原铝需求
    
    返回：原铝需求（吨）
    """
    # 1. 从config获取需求情景（如未指定，使用current_scenario）
    # 2. 从JSON文件读取对应年份和情景的需求数据（10kt）
    # 3. 转换为吨：primary_demand_tons = primary_demand_10kt * 10000
    return primary_demand_tons
```

**关键代码位置**：`scripts/scenario_utils.py:195-224`

### 3.2 需求转换为负荷

需求数据转换为网络负荷通过 `get_aluminum_load_for_network()` 函数：

```python
def get_aluminum_load_for_network(config, year, network_snapshots, nodes, 
                                  production_ratio, ...):
    """
    将原铝需求转换为网络负荷
    
    流程：
    1. 获取原铝需求（吨）
    2. 计算全国平均铝负荷：national_al_load = primary_demand_tons / 8760 (MW)
    3. 根据各省份production_ratio分配负荷
    4. 创建时间序列负荷数据
    """
    # 计算全国铝负荷 (MW)
    hours_per_year = 8760
    national_al_load = primary_demand_tons / hours_per_year
    
    # 按省份比例分配
    al_load_values = np.tile(
        national_al_load * production_ratio.values,
        (len(network_snapshots), 1)
    )
    
    return {
        'aluminum_load': aluminum_load,  # DataFrame，时间序列×省份
        'national_al_load': national_al_load,
        'primary_demand_tons': primary_demand_tons
    }
```

**关键代码位置**：`scripts/scenario_utils.py:227-274`

**转换公式**：
- 全国平均负荷 (MW) = 原铝需求 (吨) / 8760 小时
- 省份负荷 = 全国平均负荷 × 该省份的 production_ratio

## 4. 装机容量处理流程

### 4.1 容量数据读取

**代码位置**：`scripts/prepare_base_network.py:201-224`

```python
# 读取电解铝厂容量数据
al_smelter_annual_production = pd.read_csv(snakemake.input.al_smelter_p_max)
al_smelter_annual_production = al_smelter_annual_production.set_index('Province')['p_nom']

# 过滤：只保留年产量 > 0.01 10kt/year 的省份
al_smelter_annual_production = al_smelter_annual_production[
    al_smelter_annual_production > 0.01
]

# 计算生产比例
production_ratio = al_smelter_annual_production / al_smelter_annual_production.sum()
```

### 4.2 容量转换公式

**代码位置**：`scripts/prepare_base_network.py:210-219`

```python
# 转换为功率容量 (MW)
# 转换公式：
# 1. 年产量 (10kt/year) → 年产量 (吨/year) = 年产量 × 10000
# 2. 年产量 (吨/year) → 年用电量 (MWh/year) = 年产量 × 13.3
#    注：1吨铝需要约13.3 MWh电力
# 3. 年用电量 (MWh/year) → 功率容量 (MW) = 年用电量 / 8760

base_capacity = al_smelter_annual_production * 10000 * 13.3 / 8760

# 应用容量比例（可配置）
capacity_ratio = config['aluminum']['capacity_ratio']  # 默认1.0
al_smelter_p_nom = base_capacity * capacity_ratio
```

**转换公式总结**：
```
p_nom (MW) = 年产量(10kt/year) × 10000 × 13.3 / 8760 × capacity_ratio
```

**参数说明**：
- `13.3`：生产1吨铝所需的电力（MWh/吨）
- `8760`：年小时数
- `capacity_ratio`：容量比例配置（可在config.yaml中设置，默认1.0）

## 5. 模型组件集成

### 5.1 网络组件添加

**代码位置**：`scripts/prepare_base_network.py:243-283`

#### 5.1.1 电解铝冶炼设备（Link）

```python
network.madd("Link",
            production_ratio.index,  # 省份名称列表
            suffix=" aluminum smelter",
            bus0=production_ratio.index,  # 电力节点（输入）
            bus1=production_ratio.index + " aluminum",  # 铝节点（输出）
            carrier="aluminum",
            p_nom=al_smelter_p_nom,  # 装机容量 (MW)
            p_nom_extendable=False,  # 不可扩展
            efficiency=1.0/13.3,  # 效率：1 MW电力 → 1/13.3 吨铝/小时
            capital_cost=operational_params['capital_cost'],
            stand_by_cost=operational_params['stand_by_cost'],
            marginal_cost=operational_params['marginal_cost'],
            start_up_cost=0.5*operational_params['start_up_cost'],
            shut_down_cost=0.5*operational_params['start_up_cost'],
            committable=config['aluminum_commitment'],  # 是否启用启停约束
            p_min_pu=operational_params['p_min_pu'] if config['aluminum_commitment'] else 0,
)
```

**关键参数**：
- `efficiency=1.0/13.3`：表示1 MW电力输入可产生 1/13.3 吨铝/小时
- `p_nom`：装机容量，从数据文件计算得到
- `committable`：是否启用启停约束（MILP问题）

#### 5.1.2 铝存储（Store）

```python
network.madd("Store",
            production_ratio.index,
            suffix=" aluminum storage",
            bus=production_ratio.index + " aluminum",
            carrier="aluminum",
            e_nom_extendable=True,  # 存储容量可扩展
            e_cyclic=True)  # 循环存储（年末=年初）
```

**功能**：允许铝在时间上转移，提供运行灵活性

#### 5.1.3 铝负荷（Load）

```python
network.madd("Load",
            production_ratio.index,
            suffix=" aluminum",
            bus=production_ratio.index + " aluminum",
            p_set=aluminum_load[production_ratio.index])  # 时间序列负荷
```

**功能**：表示各省份的铝需求负荷

#### 5.1.4 电力负荷调整

```python
# 从电力负荷中减去铝负荷（避免重复计算）
load_minus_al = load.copy()
load_minus_al[production_ratio.index] = (
    load[production_ratio.index] - 
    aluminum_load[production_ratio.index] * 10000 * 13.3 / 8760
)
network.madd("Load", nodes, bus=nodes, p_set=load_minus_al)
```

**原因**：铝负荷已单独建模，需要从总电力负荷中减去

#### 5.1.5 全国铝枢纽（China Aluminum Hub）

```python
# 添加全国铝枢纽节点，支持省份间铝转移
network.add("Bus", "China aluminum hub", carrier="aluminum transfer")

# 添加双向转移链路
for province in production_ratio.index:
    # 省份 → 全国枢纽
    network.add("Link", f"{province} to China aluminum hub", ...)
    # 全国枢纽 → 省份
    network.add("Link", f"China aluminum hub to {province}", ...)
```

**功能**：允许省份间铝的转移，满足全国总需求约束

### 5.2 运行参数设置

**代码位置**：`scripts/scenario_utils.py:153-192`

运行参数根据配置的情景维度设置：

```python
def get_aluminum_smelter_operational_params(config, smelter_flexibility=None, 
                                            al_smelter_p_nom=None):
    """
    获取电解铝厂运行参数
    
    参数来源：config['aluminum']['scenario_dimensions']['smelter_flexibility']
    """
    smelter_params = get_smelter_params(config, smelter_flexibility)
    
    return {
        'p_min_pu': smelter_params['p_min_pu'],  # 最小出力比例
        'capital_cost': 163432.8,  # 固定资本成本
        'stand_by_cost': smelter_params['stand_by_cost'] * al_smelter_p_nom,  # $/h
        'marginal_cost': 1,  # 边际成本
        'start_up_cost': smelter_params['restart_cost'] * al_smelter_p_nom,  # 启动成本
    }
```

**情景参数**（来自config.yaml）：
- `low`：p_min_pu=0.9, restart_cost=110000 $/MW
- `mid`：p_min_pu=0.7, restart_cost=16000 $/MW
- `high`：p_min_pu=0.5, restart_cost=3200 $/MW
- `non_constrained`：p_min_pu=0.0, restart_cost=0

## 6. 迭代优化算法

### 6.1 算法概述

由于电解铝启停约束导致MILP问题，使用迭代优化算法：

**代码位置**：`scripts/solve_network_myopic.py:583-1163`

**算法流程**：

1. **初始化**：设置电解铝用能模式为空
2. **迭代求解**：
   - **步骤1**：使用连续化电解铝模型（无启停约束）求解，获得节点电价和目标函数值
   - **步骤2**：基于节点电价，运行电解铝最优运行问题（MILP），得到新的电解铝用能模式
   - **步骤3**：使用`p_set`固定电解铝用能，重新求解网络
   - **步骤4**：检查目标函数变化，判断是否收敛
3. **收敛判断**：目标函数相对变化 < 收敛阈值（默认1%）
4. **输出结果**：返回最终网络结果

### 6.2 关键实现细节

#### 6.2.1 网络重新加载

每次迭代重新加载网络，保持状态清洁：

```python
# 重新加载网络
if original_network_path:
    n_current = pypsa.Network(original_network_path, override_component_attrs=overrides)
else:
    n_current = copy.deepcopy(original_network)

# 重新应用网络准备
n_current = prepare_network(n_current, ...)
```

#### 6.2.2 电解铝用能固定

使用`p_set`固定电解铝用能，而不是添加约束：

```python
# 固定电解铝冶炼设备的出力
for smelter in aluminum_usage.columns:
    if smelter in n_current.links.index:
        fixed_aluminum_power = aluminum_usage[smelter].values
        n_current.links_t.p_set[smelter] = fixed_aluminum_power
        
        # 同时固定对应的铝负荷
        for load in aluminum_loads:
            n_current.loads_t.p_set[load] = fixed_aluminum_power
```

#### 6.2.3 节点电价提取

从连续化模型求解结果中提取节点电价：

```python
# 提取节点电价
if hasattr(n_current, 'buses_t') and hasattr(n_current.buses_t, 'marginal_price'):
    electricity_buses = n_current.buses[n_current.buses.carrier != "aluminum"].index
    current_nodal_prices = n_current.buses_t.marginal_price[electricity_buses]
```

#### 6.2.4 电解铝优化问题

基于节点电价，求解电解铝最优运行问题：

```python
def solve_aluminum_optimization(n, config, solving, opts="", 
                                nodal_prices=None, target_province=None, ...):
    """
    电解铝最优运行问题求解
    
    目标：最小化总成本（电力成本 + 启停成本）
    约束：
    1. 满足铝需求
    2. 容量约束
    3. 启停约束（如果启用）
    """
    # 创建简化的电解铝优化网络
    # 添加虚拟发电机（使用节点电价作为边际成本）
    # 求解MILP问题
    # 返回电解铝用能模式
```

## 7. 配置参数说明

### 7.1 基本配置

**文件位置**：`config.yaml`

```yaml
aluminum:
  # 电网交互开关（按年份）
  grid_interaction:
    "2020": false
    "2030": true
    "2050": true
    ...

  # 电解铝厂容量比例
  capacity_ratio: 1.0  # 1.0 = 100%, 0.9 = 90%, ...

  # 当前选择的情景组合
  current_scenario:
    smelter_flexibility: "low"    # 电解铝厂灵活性
    primary_demand: "high"        # 原铝需求情景
    market_opportunity: "high"    # 市场机会情景
```

### 7.2 情景维度参数

```yaml
aluminum:
  scenario_dimensions:
    # 电解铝厂运行灵活性
    smelter_flexibility:
      low:
        p_min_pu: 0.9          # 最小出力比例
        restart_cost: 110000   # 重启成本 ($/MW)
        stand_by_cost: 1.5     # 待机成本 ($/MW/h)
      mid:
        p_min_pu: 0.7
        restart_cost: 16000
        stand_by_cost: 1.5
      high:
        p_min_pu: 0.5
        restart_cost: 3200
        stand_by_cost: 1.5
```

### 7.3 迭代优化参数

```yaml
# 在solve_opts中设置
aluminum_max_iterations: 10              # 最大迭代次数
aluminum_convergence_tolerance: 0.01     # 收敛阈值（1%）
aluminum_commitment: false                # 是否启用启停约束
```

## 8. 数据流图

```
┌─────────────────────────────────────────────────────────────┐
│                   数据输入                                    │
├─────────────────────────────────────────────────────────────┤
│ 1. 需求数据: aluminum_demand_all_scenarios.json            │
│    └─> 原铝需求 (10kt) × 情景 × 年份                        │
│                                                              │
│ 2. 容量数据: al_smelter_p_max.csv                           │
│    └─> 各省份年产量 (10kt/year)                             │
└─────────────────────────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────┐
│                   数据处理                                    │
├─────────────────────────────────────────────────────────────┤
│ 1. 需求处理:                                                │
│    primary_demand (10kt)                                    │
│    └─> primary_demand_tons (吨)                             │
│        └─> national_al_load (MW)                            │
│            └─> 按production_ratio分配到各省份                │
│                └─> aluminum_load (时间序列 × 省份)          │
│                                                              │
│ 2. 容量处理:                                                │
│    al_smelter_annual_production (10kt/year)                 │
│    └─> base_capacity (MW) = 年产量 × 10000 × 13.3 / 8760   │
│        └─> al_smelter_p_nom (MW) = base_capacity × ratio    │
│                                                              │
│ 3. 比例计算:                                                │
│    production_ratio = 各省份产量 / 全国总产量                 │
└─────────────────────────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────┐
│                   网络组件添加                                │
├─────────────────────────────────────────────────────────────┤
│ 1. Link: 电解铝冶炼设备                                      │
│    bus0: 电力节点 → bus1: 铝节点                            │
│    p_nom: 装机容量, efficiency: 1/13.3                      │
│                                                              │
│ 2. Store: 铝存储                                             │
│    bus: 铝节点, e_nom_extendable: True                      │
│                                                              │
│ 3. Load: 铝负荷                                             │
│    bus: 铝节点, p_set: 时间序列负荷                          │
│                                                              │
│ 4. Load: 电力负荷（已减去铝负荷）                            │
│                                                              │
│ 5. Bus + Link: 全国铝枢纽（省份间转移）                      │
└─────────────────────────────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────┐
│                   迭代优化求解                                │
├─────────────────────────────────────────────────────────────┤
│ 迭代1:                                                       │
│   1. 连续化模型求解 → 节点电价                               │
│   2. 基于节点电价求解MILP → 电解铝用能模式                   │
│   3. 固定电解铝用能，重新求解                                 │
│   4. 检查收敛性                                              │
│                                                              │
│ 迭代2-N:                                                     │
│   重复上述步骤，直到收敛                                      │
└─────────────────────────────────────────────────────────────┘
```

## 9. 关键公式总结

### 9.1 容量转换

```
p_nom (MW) = 年产量(10kt/year) × 10000 × 13.3 / 8760 × capacity_ratio
```

其中：
- `10000`：10kt → 吨的转换
- `13.3`：生产1吨铝所需的电力（MWh/吨）
- `8760`：年小时数
- `capacity_ratio`：容量比例（可配置）

### 9.2 需求转换

```
全国平均负荷 (MW) = 原铝需求 (吨) / 8760
省份负荷 (MW) = 全国平均负荷 × 该省份的 production_ratio
```

### 9.3 效率参数

```
Link效率 = 1.0 / 13.3
含义：1 MW电力输入 → 1/13.3 吨铝/小时输出
```

## 10. 相关文件清单

### 10.1 数据文件
- `data/aluminum_demand/aluminum_demand_all_scenarios.json`：需求数据
- `data/p_nom/al_smelter_p_max.csv`：装机容量数据

### 10.2 代码文件
- `scripts/prepare_base_network.py`：网络准备，添加电解铝组件
- `scripts/scenario_utils.py`：情景工具函数，需求/容量处理
- `scripts/solve_network_myopic.py`：迭代优化算法实现
- `config.yaml`：配置文件

### 10.3 文档文件
- `docs/README_aluminum_iterative.md`：迭代算法说明
- `docs/aluminum_integration_guide.md`：本文档

## 11. 使用示例

### 11.1 启用电解铝功能

在`config.yaml`中设置：

```yaml
add_aluminum: true
aluminum:
  grid_interaction:
    "2030": true
  capacity_ratio: 1.0
  current_scenario:
    smelter_flexibility: "mid"
    primary_demand: "mid"
```

### 11.2 调整容量比例

```yaml
aluminum:
  capacity_ratio: 0.8  # 使用80%的容量
```

### 11.3 设置迭代参数

```yaml
solve_opts:
  aluminum_max_iterations: 10
  aluminum_convergence_tolerance: 0.01
  aluminum_commitment: true  # 启用启停约束
```

## 12. 注意事项

1. **单位转换**：注意数据单位（10kt vs 吨，MW vs MWh）
2. **时间序列**：铝负荷是恒定的时间序列（无波动）
3. **省份过滤**：只对年产量 > 0.01 10kt/year 的省份添加组件
4. **电力负荷调整**：需要从总电力负荷中减去铝负荷，避免重复计算
5. **迭代收敛**：如果迭代不收敛，检查收敛阈值和最大迭代次数设置
6. **MILP求解**：启用启停约束时，需要支持MILP的求解器（如Gurobi）

## 13. 常见问题

### Q1: 如何修改铝需求情景？
A: 在`config.yaml`中修改`aluminum.current_scenario.primary_demand`，可选值：`low`、`mid`、`high`

### Q2: 如何调整电解铝厂容量？
A: 修改`config.yaml`中的`aluminum.capacity_ratio`参数（0.0-1.0）

### Q3: 迭代不收敛怎么办？
A: 检查：
- 收敛阈值是否过小（建议0.01-0.05）
- 最大迭代次数是否足够（建议10-20）
- 求解器设置是否合理

### Q4: 如何禁用电解铝功能？
A: 设置`add_aluminum: false`或`aluminum.grid_interaction[年份]: false`

---

**文档版本**：v1.0  
**最后更新**：2025年  
**维护者**：Ruike Lyu

