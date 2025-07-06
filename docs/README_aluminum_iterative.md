# 电解铝迭代优化算法改写说明

## 改写概述

根据 `capacity_expansion_planning_aluminum_iterative.py` 的框架，对 `solve_network_aluminum_iterative.py` 进行了重构，主要改进包括：

## 主要改进

### 1. 迭代逻辑优化
- **收敛判断方式**：从基于电解铝用能变化改为基于目标函数变化的相对值
- **收敛阈值**：默认设置为 1%（0.01），比原来的电解铝用能变化阈值更合理
- **时间统计**：添加了详细的迭代时间统计，包括总时间、平均时间、最快/最慢迭代时间

### 2. 网络重新创建方式 ⭐
- **重新加载网络**：参考 `capacity_expansion_planning_aluminum_iterative.py` 的实现，每次迭代都重新加载网络，而不是复制
- **保持网络状态清洁**：避免网络对象状态的累积和污染
- **更稳定的求解**：每次迭代都从干净的网络状态开始

### 3. 电解铝用能固定方式改进 ⭐
- **使用p_set固定**：参考 `capacity_expansion_planning_aluminum_iterative.py` 的实现，使用 `p_set` 来固定电解铝用能，而不是通过约束
- **同时固定负荷**：不仅固定电解铝冶炼设备的出力，还同时固定对应的电解铝负荷
- **更简洁的实现**：删除了复杂的约束添加逻辑，代码更简洁高效

### 4. 虚拟发电机边际成本改进 ⭐
- **使用节点边际电价**：虚拟发电机的边际成本设置为节点边际电价，而不是零边际成本
- **节点对应关系**：确保虚拟发电机与对应的电力节点正确关联
- **更准确的优化**：基于真实的节点电价进行电解铝优化，提高优化精度

### 5. 电解铝优化网络简化
- **组件识别**：简化了电解铝相关组件的识别和保留逻辑
- **虚拟发电机**：改进了虚拟发电机的添加方式，确保为电力节点正确添加
- **Carrier管理**：确保虚拟carrier存在，避免运行时错误

### 6. 优化方法简化
- **删除safe_optimize**：移除了复杂的 `safe_optimize` 函数，直接使用 PyPSA 的标准优化方法
- **错误处理简化**：简化了错误处理逻辑，减少了代码复杂度
- **状态检查简化**：移除了复杂的状态和条件检查

### 7. 配置参数支持
- **迭代参数**：支持通过配置文件设置最大迭代次数和收敛阈值
- **功能开关**：通过 `add_aluminum` 配置项控制是否启用电解铝迭代优化

## 配置参数

在配置文件中可以设置以下参数：

```yaml
# 电解铝功能开关
add_aluminum: true

# 迭代优化参数
aluminum_max_iterations: 10          # 最大迭代次数
aluminum_convergence_tolerance: 0.01 # 收敛阈值（1%）

# 电解铝启停约束
aluminum_commitment: false           # 是否启用启停约束
```

## 主要函数说明

### `solve_network_iterative()`
- **功能**：电解铝迭代优化算法主函数
- **收敛条件**：目标函数相对变化 < 收敛阈值
- **网络处理**：每次迭代重新加载网络，保持状态清洁
- **输出**：时间统计和最终网络结果

### `solve_aluminum_optimization()`
- **功能**：电解铝最优运行问题求解
- **输入**：原始网络、配置、求解参数
- **输出**：电解铝用能模式

### `extra_functionality()`
- **功能**：添加额外的约束条件
- **包含**：CHP约束、传输约束、改造约束等
- **注意**：不再包含电解铝用能约束，因为现在使用p_set固定

## 使用方式

1. **启用电解铝功能**：在配置文件中设置 `add_aluminum: true`
2. **设置迭代参数**：根据需要调整 `aluminum_max_iterations` 和 `aluminum_convergence_tolerance`
3. **运行求解**：使用标准的 PyPSA-China 工作流程运行

## 算法流程

1. **初始化**：设置电解铝用能模式为空，目标函数值为空
2. **迭代求解**：
   - 步骤1：重新加载网络，使用连续化电解铝模型求解，获得节点电价
   - 步骤2：基于节点电价，运行电解铝最优运行问题
   - 步骤3：检查目标函数变化，判断是否收敛
3. **收敛判断**：目标函数相对变化 < 收敛阈值时停止
4. **输出结果**：返回最终网络结果和时间统计

## 网络重新创建方式

### 旧方式（复制网络）
```python
# 复制网络对象
n_current = copy.deepcopy(n)
n_current.config = config
n_current.opts = opts
```

### 新方式（重新加载网络）⭐
```python
# 重新加载网络，保持状态清洁
if original_network_path:
    # 从文件重新加载网络
    if "overrides" in kwargs:
        overrides = kwargs["overrides"]
        n_current = pypsa.Network(original_network_path, override_component_attrs=overrides)
    else:
        n_current = pypsa.Network(original_network_path)
else:
    # 复制原始网络
    n_current = copy.deepcopy(original_network)

# 重新应用网络准备
n_current = prepare_network(
    n_current,
    kwargs.get("solve_opts", {}),
    kwargs.get("using_single_node", False),
    kwargs.get("single_node_province", "Shandong")
)
```

## 电解铝用能固定方式

### 旧方式（约束方式）
```python
# 通过添加约束来固定电解铝用能
def add_aluminum_usage_constraints(n, fixed_aluminum_usage):
    p = n.model["Link-p"]
    for smelter in aluminum_smelters:
        lhs = p.loc[:, smelter]
        rhs = fixed_aluminum_usage[smelter].values
        n.model.add_constraints(lhs == rhs, name=f"aluminum-fixed-usage-{smelter}")
```

### 新方式（p_set方式）⭐
```python
# 通过设置p_set来固定电解铝用能
for smelter in aluminum_usage.columns:
    if smelter in n_current.links.index:
        fixed_aluminum_power = aluminum_usage[smelter].values
        # 确保p_set存在
        if not hasattr(n_current.links_t, 'p_set'):
            n_current.links_t.p_set = pd.DataFrame(index=n_current.snapshots, columns=n_current.links.index)
        # 修改smelter的link出力
        n_current.links_t.p_set[smelter] = fixed_aluminum_power
        # 同时修改对应的负荷设定
        for load in aluminum_loads:
            if hasattr(n_current.loads_t, 'p_set'):
                n_current.loads_t.p_set[load] = fixed_aluminum_power
```

## 虚拟发电机边际成本设置

### 旧方式（零边际成本）
```python
# 虚拟发电机使用零边际成本
n_al.add("Generator",
         f"virtual_gen_{bus}",
         bus=bus,
         carrier="virtual",
         p_nom=1e6,  # 大容量
         marginal_cost=0.0)  # 零边际成本
```

### 新方式（节点边际电价）⭐
```python
# 虚拟发电机使用节点边际电价
if nodal_prices is not None and bus in nodal_prices.index:
    # 使用节点边际电价作为边际成本
    marginal_cost = nodal_prices[bus]
    logger.info(f"为节点 {bus} 添加虚拟发电机，使用节点边际电价")
else:
    # 如果没有节点电价，使用零边际成本
    marginal_cost = 0.0
    logger.info(f"为节点 {bus} 添加虚拟发电机，使用零边际成本（无节点电价数据）")
    
n_al.add("Generator",
         f"virtual_gen_{bus}",
         bus=bus,
         carrier="virtual",
         p_nom=1e6,  # 大容量
         marginal_cost=marginal_cost)  # 使用节点边际电价作为边际成本
```

### 节点电价提取
```python
# 提取当前迭代的节点电价（用于电解铝优化）
current_nodal_prices = None
if hasattr(n_current, 'buses_t') and hasattr(n_current.buses_t, 'marginal_price'):
    # 获取所有电力节点的边际电价
    electricity_buses = n_current.buses[n_current.buses.carrier != "aluminum"].index
    if len(electricity_buses) > 0:
        current_nodal_prices = n_current.buses_t.marginal_price[electricity_buses]
        logger.info(f"提取节点电价: {list(electricity_buses)}")
```

## 优势

1. **更稳定的收敛**：基于目标函数变化判断收敛，比电解铝用能变化更稳定
2. **更清洁的网络状态**：每次迭代重新加载网络，避免状态累积
3. **更准确的优化**：虚拟发电机使用节点边际电价，提高电解铝优化精度
4. **更简洁的实现**：使用p_set固定电解铝用能，代码更简洁高效
5. **更好的性能监控**：详细的时间统计帮助分析算法性能
6. **更灵活的配置**：支持通过配置文件调整算法参数
7. **更一致的实现**：与 `capacity_expansion_planning_aluminum_iterative.py` 保持一致 