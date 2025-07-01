# 电解铝迭代优化算法使用说明

## 概述

本文档介绍了PyPSA-China中电解铝迭代优化算法的实现和使用方法。该算法实现了您要求的四步迭代优化过程：

1. 使用连续化电解铝模型(aluminum_commitment: False)，求解模型，得到每个小时的节点电价
2. 基于节点电价，运行以满足铝需求为约束的电解铝最优运行问题(aluminum_commitment: True)，得到电解铝的用能
3. 固定电解铝的用能，求解剩下的优化问题
4. 重复2-3，直到收敛

## 算法原理

### 迭代过程

1. **初始化**: 设置电解铝启停约束为False，使用连续化模型
2. **网络求解**: 求解整个网络，获得节点电价
3. **电解铝优化**: 基于节点电价，单独优化电解铝运行模式
4. **固定用能**: 将电解铝用能固定，重新求解网络
5. **收敛检查**: 检查电解铝用能是否收敛
6. **迭代**: 重复步骤2-5直到收敛或达到最大迭代次数

### 收敛条件

算法在以下条件下停止：
- 电解铝用能的最大变化小于收敛容差（默认1e-6）
- 达到最大迭代次数（默认10次）

## 文件结构

### 新增文件

- `scripts/solve_network_aluminum_iterative.py`: 电解铝迭代优化求解脚本
- `config-aluminum-iterative.yaml`: 电解铝迭代优化配置文件
- `docs/aluminum_iterative_optimization.md`: 本文档

### 修改文件

- `Snakefile`: 添加了新的求解规则 `solve_network_aluminum_iterative`

## 使用方法

### 1. 配置参数

在配置文件中设置电解铝迭代优化参数：

```yaml
# 电解铝迭代优化参数
aluminum_iterative:
  max_iterations: 10  # 最大迭代次数
  convergence_tolerance: 1e-6  # 收敛容差
  enable_iterative_solving: True  # 是否启用迭代求解

# 电解铝基本参数
add_aluminum: True
aluminum_commitment: False  # 初始设置为False，在迭代过程中会动态调整
```

### 2. 运行优化

#### 方法1: 使用专用配置文件

```bash
# 使用电解铝迭代优化配置文件
snakemake --configfile config-aluminum-iterative.yaml solve_network_aluminum_iterative
```

#### 方法2: 在现有配置中启用

在您的配置文件中添加电解铝迭代优化参数，然后运行：

```bash
snakemake solve_network_aluminum_iterative
```

### 3. 输出结果

迭代优化会生成带有 `-iterative` 后缀的结果文件：

```
results/version-0628.24H.1/postnetworks/positive/postnetwork-ll-current+Neighbor-linear2050-2050-iterative.nc
```

## 配置参数详解

### aluminum_iterative 部分

- `max_iterations`: 最大迭代次数，默认10
- `convergence_tolerance`: 收敛容差，默认1e-6
- `enable_iterative_solving`: 是否启用迭代求解，默认True

### aluminum 部分

- `al_demand_ratio`: 铝需求比例，默认0.07（7%负荷）
- `al_excess_rate`: 铝过剩率，按年份设置
- `al_start_up_cost`: 电解铝启动成本
- `al_p_min_pu`: 电解铝最小出力比例
- `al_marginal_cost_storage`: 铝存储边际成本

## 算法特点

### 优势

1. **分离优化**: 将电解铝优化与网络优化分离，降低问题复杂度
2. **迭代收敛**: 通过迭代实现电解铝用能与节点电价的协调
3. **灵活配置**: 支持不同的收敛条件和迭代参数
4. **兼容性**: 与现有PyPSA-China框架完全兼容

### 注意事项

1. **计算时间**: 迭代算法会增加计算时间，特别是迭代次数较多时
2. **收敛性**: 算法可能不会完全收敛，需要根据实际情况调整参数
3. **内存使用**: 每次迭代都会创建网络副本，需要注意内存使用

## 调试和监控

### 日志输出

算法会在日志中输出详细的迭代信息：

```
开始电解铝迭代优化算法
开始第 1 次迭代
步骤1: 使用连续化电解铝模型求解
网络求解成功，获得节点电价
步骤2: 运行电解铝最优运行问题
电解铝优化问题求解成功
电解铝用能最大变化: 0.001234
第 1 次迭代完成
...
算法收敛，在第 3 次迭代后停止
电解铝迭代优化算法完成，共进行 3 次迭代
```

### 监控指标

- 迭代次数
- 电解铝用能变化量
- 求解状态和条件
- 收敛情况

## 故障排除

### 常见问题

1. **求解失败**: 检查求解器配置和网络参数
2. **不收敛**: 调整收敛容差或增加最大迭代次数
3. **内存不足**: 减少线程数或使用更小的网络

### 调试建议

1. 首先使用较小的网络测试算法
2. 检查电解铝相关参数设置
3. 监控日志输出，了解算法运行状态
4. 必要时调整求解器参数

## 扩展和定制

### 自定义约束

可以在 `extra_functionality` 函数中添加自定义约束：

```python
def add_custom_aluminum_constraints(n, fixed_aluminum_usage):
    # 添加自定义电解铝约束
    pass
```

### 修改收敛条件

可以修改收敛检查逻辑：

```python
# 在 solve_network_iterative 函数中修改
if change < convergence_tolerance:
    logger.info(f"算法收敛，在第 {iteration} 次迭代后停止")
    break
```

## 参考文献

1. PyPSA-China 项目文档
2. 电解铝行业技术规范
3. 电力系统优化理论

## 联系方式

如有问题或建议，请联系项目维护者。 