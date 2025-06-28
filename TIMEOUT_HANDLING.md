# Gurobi 超时处理功能说明

## 问题描述

当 Gurobi 求解器达到时间限制时，会出现以下错误：
```
AttributeError: Unable to retrieve attribute 'x'. Did you mean: 'X'?
```

这是因为 Gurobi 在超时后，变量属性访问方式发生了变化，但 linopy 库没有正确处理这种情况。

## 解决方案

我们已经在 `scripts/solve_network_myopic.py` 中添加了专门的超时处理功能：

### 1. 超时检测和处理

- 自动检测 Gurobi 超时情况
- 尝试从 Gurobi 模型中提取次优解
- 处理属性访问错误

### 2. 三级求解策略

在 `config.yaml` 中配置了三个级别的求解器选项：

#### 默认选项 (default)
- 时间限制：8小时
- 标准精度参数
- 适用于大多数情况

#### 保守选项 (conservative)
- 时间限制：4小时
- 放宽的数值参数
- 当默认选项失败时使用

#### 超时回退选项 (timeout_fallback)
- 时间限制：2小时
- 进一步放宽的参数
- 最后的安全网

### 3. 次优解的使用

**重要：次优解是可以使用的！**

当 Gurobi 超时时，通常会返回一个次优解，这个解：
- 满足所有约束条件
- 虽然不是全局最优，但通常是很好的解
- 可以用于后续分析

## 使用方法

### 1. 配置时间限制

在 `config.yaml` 中调整时间限制：

```yaml
solving:
  solver_options:
    default:
      TimeLimit: 28800  # 8小时
    conservative:
      TimeLimit: 14400  # 4小时
    timeout_fallback:
      TimeLimit: 7200   # 2小时
```

### 2. 监控求解过程

启用输出标志来监控进度：

```yaml
solving:
  solver_options:
    default:
      OutputFlag: 1  # 启用输出
```

### 3. 处理结果

代码会自动处理不同的求解状态：

- `optimal`: 找到最优解
- `time_limit`: 超时但获得次优解
- `infeasible`: 问题无解

## 日志输出示例

```
INFO:检测到 Gurobi 超时后的属性访问错误
INFO:尝试从 Gurobi 模型中获取次优解...
INFO:Gurobi 模型状态: 9
INFO:Gurobi 返回了次优解，尝试提取解值...
INFO:成功提取了 1234 个变量的值
WARNING:求解达到时间限制，返回次优解
```

## 建议

1. **合理设置时间限制**：根据问题规模和计算资源调整
2. **监控求解进度**：使用 `OutputFlag: 1` 查看求解状态
3. **接受次优解**：在时间有限的情况下，次优解通常足够好
4. **分析解的质量**：检查目标函数值和约束违反情况

## 故障排除

如果仍然遇到问题：

1. 检查 Gurobi 许可证是否有效
2. 确认内存是否足够
3. 尝试进一步放宽数值参数
4. 考虑简化问题规模

## 测试

运行测试脚本验证功能：

```bash
python test_timeout_handling.py
```

这将创建一个简单的测试网络并验证超时处理功能是否正常工作。 