# 场景结果可视化指南

本指南介绍如何使用 `plot_scenario_comparison.py` 脚本来可视化 `generate_value_test_configs` 生成的场景结果。

## 概述

`plot_scenario_comparison.py` 脚本可以分析 `generate_value_test_configs` 生成的72个场景配置的结果，并生成：

1. **9个子图的对比图表**：每个子图对应一个market-demand组合（3×3），显示不同flexibility场景下的成本或容量变化
2. **摘要表格**：汇总所有场景的关键指标
3. **详细数据文件**：便于进一步分析

## 脚本功能

### 主要特性

- **基于配置文件**：从config.yaml读取基准版本号，构建完整的场景矩阵
- **智能场景构建**：自动构建所有demand-market组合的场景（3×3=9个）
- **容错数据处理**：如果某个场景数据缺失，会显示相应的提示信息
- **多维度对比**：支持costs和capacities两种文件类型的分析
- **智能分类**：将成本按技术类型和成本类型进行智能分类
- **可视化输出**：生成清晰的对比图表和表格
- **数据完整性统计**：提供详细的数据缺失统计信息
- **中英文支持**：支持中文显示和标签

### 场景维度

脚本处理的场景维度包括：

- **Demand（需求）**：L（Low）、M（Mid）、H（High）
- **Market（市场机会）**：L（Low）、M（Mid）、H（High）
- **配置类型**：100p（100%容量）、non_flexible（基准组）

脚本会自动构建3×3=9个场景组合，每个场景都包含100p和non_flexible两种配置的对比。

## 使用方法

### 基本用法

```bash
# 分析成本变化（默认）
python scripts/plot_scenario_comparison.py

# 分析容量变化
python scripts/plot_scenario_comparison.py --file-type capacities

# 指定结果目录
python scripts/plot_scenario_comparison.py --results-dir /path/to/results

# 指定输出目录
python scripts/plot_scenario_comparison.py --output /path/to/output

# 指定配置文件
python scripts/plot_scenario_comparison.py --config /path/to/config.yaml

# 详细输出
python scripts/plot_scenario_comparison.py --verbose
```

### 参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--results-dir` | `results` | 结果目录路径，包含所有version-*目录 |
| `--output` | `results/scenario_analysis` | 输出目录路径 |
| `--file-type` | `costs` | 分析的文件类型：`costs` 或 `capacities` |
| `--config` | `config.yaml` | 配置文件路径，用于获取基准版本号 |
| `--verbose` | False | 启用详细输出模式 |

## 输出结果

### 1. 场景对比图表

生成一个包含9个子图的大图，每个子图对应一个market-demand组合：

- **行**：Market级别（Low, Mid, High）
- **列**：Demand级别（Low, Mid, High）
- **内容**：显示100p配置相对于non_flexible配置的成本/容量变化

每个子图包含：
- 水平条形图显示各成本/容量分类的变化
- 绿色条形表示节约（正值）
- 红色条形表示增加（负值）
- 数值标签显示具体变化量（十亿人民币单位）

### 2. 摘要表格

生成CSV格式的摘要表格，包含：

- 场景代码和描述
- 总变化量（十亿人民币）
- 变化最大的前3个分类及其变化量

### 3. 目录结构

```
results/scenario_analysis/
├── scenario_plots/
│   └── scenario_comparison_costs.png
├── summary_tables/
│   └── scenario_summary_costs.csv
└── 其他分析文件...
```

### 4. 数据完整性统计

脚本会输出详细的数据完整性统计信息，包括：
- 总场景数和数据完整率
- 各场景组合的数据状态
- 缺失数据的类型和数量

## 数据要求

### 目录结构

脚本期望的结果目录结构：

```
results/
├── version-base_version-LLL-2050-100p/
│   └── summary/
│       └── postnetworks/
│           └── positive/
│               └── postnetwork-ll-current+Neighbor-linear2050-2050/
│                   ├── costs.csv
│                   └── capacities.csv
├── version-base_version-LLL-2050-non_flexible/
│   └── summary/
│       └── postnetworks/
│           └── positive/
│               └── postnetwork-ll-current+Neighbor-linear2050-2050/
│                   ├── costs.csv
│                   └── capacities.csv
└── ... (其他场景)
```

### 文件格式

- **costs.csv**：成本数据，包含多级索引（component_type, cost_type, carrier）
- **capacities.csv**：容量数据，格式类似

## 成本分类说明

脚本将成本按以下分类进行组织：

1. **variable cost-non-renewable**：非可再生能源可变成本
2. **capital-non-renewable**：非可再生能源资本成本
3. **capital–demand side**：需求侧资本成本
4. **capital–renewable**：可再生能源资本成本
5. **transmission lines**：输电线路
6. **batteries**：电池储能
7. **long-duration storages**：长时储能
8. **carbon capture**：碳捕获
9. **synthetic fuels**：合成燃料
10. **carbon management**：碳管理

## 使用示例

### 示例1：分析所有场景的成本变化

```bash
cd PyPSA-China
python scripts/plot_scenario_comparison.py --verbose
```

### 示例2：使用自定义配置文件

```bash
python scripts/plot_scenario_comparison.py \
    --config configs/config_HHH_2050_100p.yaml \
    --verbose
```

### 示例3：分析容量变化并指定输出目录

```bash
python scripts/plot_scenario_comparison.py \
    --file-type capacities \
    --output results/capacity_analysis \
    --verbose
```

### 示例4：分析特定结果目录

```bash
python scripts/plot_scenario_comparison.py \
    --results-dir /path/to/specific/results \
    --output /path/to/output \
    --verbose
```

## 故障排除

### 常见问题

1. **没有找到场景结果**
   - 检查结果目录是否存在
   - 确认版本目录命名格式正确
   - 验证summary文件是否完整

2. **数据加载失败**
   - 检查CSV文件格式
   - 确认文件路径正确
   - 查看详细错误日志

3. **图表生成失败**
   - 检查matplotlib配置
   - 确认有足够的内存
   - 验证输出目录权限

### 调试模式

使用 `--verbose` 参数启用详细输出，可以：

- 查看每个场景的数据加载状态
- 监控成本分类过程
- 检查数据转换和计算过程

## 扩展功能

### 自定义成本分类

可以修改脚本中的 `cost_category_mapping` 字典来自定义成本分类：

```python
cost_category_mapping = {
    ('marginal', 'coal'): '自定义分类名称',
    # 添加更多映射...
}
```

### 添加新的分析维度

可以扩展脚本来分析其他维度，如：

- 时间序列分析
- 地理分布分析
- 技术组合分析

## 注意事项

1. **内存使用**：处理大量场景时可能需要较大内存
2. **运行时间**：首次运行可能需要较长时间来加载所有数据
3. **文件权限**：确保有读取结果目录和写入输出目录的权限
4. **数据完整性**：建议在运行脚本前确保所有场景都已完成计算

## 联系支持

如果遇到问题或需要帮助，请：

1. 查看详细错误日志（使用 `--verbose` 参数）
2. 检查数据文件格式和完整性
3. 参考现有的示例和文档
4. 联系开发团队获取支持
