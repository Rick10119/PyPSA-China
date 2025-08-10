# PyPSA-China SLURM作业文件使用指南

## 概述

本指南介绍如何使用自动生成的SLURM作业文件来运行PyPSA-China的不同场景模拟。

## 自动生成SLURM作业文件

### 方法1：使用基础生成器

```bash
# 在项目根目录下运行
python scripts/generate_slurm_jobs.py
```

这将生成以下SLURM作业文件：
- `job_no_aluminum.slurm` - 不包含电解铝厂的基准场景
- `job_100p.slurm` - 100%容量比例
- `job_90p.slurm` - 90%容量比例
- `job_80p.slurm` - 80%容量比例
- `job_70p.slurm` - 70%容量比例
- `job_60p.slurm` - 60%容量比例
- `job_55p.slurm` - 55%容量比例

### 方法2：使用高级生成器

```bash
# 在项目根目录下运行
python scripts/generate_slurm_jobs_advanced.py
```

高级生成器会自动发现所有可用的配置文件，并生成对应的SLURM作业文件。

## 批量提交作业

### 提交所有场景

```bash
./submit_multiple_jobs.sh
```

这将按顺序提交所有场景的作业：
1. 不包含电解铝厂的基准场景
2. 100%容量比例
3. 90%容量比例
4. 80%容量比例
5. 70%容量比例
6. 60%容量比例
7. 55%容量比例

### 只提交容量比例作业

```bash
./submit_capacity_jobs.sh
```

这将只提交容量比例相关的作业（不包括不包含电解铝厂的场景）。

### 手动提交单个作业

```bash
# 提交100%容量比例作业
sbatch job_100p.slurm

# 提交不包含电解铝厂的作业
sbatch job_no_aluminum.slurm
```

## 监控作业状态

### 查看所有作业状态

```bash
squeue -u $USER
```

### 实时监控作业状态

```bash
watch -n 10 'squeue -u $USER'
```

### 查看特定作业输出

```bash
# 查看标准输出
tail -f slurm_100p_<作业ID>.out

# 查看错误输出
tail -f slurm_100p_<作业ID>.err
```

## 作业管理

### 取消特定作业

```bash
scancel <作业ID>
```

### 取消所有用户作业

```bash
scancel $(squeue -u $USER -h -o %i)
```

### 查看作业历史

```bash
sacct -u $USER --format=JobID,JobName,State,Start,End,Elapsed
```

## 自定义SLURM作业

### 使用Python脚本生成自定义作业

```python
from scripts.generate_slurm_jobs_advanced import SlurmJobGenerator

# 创建生成器实例
generator = SlurmJobGenerator()

# 生成自定义参数的作业文件
generator.generate_custom_job(
    'test_scenario',           # 场景名称
    'config_test.yaml',        # 配置文件
    '测试场景描述',            # 描述
    cpus_per_task=32,         # CPU核心数
    time_limit='6:00:00',     # 时间限制
    mem_per_cpu='20G'         # 内存设置
)
```

### 修改现有作业文件

每个生成的SLURM作业文件都可以手动编辑，常见的修改包括：

- 调整CPU核心数：`#SBATCH --cpus-per-task=32`
- 调整内存设置：`#SBATCH --mem-per-cpu=20G`
- 调整时间限制：`#SBATCH --time=24:00:00`
- 修改邮件通知：`#SBATCH --mail-user=your.email@example.com`

## 输出文件说明

### SLURM输出文件

- `slurm_<场景名>_<作业ID>.out` - 标准输出
- `slurm_<场景名>_<作业ID>.err` - 标准错误
- `job_<场景名>_<时间戳>.log` - 作业日志

### 结果文件

模拟完成后，结果文件将保存在：
```
results/version-<版本号>-<场景名>/
```

## 故障排除

### 常见问题

1. **作业提交失败**
   - 检查SLURM作业文件语法是否正确
   - 确认配置文件存在且格式正确
   - 检查用户权限和队列设置

2. **作业运行失败**
   - 查看错误输出文件：`slurm_<场景名>_<作业ID>.err`
   - 检查配置文件中的参数设置
   - 确认所有依赖模块已正确加载

3. **内存不足**
   - 增加内存设置：`#SBATCH --mem-per-cpu=20G`
   - 减少CPU核心数：`#SBATCH --cpus-per-task=20`

4. **时间超限**
   - 增加时间限制：`#SBATCH --time=24:00:00`
   - 优化模拟参数以提高运行效率

### 获取帮助

如果遇到问题，可以：

1. 查看SLURM文档：`man sbatch`
2. 联系系统管理员
3. 检查PyPSA-China项目文档

## 最佳实践

1. **作业提交顺序**：建议先运行基准场景（不包含电解铝厂），再运行其他容量比例场景
2. **资源管理**：根据实际需求调整CPU核心数和内存设置
3. **监控**：定期检查作业状态，及时处理失败的作业
4. **备份**：保存重要的配置文件和结果文件
5. **日志**：保留作业日志以便问题排查 