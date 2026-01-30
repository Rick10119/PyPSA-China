#!/bin/bash
# 批量提交多个独立的SLURM作业，基于snakemake dry-run判断是否提交

echo "=== 批量提交多个PyPSA-China SLURM作业（基于snakemake dry-run智能检测）==="
echo "开始时间: $(date)"
echo "注意: 使用 --rerun-incomplete --ignore-incomplete --rerun-triggers mtime 参数忽略配置文件和参数变化"
echo

# 读取基准版本号（来自config.yaml）
BASE_VERSION=$(grep -m 1 '^version:' config.yaml | sed 's/^version: //')
if [ -n "$BASE_VERSION" ]; then
    echo "基准版本号: $BASE_VERSION"
else
    echo "⚠️  未能从config.yaml读取基准版本号，将跳过版本一致性校验"
fi
echo

# 检查jobs文件夹是否存在
if [ ! -d "jobs" ]; then
    echo "✗ jobs文件夹不存在，请先运行 generate_slurm_jobs_advanced.py 生成SLURM作业文件"
    exit 1
fi

# 检查results文件夹是否存在
if [ ! -d "results" ]; then
    echo "✗ results文件夹不存在，将提交所有作业"
    CHECK_RESULTS=false
else
    CHECK_RESULTS=true
    echo "✓ 发现results文件夹，将基于snakemake dry-run检查已完成的配置"
fi

# 函数：从job文件名解析配置参数
parse_job_filename() {
    local job_file="$1"
    local basename_job=$(basename "$job_file" .slurm)
    
    # 移除 "job_" 前缀
    local config_name=${basename_job#job_}
    
    # 解析配置名称，格式如: HMH_2030_100p, LMH_2040_50p, NMM_2050_100p 等
    # 格式: {scenario}_{year}_{capacity_ratio}
    if [[ "$config_name" =~ ^([A-Z]{3,4})_([0-9]{4})_(.+)$ ]]; then
        local scenario="${BASH_REMATCH[1]}"
        local year="${BASH_REMATCH[2]}"
        local capacity_ratio="${BASH_REMATCH[3]}"
        
        echo "$scenario|$year|$capacity_ratio"
    else
        echo ""
    fi
}

# 函数：从config文件读取version信息
get_config_version() {
    local scenario="$1"
    local year="$2"
    local capacity_ratio="$3"
    
    local config_file="configs/config_${scenario}_${year}_${capacity_ratio}.yaml"
    
    if [ -f "$config_file" ]; then
        # 读取第一行的version信息
        local version=$(head -n 1 "$config_file" | sed 's/^version: //')
        if [ -n "$version" ]; then
            echo "$version"
            return 0
        fi
    fi
    
    echo ""
    return 1
}

# 函数：从完整version中提取基准版本号（第一个'-'前）
get_base_version_from_full() {
    local full_version="$1"
    echo "${full_version%%-*}"
}

# 函数：基于snakemake dry-run检查结果文件是否存在
check_results_by_snakemake() {
    local scenario="$1"
    local year="$2"
    local capacity_ratio="$3"
    
    local config_file="configs/config_${scenario}_${year}_${capacity_ratio}.yaml"
    
    if [ ! -f "$config_file" ]; then
        return 1  # 配置文件不存在，认为需要处理
    fi
    
    # 使用snakemake dry-run模式检查是否需要运行
    # -np 表示dry-run，不实际执行
    # --rerun-incomplete --ignore-incomplete 忽略配置文件变化
    # --rerun-triggers mtime 只基于文件修改时间判断，忽略参数变化
    # 如果输出包含"Nothing to be done"，说明所有目标都已满足
    local snakemake_output
    snakemake_output=$(snakemake --configfile "$config_file" -np --rerun-incomplete --ignore-incomplete --rerun-triggers mtime 2>/dev/null)
    
    if [ $? -eq 0 ] && echo "$snakemake_output" | grep -q "Nothing to be done"; then
        return 0  # 所有目标都已满足，不需要运行
    else
        return 1  # 需要运行
    fi
}

# 自动发现所有可用的SLURM作业文件
echo "正在发现可用的SLURM作业文件..."
JOBS=()
for job_file in jobs/job_*.slurm; do
    if [ -f "$job_file" ]; then
        JOBS+=("$job_file")
        echo "  ✓ 发现作业文件: $(basename "$job_file")"
    fi
done

if [ ${#JOBS[@]} -eq 0 ]; then
    echo "✗ 未发现任何SLURM作业文件，请先运行 generate_slurm_jobs_advanced.py 生成作业文件"
    exit 1
fi

echo "共发现 ${#JOBS[@]} 个作业文件"
echo

# 检查每个作业的结果文件状态
echo "基于snakemake dry-run检查已完成的配置..."
PENDING_JOBS=()
COMPLETED_JOBS=()

for job_file in "${JOBS[@]}"; do
    config_info=$(parse_job_filename "$job_file")
    
    if [ -n "$config_info" ]; then
        IFS='|' read -r scenario year capacity_ratio <<< "$config_info"
        
        # 获取config文件中的version信息
        version=$(get_config_version "$scenario" "$year" "$capacity_ratio")
        
        if [ -n "$version" ]; then
            # 校验版本号是否与基准版本一致
            if [ -n "$BASE_VERSION" ]; then
                version_base=$(get_base_version_from_full "$version")
                if [ "$version_base" != "$BASE_VERSION" ]; then
                    PENDING_JOBS+=("$job_file")
                    echo "  ⚠️  $(basename "$job_file") -> 版本: $version (与基准版本 $BASE_VERSION 不一致，标记为待处理)"
                    continue
                fi
            fi

            echo "  📋 $(basename "$job_file") -> 版本: $version"
            
            if [ "$CHECK_RESULTS" = true ] && check_results_by_snakemake "$scenario" "$year" "$capacity_ratio"; then
                COMPLETED_JOBS+=("$job_file")
                echo "    ✓ 已完成 (snakemake dry-run显示无需运行)"
            else
                PENDING_JOBS+=("$job_file")
                if [ "$CHECK_RESULTS" = true ]; then
                    echo "    ⏳ 待处理 (snakemake dry-run显示需要运行)"
                else
                    echo "    ⏳ 待处理 (未检查结果文件)"
                fi
            fi
        else
            # 无法读取version信息，当作待处理
            PENDING_JOBS+=("$job_file")
            echo "  ⚠️  $(basename "$job_file") - 待处理 (无法读取config版本信息)"
        fi
    else
        # 无法解析的文件名，当作待处理
        PENDING_JOBS+=("$job_file")
        echo "  ⚠️  $(basename "$job_file") - 待处理 (无法解析配置参数)"
    fi
done

echo
echo "检查结果:"
echo "  已完成: ${#COMPLETED_JOBS[@]} 个配置"
echo "  待处理: ${#PENDING_JOBS[@]} 个配置"
echo

if [ ${#PENDING_JOBS[@]} -eq 0 ]; then
    echo "🎉 所有配置都已完成！无需提交新作业。"
    exit 0
fi

# 按优先级排序待处理的作业（non_flexible优先，然后100p，然后按字母顺序）
echo "按优先级排序待处理的作业文件..."
JOBS_SORTED=()

# 先添加non_flexible作业（基准组）
for job_file in "${PENDING_JOBS[@]}"; do
    if [[ "$job_file" == *"non_flexible"* ]]; then
        JOBS_SORTED+=("$job_file")
    fi
done

# 再添加100p作业
for job_file in "${PENDING_JOBS[@]}"; do
    if [[ "$job_file" == *"100p"* ]]; then
        JOBS_SORTED+=("$job_file")
    fi
done

# 再添加no_aluminum作业
for job_file in "${PENDING_JOBS[@]}"; do
    if [[ "$job_file" == *"no_aluminum"* ]]; then
        JOBS_SORTED+=("$job_file")
    fi
done

# 最后添加其他容量比例的作业
for job_file in "${PENDING_JOBS[@]}"; do
    if [[ "$job_file" != *"non_flexible"* ]] && [[ "$job_file" != *"100p"* ]] && [[ "$job_file" != *"no_aluminum"* ]]; then
        JOBS_SORTED+=("$job_file")
    fi
done

echo "排序后的待处理作业顺序:"
for i in "${!JOBS_SORTED[@]}"; do
    echo "  $((i+1)). $(basename "${JOBS_SORTED[$i]}")"
done
echo

# 批量提交作业
echo "开始批量提交待处理的作业..."
echo

declare -A JOB_IDS
SUBMITTED_COUNT=0
FAILED_COUNT=0

for job_file in "${JOBS_SORTED[@]}"; do
    echo "提交作业: $(basename "$job_file")"
    
    # 提交作业并获取作业ID
    # 只提交待处理作业；既然需要运行，就强制从头跑
    FORCE_RESTART_FLAG=0
    JOB_ID=$(sbatch --export=ALL,FORCE_RESTART="$FORCE_RESTART_FLAG" "$job_file" | awk '{print $4}')
    
    if [ $? -eq 0 ] && [ -n "$JOB_ID" ]; then
        echo "  ✓ 提交成功，作业ID: $JOB_ID"
        JOB_IDS["$job_file"]="$JOB_ID"
        SUBMITTED_COUNT=$((SUBMITTED_COUNT + 1))
        
        # 等待一小段时间再提交下一个作业，避免系统负载过高
        sleep 2
    else
        echo "  ✗ 提交失败"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
    echo
done

# 显示提交结果
echo "=== 提交结果汇总 ==="
echo "成功提交: $SUBMITTED_COUNT 个作业"
echo "提交失败: $FAILED_COUNT 个作业"
echo "跳过已完成: ${#COMPLETED_JOBS[@]} 个配置"
echo

if [ $SUBMITTED_COUNT -gt 0 ]; then
    echo "成功提交的作业:"
    for job_file in "${JOBS_SORTED[@]}"; do
        if [[ -n "${JOB_IDS[$job_file]}" ]]; then
            echo "  $(basename "$job_file") -> 作业ID: ${JOB_IDS[$job_file]}"
        fi
    done
    echo
    
    # 保存作业信息到文件
    echo "批量作业信息 - $(date)" > "batch_jobs_info.txt"
    echo "========================" >> "batch_jobs_info.txt"
    echo "作业提交顺序:" >> "batch_jobs_info.txt"
    for job_file in "${JOBS_SORTED[@]}"; do
        if [[ -n "${JOB_IDS[$job_file]}" ]]; then
            echo "$(basename "$job_file") -> ${JOB_IDS[$job_file]}" >> "batch_jobs_info.txt"
        fi
    done
    
    echo "作业信息已保存到: batch_jobs_info.txt"
    echo
    
    # 显示所有作业状态
    echo "当前所有作业状态:"
    squeue -u $USER
    
    echo
    echo "监控命令:"
    echo "  squeue -u $USER                    # 查看用户所有作业"
    echo "  scancel <作业ID>                   # 取消特定作业"
    echo "  tail -f slurm_<作业ID>.out         # 查看特定作业输出"
    echo "  tail -f slurm_<作业ID>.err         # 查看特定作业错误"
    echo
    echo "批量操作命令:"
    echo "  scancel \$(squeue -u $USER -h -o %i)  # 取消所有用户作业"
    echo "  watch -n 10 'squeue -u $USER'         # 每10秒监控作业状态"
fi

echo
echo "=== 批量提交完成 ==="
echo "完成时间: $(date)" 