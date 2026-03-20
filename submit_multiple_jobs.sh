#!/bin/bash
# Submit multiple independent SLURM jobs in batches，Based on snakemake dry-runDetermine whether to submit

echo "=== Submit multiple PyPSAs in batches-China SLURMOperation（Based on snakemake dry-runIntelligent detection）==="
echo "start time: $(date)"
echo "Notice: use --rerun-incomplete --ignore-incomplete --rerun-triggers mtime Parameters ignore configuration files and parameter changes"
echo

# Read the baseline version number（from config.yaml）
BASE_VERSION=$(grep -m 1 '^version:' config.yaml | sed 's/^version: //')
if [ -n "$BASE_VERSION" ]; then
    echo "Baseline version number: $BASE_VERSION"
else
    echo "⚠️  Failed to get from config.yamlRead the baseline version number，Version consistency check will be skipped"
fi
echo

# Check if the jobs folder exists
if [ ! -d "jobs" ]; then
    echo "✗ jobsFolder does not exist，Please run first generate_slurm_jobs_advanced.py Generate SLURM job files"
    exit 1
fi

# Check if the results folder exists
if [ ! -d "results" ]; then
    echo "✗ resultsFolder does not exist，All assignments will be submitted"
    CHECK_RESULTS=false
else
    CHECK_RESULTS=true
    echo "✓ Found the results folder，will be based on snakemake dry-runCheck completed configuration"
fi

# function：Parse configuration parameters from job file name
parse_job_filename() {
    local job_file="$1"
    local basename_job=$(basename "$job_file" .slurm)
    
    # Remove "job_" prefix
    local config_name=${basename_job#job_}
    
    # Parse configuration name，The format is like: HMH_2030_100p, LMH_2040_50p, NMM_2050_100p wait
    # Format: {scenario}_{year}_{capacity_ratio}
    if [[ "$config_name" =~ ^([A-Z]{3,4})_([0-9]{4})_(.+)$ ]]; then
        local scenario="${BASH_REMATCH[1]}"
        local year="${BASH_REMATCH[2]}"
        local capacity_ratio="${BASH_REMATCH[3]}"
        
        echo "$scenario|$year|$capacity_ratio"
    else
        echo ""
    fi
}

# function：Read version information from config file
get_config_version() {
    local scenario="$1"
    local year="$2"
    local capacity_ratio="$3"
    
    local config_file="configs/config_${scenario}_${year}_${capacity_ratio}.yaml"
    
    if [ -f "$config_file" ]; then
        # Read the version information of the first line
        local version=$(head -n 1 "$config_file" | sed 's/^version: //')
        if [ -n "$version" ]; then
            echo "$version"
            return 0
        fi
    fi
    
    echo ""
    return 1
}

# function：Extract base version number from full version（first one'-'forward）
get_base_version_from_full() {
    local full_version="$1"
    echo "${full_version%%-*}"
}

# function：Based on snakemake dry-runCheck if the result file exists
check_results_by_snakemake() {
    local scenario="$1"
    local year="$2"
    local capacity_ratio="$3"
    
    local config_file="configs/config_${scenario}_${year}_${capacity_ratio}.yaml"
    
    if [ ! -f "$config_file" ]; then
        return 1  # Configuration file does not exist，Think it needs to be dealt with
    fi
    
    # Use snakemake dry-runMode checks need to be run
    # -np means dry-run，Not actually implemented
    # --rerun-incomplete --ignore-incomplete Ignore configuration file changes
    # --rerun-triggers mtime Judgment based on file modification time only，Ignore parameter changes
    # If the output contains"Nothing to be done"，Indicates that all goals have been met
    local snakemake_output
    snakemake_output=$(snakemake --configfile "$config_file" -np --rerun-incomplete --ignore-incomplete --rerun-triggers mtime 2>/dev/null)
    
    if [ $? -eq 0 ] && echo "$snakemake_output" | grep -q "Nothing to be done"; then
        return 0  # All goals met，No need to run
    else
        return 1  # need to run
    fi
}

# Automatic discovery of all available SLURM job files
echo "Discovering available SLURM job files..."
JOBS=()
for job_file in jobs/job_*.slurm; do
    if [ -f "$job_file" ]; then
        JOBS+=("$job_file")
        echo "  ✓ Discover job files: $(basename "$job_file")"
    fi
done

if [ ${#JOBS[@]} -eq 0 ]; then
    echo "✗ No SLURM job files found，Please run first generate_slurm_jobs_advanced.py Generate job file"
    exit 1
fi

echo "Found in total ${#JOBS[@]} job files"
echo

# Check the results file status of each job
echo "Based on snakemake dry-runCheck completed configuration..."
PENDING_JOBS=()
COMPLETED_JOBS=()

for job_file in "${JOBS[@]}"; do
    config_info=$(parse_job_filename "$job_file")
    
    if [ -n "$config_info" ]; then
        IFS='|' read -r scenario year capacity_ratio <<< "$config_info"
        
        # Get the version information in the config file
        version=$(get_config_version "$scenario" "$year" "$capacity_ratio")
        
        if [ -n "$version" ]; then
            # Verify whether the version number is consistent with the baseline version
            if [ -n "$BASE_VERSION" ]; then
                version_base=$(get_base_version_from_full "$version")
                if [ "$version_base" != "$BASE_VERSION" ]; then
                    PENDING_JOBS+=("$job_file")
                    echo "  ⚠️  $(basename "$job_file") -> Version: $version (with baseline version $BASE_VERSION inconsistent，Mark as pending)"
                    continue
                fi
            fi

            echo "  📋 $(basename "$job_file") -> Version: $version"
            
            if [ "$CHECK_RESULTS" = true ] && check_results_by_snakemake "$scenario" "$year" "$capacity_ratio"; then
                COMPLETED_JOBS+=("$job_file")
                echo "    ✓ Completed (snakemake dry-runShows no need to run)"
            else
                PENDING_JOBS+=("$job_file")
                if [ "$CHECK_RESULTS" = true ]; then
                    echo "    ⏳ Pending (snakemake dry-runShows need to run)"
                else
                    echo "    ⏳ Pending (Results file not checked)"
                fi
            fi
        else
            # Unable to read version information，treat as pending
            PENDING_JOBS+=("$job_file")
            echo "  ⚠️  $(basename "$job_file") - Pending (Unable to read config version information)"
        fi
    else
        # Unresolved file name，treat as pending
        PENDING_JOBS+=("$job_file")
        echo "  ⚠️  $(basename "$job_file") - Pending (Unable to parse configuration parameters)"
    fi
done

echo
echo "Check results:"
echo "  Completed: ${#COMPLETED_JOBS[@]} configuration"
echo "  Pending: ${#PENDING_JOBS[@]} configuration"
echo

if [ ${#PENDING_JOBS[@]} -eq 0 ]; then
    echo "🎉 All configurations have been completed！No need to submit new assignments。"
    exit 0
fi

# Sort pending jobs by priority（non_flexiblepriority，then 100p，then in alphabetical order）
echo "Sort pending job files by priority..."
JOBS_SORTED=()

# Add non first_flexibleOperation（Baseline group）
for job_file in "${PENDING_JOBS[@]}"; do
    if [[ "$job_file" == *"non_flexible"* ]]; then
        JOBS_SORTED+=("$job_file")
    fi
done

# Add another 100p job
for job_file in "${PENDING_JOBS[@]}"; do
    if [[ "$job_file" == *"100p"* ]]; then
        JOBS_SORTED+=("$job_file")
    fi
done

# Add no_aluminumOperation
for job_file in "${PENDING_JOBS[@]}"; do
    if [[ "$job_file" == *"no_aluminum"* ]]; then
        JOBS_SORTED+=("$job_file")
    fi
done

# Finally add jobs with other capacity ratios
for job_file in "${PENDING_JOBS[@]}"; do
    if [[ "$job_file" != *"non_flexible"* ]] && [[ "$job_file" != *"100p"* ]] && [[ "$job_file" != *"no_aluminum"* ]]; then
        JOBS_SORTED+=("$job_file")
    fi
done

echo "Sorted order of pending jobs:"
for i in "${!JOBS_SORTED[@]}"; do
    echo "  $((i+1)). $(basename "${JOBS_SORTED[$i]}")"
done
echo

# Submit jobs in batches
echo "Start bulk submission of pending jobs..."
echo

declare -A JOB_IDS
SUBMITTED_COUNT=0
FAILED_COUNT=0

for job_file in "${JOBS_SORTED[@]}"; do
    echo "Submit assignment: $(basename "$job_file")"
    
    # Submit a job and get the job ID
    # Submit only pending jobs；Since it needs to be run，Just force you to run from the beginning
    FORCE_RESTART_FLAG=1
    JOB_ID=$(sbatch --export=ALL,FORCE_RESTART="$FORCE_RESTART_FLAG" "$job_file" | awk '{print $4}')
    
    if [ $? -eq 0 ] && [ -n "$JOB_ID" ]; then
        echo "  ✓ Submission successful，Job ID: $JOB_ID"
        JOB_IDS["$job_file"]="$JOB_ID"
        SUBMITTED_COUNT=$((SUBMITTED_COUNT + 1))
        
        # Wait a short period of time before submitting the next job，Avoid excessive system load
        sleep 2
    else
        echo "  ✗ Submission failed"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
    echo
done

# Show submission results
echo "=== Submit results summary ==="
echo "Submitted successfully: $SUBMITTED_COUNT homework"
echo "Submission failed: $FAILED_COUNT homework"
echo "skip completed: ${#COMPLETED_JOBS[@]} configuration"
echo

if [ $SUBMITTED_COUNT -gt 0 ]; then
    echo "Successfully submitted assignment:"
    for job_file in "${JOBS_SORTED[@]}"; do
        if [[ -n "${JOB_IDS[$job_file]}" ]]; then
            echo "  $(basename "$job_file") -> Job ID: ${JOB_IDS[$job_file]}"
        fi
    done
    echo
    
    # Save job information to file
    echo "Batch job information - $(date)" > "batch_jobs_info.txt"
    echo "========================" >> "batch_jobs_info.txt"
    echo "Job submission order:" >> "batch_jobs_info.txt"
    for job_file in "${JOBS_SORTED[@]}"; do
        if [[ -n "${JOB_IDS[$job_file]}" ]]; then
            echo "$(basename "$job_file") -> ${JOB_IDS[$job_file]}" >> "batch_jobs_info.txt"
        fi
    done
    
    echo "Job information has been saved to: batch_jobs_info.txt"
    echo
    
    # Show all job status
    echo "Current status of all jobs:"
    squeue -u $USER
    
    echo
    echo "Monitoring commands:"
    echo "  squeue -u $USER                    # View all user jobs"
    echo "  scancel <Job ID>                   # Cancel a specific job"
    echo "  tail -f slurm_<Job ID>.out         # View specific job output"
    echo "  tail -f slurm_<Job ID>.err         # View specific job errors"
    echo
    echo "Batch operation commands:"
    echo "  scancel \$(squeue -u $USER -h -o %i)  # Cancel all user jobs"
    echo "  watch -n 10 'squeue -u $USER'         # Monitor job status every 10 seconds"
fi

echo
echo "=== Batch submission completed ==="
echo "completion time: $(date)" 