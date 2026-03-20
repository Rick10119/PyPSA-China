# PyPSA-China SLURM Job File Guide

## Overview

This guide explains how to use the auto-generated SLURM job files to run different scenario simulations of PyPSA-China.

## Auto-generating SLURM Job Files

```bash
# Run from the project root directory
python scripts/generate_slurm_jobs_advanced.py
```

The advanced generator automatically discovers all available configuration files and produces the corresponding SLURM job files.

## Batch Submission

### Submit all scenarios

```bash
./submit_multiple_jobs.sh
```

This submits all scenario jobs in order:
1. Baseline scenario without aluminum smelters
2. 100 % capacity ratio
3. 90 % capacity ratio
4. 80 % capacity ratio
5. 70 % capacity ratio
6. 60 % capacity ratio
7. 55 % capacity ratio

### Submit capacity-ratio jobs only

```bash
./submit_capacity_jobs.sh
```

This submits only the capacity-ratio jobs (excluding the no-aluminum baseline).

### Submit a single job manually

```bash
# Submit the 100 % capacity-ratio job
sbatch job_100p.slurm

# Submit the no-aluminum job
sbatch job_no_aluminum.slurm
```

## Monitoring Job Status

### View all jobs

```bash
squeue -u $USER
```

### Real-time monitoring

```bash
watch -n 10 'squeue -u $USER'
```

### View output of a specific job

```bash
# Standard output
tail -f slurm_100p_<JOB_ID>.out

# Standard error
tail -f slurm_100p_<JOB_ID>.err
```

## Job Management

### Cancel a specific job

```bash
scancel <JOB_ID>
```

### Cancel all jobs for the current user

```bash
scancel $(squeue -u $USER -h -o %i)
```

### View job history

```bash
sacct -u $USER --format=JobID,JobName,State,Start,End,Elapsed
```

## Custom SLURM Jobs

### Generate a custom job with Python

```python
from scripts.generate_slurm_jobs_advanced import SlurmJobGenerator

# Create a generator instance
generator = SlurmJobGenerator()

# Generate a custom job file
generator.generate_custom_job(
    'test_scenario',           # Scenario name
    'config_test.yaml',        # Config file
    'Test scenario description',  # Description
    cpus_per_task=32,         # CPU cores
    time_limit='6:00:00',     # Time limit
    mem_per_cpu='20G'         # Memory per CPU
)
```

### Editing existing job files

Every generated SLURM job file can be edited manually. Common modifications include:

- Adjust CPU cores: `#SBATCH --cpus-per-task=32`
- Adjust memory: `#SBATCH --mem-per-cpu=20G`
- Adjust time limit: `#SBATCH --time=24:00:00`
- Set email notifications: `#SBATCH --mail-user=your.email@example.com`

## Output Files

### SLURM output files

- `slurm_<scenario>_<JOB_ID>.out` – standard output
- `slurm_<scenario>_<JOB_ID>.err` – standard error
- `job_<scenario>_<timestamp>.log` – job log

### Result files

After a simulation completes, results are saved under:
```
results/version-<version>-<scenario>/
```

## Troubleshooting

### Common issues

1. **Job submission fails**
   - Check the SLURM file syntax.
   - Verify that the config file exists and is well-formed.
   - Check user permissions and queue settings.

2. **Job run fails**
   - Inspect the error file: `slurm_<scenario>_<JOB_ID>.err`
   - Review parameters in the config file.
   - Make sure all required modules are loaded.

3. **Out of memory**
   - Increase memory: `#SBATCH --mem-per-cpu=20G`
   - Reduce CPU cores: `#SBATCH --cpus-per-task=20`

4. **Time limit exceeded**
   - Increase the time limit: `#SBATCH --time=24:00:00`
   - Tune simulation parameters for faster runs.

### Getting help

If you encounter problems:

1. Consult the SLURM documentation: `man sbatch`
2. Contact the system administrator.
3. Refer to the PyPSA-China project documentation.

## Best Practices

1. **Submission order**: run the baseline (no-aluminum) scenario first, then the capacity-ratio scenarios.
2. **Resource management**: adjust CPU cores and memory to match actual requirements.
3. **Monitoring**: check job status regularly and handle failures promptly.
4. **Backups**: keep copies of important config files and results.
5. **Logging**: retain job logs for debugging.
