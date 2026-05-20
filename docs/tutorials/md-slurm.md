# SLURM for Molecular Simulation Workflows

This tutorial explains how to submit and monitor jobs with **SLURM** using files prepared by **CHARMM-GUI**. The focus is correct use of CPU, GPU, memory, and the essential SLURM commands for GROMACS, AMBER, and MMPBSA workflows.

## Example Environment

- SLURM is configured as a single-node server named `node03`.
- `node03` has 32 CPUs, about 128 GB RAM, and 2 GPUs.
- Partitions:
  - `gromacs`: expected up to 12 CPUs + 1 GPU
  - `amber`: expected 1 CPU + 1 GPU
  - `MMPBSA`: CPU only

Adjust partition names, CPU counts, GPU counts, and memory to match your actual cluster.

## 1. Quick Concepts

### CPU: `--ntasks` vs `--cpus-per-task`

- `--ntasks`: number of processes, usually MPI ranks
- `--cpus-per-task`: number of threads per process, such as OpenMP threads

Recommended patterns:

- **GROMACS GPU**: one task, `--ntasks=1`, with 12 threads, `--cpus-per-task=12`; use `mdrun -ntomp 12`
- **AMBER GPU with `pmemd.cuda`**: usually one task and one CPU
- **MMPBSA.py**: CPU-based and may be parallelized depending on installation and flags; reserve a reasonable CPU count, such as 8

### GPU: `--gres=gpu:X` and `CUDA_VISIBLE_DEVICES`

Request a GPU with:

```bash
#SBATCH --gres=gpu:1
```

SLURM restricts which GPUs the job can see through `CUDA_VISIBLE_DEVICES`.

Best practices:

- Avoid forcing `-gpu_id 0` blindly, because GPU 0 inside the job may not be physical GPU 0.
- If SLURM sets `CUDA_VISIBLE_DEVICES`, you usually do not need manual GPU selection.

### Memory: `--mem`

`--mem=XXXX` sets the maximum RAM for the job. If the job exceeds this limit, SLURM may kill it with an out-of-memory error.

## 2. Essential SLURM Commands

### Show Partitions and Resources

```bash
sinfo
sinfo -N -l
```

### Show Job Queue

```bash
squeue
squeue -u $USER
```

### Submit a Job

```bash
sbatch my_job.sbatch
```

### Cancel a Job

```bash
scancel <JOBID>
```

### Show Job Details

```bash
scontrol show job <JOBID>
```

### Show Resource Usage after Completion

```bash
sacct -j <JOBID> --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,AllocTRES%50
```

### Follow Job Output

If using `--output=slurm-%x-%j.out`:

```bash
tail -f slurm-<JOB_NAME>-<JOBID>.out
```

## 3. Recommended Working Directory Structure

Example for GROMACS:

- `step3_input.gro`
- `step4.0_minimization.mdp`
- `step4.1_equilibration.mdp`
- `step5_production.mdp`
- `topol.top`
- `index.ndx`
- generated outputs: `.tpr`, `.gro`, `.cpt`, `.xtc`, `.edr`, logs

Example for AMBER:

- `step3_input.parm7`
- `step3_input.rst7`
- `step4.0_minimization.mdin`
- `step4.1_equilibration.mdin`
- `step5_production.mdin`
- `dihe.restraint`, if applicable

## 4. Practical Resource Rules

### GROMACS GPU

Recommended starting point:

```bash
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
```

Inside the script, use:

```bash
-ntomp ${SLURM_CPUS_PER_TASK}
```

For GPU acceleration, use `-nb gpu`.

### AMBER `pmemd.cuda`

Recommended starting point:

- CPU: 1
- GPU: 1
- memory: 8-32 GB depending on system size

### MMPBSA

Recommended starting point:

- GPU: none
- CPU: 4-16 CPUs
- memory: start with 16 GB and adjust for large trajectories

## 5. Ready-to-Use `.sbatch` Templates

Use:

- `#!/bin/bash`
- `set -euo pipefail`
- useful debug prints: hostname, date, SLURM variables, and visible GPU
- unique job names such as `gmx_min_<user>` when multiple people use the same cluster

### 5.1 GROMACS Minimization

Create `minim_gromacs.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=gmx_min
#SBATCH --partition=gromacs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=5-00:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

echo "Job: ${SLURM_JOB_NAME}  ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-N/A}"
date

module purge 2>/dev/null || true

gmx grompp -f step4.0_minimization.mdp \
  -o step4.0_minimization.tpr \
  -c step3_input.gro -r step3_input.gro \
  -p topol.top -n index.ndx

gmx mdrun -v -deffnm step4.0_minimization \
  -ntomp ${SLURM_CPUS_PER_TASK} -nb gpu
```

Submit:

```bash
sbatch minim_gromacs.sbatch
```

### 5.2 GROMACS Equilibration

Create `equil_gromacs.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=gmx_equil
#SBATCH --partition=gromacs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

echo "Job: ${SLURM_JOB_NAME}  ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-N/A}"
date

module purge 2>/dev/null || true

gmx grompp -f step4.1_equilibration.mdp \
  -o step4.1_equilibration.tpr \
  -c step4.0_minimization.gro -r step3_input.gro \
  -p topol.top -n index.ndx

gmx mdrun -v -deffnm step4.1_equilibration \
  -ntomp ${SLURM_CPUS_PER_TASK} -nb gpu
```

### 5.3 GROMACS Production with Checkpoint Restart

Create `prod_gromacs.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=gmx_prod
#SBATCH --partition=gromacs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=5-00:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

echo "Job: ${SLURM_JOB_NAME}  ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-N/A}"
date

module purge 2>/dev/null || true

if [ ! -f step5_production.tpr ]; then
  gmx grompp -f step5_production.mdp \
    -o step5_production.tpr \
    -c step4.1_equilibration.gro \
    -p topol.top -n index.ndx
fi

if [ -f step5_production.cpt ]; then
  echo "Checkpoint found: resuming..."
  gmx mdrun -v -deffnm step5_production \
    -ntomp ${SLURM_CPUS_PER_TASK} -cpi -nb gpu
else
  echo "No checkpoint found: starting from scratch..."
  gmx mdrun -v -deffnm step5_production \
    -ntomp ${SLURM_CPUS_PER_TASK} -nb gpu
fi
```

Checkpoint restart avoids losing the full simulation if the job stops because of wall time or maintenance.

### 5.4 AMBER Minimization with `pmemd.cuda`

Create `minim_amber.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=amber_min
#SBATCH --partition=amber
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=05:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

echo "Job: ${SLURM_JOB_NAME}  ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-N/A}"
date

module purge 2>/dev/null || true

sed -e "s/FC/1.0/g" dihe.restraint > step4.0_minimization.rest

pmemd.cuda -O \
  -i step4.0_minimization.mdin \
  -p step3_input.parm7 \
  -c step3_input.rst7 \
  -o step4.0_minimization.mdout \
  -r step4.0_minimization.rst7 \
  -inf step4.0_minimization.mdinfo \
  -ref step3_input.rst7
```

### 5.5 AMBER Equilibration with `pmemd.cuda`

Create `equil_amber.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=amber_equil
#SBATCH --partition=amber
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

echo "Job: ${SLURM_JOB_NAME}  ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-N/A}"
date

module purge 2>/dev/null || true

sed -e "s/FC/1.0/g" dihe.restraint > step4.1_equilibration.rest

pmemd.cuda -O \
  -i step4.1_equilibration.mdin \
  -p step3_input.parm7 \
  -c step4.0_minimization.rst7 \
  -o step4.1_equilibration.mdout \
  -r step4.1_equilibration.rst7 \
  -inf step4.1_equilibration.mdinfo \
  -ref step3_input.rst7 \
  -x step4.1_equilibration.nc
```

### 5.6 AMBER Production with `pmemd.cuda`

Create `prod_amber.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=amber_prod
#SBATCH --partition=amber
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=5-00:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

echo "Job: ${SLURM_JOB_NAME}  ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-N/A}"
date

module purge 2>/dev/null || true

pmemd.cuda -O \
  -i step5_production.mdin \
  -p step3_input.parm7 \
  -c step4.1_equilibration.rst7 \
  -o step5_production.mdout \
  -r step5_production.rst7 \
  -inf step5_production.mdinfo \
  -x step5_production.nc
```

### 5.7 MMPBSA on CPU

If the partition is named `MMPBSA`, use:

```bash
#SBATCH --partition=MMPBSA
```

Create `mmpbsa.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=mmpbsa
#SBATCH --partition=MMPBSA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=5-00:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail

echo "Job: ${SLURM_JOB_NAME}  ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
date

module purge 2>/dev/null || true

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

MMPBSA.py -O \
  -i mmpbsa.in \
  -cp complex.parm7 \
  -rp receptor.parm7 \
  -lp ligand.parm7 \
  -y step5_production.nc \
  -o FINAL_RESULTS_MMPBSA.dat
```

## 6. Debug Checklist

### Job Does Not Start or Stays Pending

```bash
squeue -j <JOBID> -o "%.18i %.9P %.20j %.8u %.2t %.10M %.6D %R"
```

### Inspect Allocated Resources

```bash
scontrol show job <JOBID>
```

### Confirm Visible GPU

Inside the `.sbatch`, print:

```bash
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
```

If available in the job environment:

```bash
nvidia-smi
```

### Memory Problems

Increase `#SBATCH --mem=...` and inspect:

```bash
sacct -j <JOBID> --format=JobID,State,ReqMem,MaxRSS,ExitCode%20
```

## 7. Laboratory Practices

- Use one job per stage: minimization, equilibration, and production.
- Use `--output=slurm-%x-%j.out` to avoid overwriting logs.
- Use checkpointing for long production runs.
- Avoid manually forcing GPU IDs unless your cluster policy requires it.
- Adjust memory for large systems, large trajectories, and MMPBSA jobs.

## 8. Execution Summary

GROMACS:

```bash
sbatch minim_gromacs.sbatch
sbatch equil_gromacs.sbatch
sbatch prod_gromacs.sbatch
```

AMBER:

```bash
sbatch minim_amber.sbatch
sbatch equil_amber.sbatch
sbatch prod_amber.sbatch
```

MMPBSA:

```bash
sbatch mmpbsa.sbatch
```

Monitor:

```bash
squeue -u $USER
tail -f slurm-<jobname>-<jobid>.out
```
