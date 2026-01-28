# Tutorial SLURM para Dinâmica Molecular (GROMACS/AMBER/MMPBSA)

O objetivo deste tutorial é explicar **como submeter e acompanhar jobs no SLURM** usando os arquivos já preparados no **CHARMM-GUI**.  
O foco aqui é **uso correto de CPU, GPU e memória**, além dos **comandos essenciais** do SLURM.

> **Contexto do seu ambiente**
- SLURM configurado como **single-node** (um servidor) chamado `node03`.
- `node03` tem **32 CPUs**, **128000 MB (~128 GB) de RAM** e **2 GPUs**.
- Partições:
  - `gromacs` → esperado **até 12 CPUs + 1 GPU**
  - `amber` → esperado **1 CPU + 1 GPU**
  - `MMPBSA` → **somente CPU**

---

## 1) Conceitos rápidos (o que você precisa entender)

### CPU: `--ntasks` vs `--cpus-per-task`
- `--ntasks` = número de **processos** (MPI).
- `--cpus-per-task` = número de **threads** por processo (OpenMP/pthreads).

No seu padrão:
- **GROMACS GPU**: 1 tarefa (`--ntasks=1`) com **12 threads** (`--cpus-per-task=12`), e o `mdrun` usa `-ntomp 12`.
- **AMBER GPU (pmemd.cuda)**: normalmente **1 tarefa e 1 CPU** (`--cpus-per-task=1`).
- **MMPBSA.py**: usa CPU e pode paralelizar (depende da instalação/flags); aqui vamos reservar, por exemplo, 8 CPUs.

### GPU: `--gres=gpu:X` e `CUDA_VISIBLE_DEVICES`
- Você pede GPU com:
  - `#SBATCH --gres=gpu:1`
- O SLURM então restringe quais GPUs o job pode ver via:
  - `CUDA_VISIBLE_DEVICES`
- **Boas práticas**:
  - Não force `-gpu_id 0` “no escuro”, porque no SLURM **GPU 0 dentro do job** pode ser *outra* GPU física.
  - Em geral, se o SLURM já setou `CUDA_VISIBLE_DEVICES`, você não precisa escolher manualmente.

### Memória: `--mem`
- `--mem=XXXX` define o máximo de RAM para o job.
- Se o job passar do limite, ele pode ser morto por OOM (out-of-memory).
- Como o nó tem `RealMemory=120000` MB, o total disponível é ~120 GB (esse valor foi setado justamente para deixar uma folga, sobrando assim 8GB para o sistema).

---

## 2) Comandos essenciais do SLURM (dia a dia)

### Ver partições e recursos disponíveis
```bash
sinfo
sinfo -N -l
````

### Ver sua fila de jobs

```bash
squeue
squeue -u $USER
```

### Submeter um job

```bash
sbatch meu_job.sbatch
```

### Cancelar um job

```bash
scancel <JOBID>
```

### Ver detalhes do job (muito útil para debug)

```bash
scontrol show job <JOBID>
```

### Ver histórico/uso de recursos após terminar (accounting)

```bash
sacct -j <JOBID> --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,AllocTRES%50
```

### Acompanhar o output

Pelo seu padrão de saída: `--output=slurm-%x-%j.out`

```bash
tail -f slurm-<NOME_DO_JOB>-<JOBID>.out
```

---

## 3) Estrutura recomendada do seu diretório de trabalho

Exemplo (GROMACS):

* `step3_input.gro`, `step4.0_minimization.mdp`, `step4.1_equilibration.mdp`, `step5_production.mdp`
* `topol.top`, `index.ndx`
* (saídas geradas pelo job: `.tpr`, `.gro`, `.cpt`, `.xtc`, `.edr`, logs)

Exemplo (AMBER):

* `step3_input.parm7`, `step3_input.rst7`
* `step4.0_minimization.mdin`, `step4.1_equilibration.mdin`, `step5_production.mdin`
* `dihe.restraint` (se aplicável), etc.

---

## 4) Regras práticas para pedir recursos (CPU/GPU/MEM)

### A) GROMACS (GPU) — padrão recomendado

* CPU: **12 threads**
* GPU: **1**
* Memória: normalmente 4–16 GB basta, mas depende do sistema (comece com **16G** e ajuste).

Exemplos de diretivas:

```bash
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
```

Dentro do script:

* Use `-ntomp ${SLURM_CPUS_PER_TASK}`
* Para GPU: `-nb gpu` (como você já faz)

### B) AMBER (pmemd.cuda)

* CPU: **1** (normal)
* GPU: **1**
* Memória: em geral **8G–32G** dependendo do sistema/tamanho do solvente.

### C) MMPBSA

* GPU: **nenhuma**
* CPU: 4–16 CPUs (conforme pressa vs fila)
* Memória: pode variar (trajetórias grandes podem exigir mais). Comece com **16G**.

---

## 5) Templates `.sbatch` prontos (GROMACS / AMBER / MMPBSA)

> **Importante**

* Use `#!/bin/bash`
* Adicione `set -euo pipefail` para falhar corretamente quando algo dá errado
* Sempre imprima informações úteis (hostname, data, variáveis SLURM, GPU visível)

---

# 5.1) GROMACS — Minimização

Crie `minim_gromacs.sbatch`:

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

# Pré-processamento
gmx grompp -f step4.0_minimization.mdp \
  -o step4.0_minimization.tpr \
  -c step3_input.gro -r step3_input.gro \
  -p topol.top -n index.ndx

# Execução
gmx mdrun -v -deffnm step4.0_minimization \
  -ntomp ${SLURM_CPUS_PER_TASK} -nb gpu
```

Submeter:

```bash
sbatch minim_gromacs.sbatch
```

---

# 5.2) GROMACS — Equilíbrio

Crie `equil_gromacs.sbatch`:

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

---

# 5.3) GROMACS — Produção (com restart automático por checkpoint)

Crie `prod_gromacs.sbatch`:

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

# Gere o .tpr (se não existir)
if [ ! -f step5_production.tpr ]; then
  gmx grompp -f step5_production.mdp \
    -o step5_production.tpr \
    -c step4.1_equilibration.gro \
    -p topol.top -n index.ndx
fi

# Se existir checkpoint, retoma; senão inicia do zero
if [ -f step5_production.cpt ]; then
  echo "Checkpoint encontrado: retomando..."
  gmx mdrun -v -deffnm step5_production \
    -ntomp ${SLURM_CPUS_PER_TASK} -cpi -nb gpu
else
  echo "Sem checkpoint: iniciando do zero..."
  gmx mdrun -v -deffnm step5_production \
    -ntomp ${SLURM_CPUS_PER_TASK} -nb gpu
fi
```

> **Por que isso é importante?**
> Se o job cair por tempo (`--time`) ou manutenção, o `.cpt` permite retomar sem perder tudo.

---

# 5.4) AMBER — Minimização (pmemd.cuda)

Crie `minim_amber.sbatch`:

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

# Se você usa um placeholder FC no arquivo de restrição:
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

---

# 5.5) AMBER — Equilíbrio (pmemd.cuda)

Crie `equil_amber.sbatch`:

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

---

# 5.6) AMBER — Produção (pmemd.cuda)

Crie `prod_amber.sbatch`:

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

---

# 5.7) MMPBSA (CPU)

## Atenção ao nome da partição

No seu `slurm.conf`, a partição chama **`MMPBSA`**:

```conf
PartitionName=MMPBSA Nodes=node03 ...
```

Então o correto é:

```bash
#SBATCH --partition=MMPBSA
```

Crie `mmpbsa.sbatch`:

```bash
#!/bin/bash
#SBATCH --job-name=mmpbsa
#SBATCH --partition=MMPBSA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time5-00:00:00
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

---

## 6) “Checklist” de debug (quando algo dá errado)

### Job não inicia / fica PENDING

* Ver motivo:

```bash
squeue -j <JOBID> -o "%.18i %.9P %.20j %.8u %.2t %.10M %.6D %R"
```

### Ver como o SLURM alocou recursos (CPU/GPU/MEM)

```bash
scontrol show job <JOBID>
```

### Confirmar GPU visível dentro do job

No `.sbatch`, já imprimimos:

```bash
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
```

Se quiser checar com `nvidia-smi` (se disponível no ambiente do job):

```bash
nvidia-smi
```

### Memória insuficiente (job morre “do nada”)

* Aumente `#SBATCH --mem=...` e consulte:

```bash
sacct -j <JOBID> --format=JobID,State,ReqMem,MaxRSS,ExitCode%20
```

---

## 7) Boas práticas para usuários (padrão “laboratório”)

* **Um job por etapa** (min, equil, prod) → facilita re-submissão e organização.
* Sempre use `--output=slurm-%x-%j.out` para não sobrescrever logs.
* Em produção longa, **use checkpoint** (GROMACS: `.cpt`) e tempo suficiente.
* Evite fixar GPU manualmente (`-gpu_id`), deixe o SLURM controlar via `CUDA_VISIBLE_DEVICES`.
* Ajuste `--mem` quando necessário (trajetórias grandes, sistemas grandes, MMPBSA pesado).

---

## 8) Como executar (resumo)

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

Monitorar:

```bash
squeue -u $USER
tail -f slurm-<jobname>-<jobid>.out
```
