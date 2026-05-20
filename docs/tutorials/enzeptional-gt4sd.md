# Enzyme Prediction and Optimization with Enzeptional (GT4SD)

This tutorial describes how to predict **feasibility** and **kcat** for enzymes using the **Enzeptional** framework from the **Generative Toolkit for Scientific Discovery (GT4SD)**.

The example focuses on polysaccharide hydrolysis reactions:

- exo-inulinase: inulin to fructose
- exo-levanase: levan to fructose
- invertase: sucrose to glucose + fructose
- endo-inulinase: inulin to inulooligosaccharides such as F4

## Technologies

- **GT4SD v1.5.0** with the Enzeptional framework for enzyme optimization
- **ESM2 650M** (`facebook/esm2_t33_650M_UR50D`) for protein embeddings
- **ChemBERTa** (`seyonec/ChemBERTa-zinc-base-v1`) for substrate/product representation
- **PyTorch 2.5.1** with GPU support, tested on RTX 4090
- Docker container execution with NVIDIA GPU support

## Requirements

- Docker
- Base image: `drugilsberg/gt4sd-base:v1.4.2-gpu`
- Working-directory files:
  - `all_mapping_active_sites.fasta`: enzyme sequences with mapped active sites
  - SMILES files such as `inulin_substrato.smi`, `produto_exo_inulin.smi`, and related reaction files
  - `run_pred.py`: adapted Enzeptional prediction script
  - `run_all_predictions.sh`: Bash script for multiple predictions

## Container Setup

Create the GPU-enabled container:

```bash
docker run --gpus all -dit --name gt4sd-gpu -v /home/node03/Alunos/lucas/gt4sd-core-gpu:/workspace drugilsberg/gt4sd-base:v1.4.2-gpu bash
```

To use a specific GPU, for example GPU device 1:

```bash
docker run --gpus '"device=1"' -dit --name gt4sd-gpu -v /home/node03/Alunos/lucas/gt4sd-core-gpu:/workspace drugilsberg/gt4sd-base:v1.4.2-gpu bash
```

Enter the container:

```bash
docker exec -it gt4sd-gpu bash
```

Run the setup commands inside the container:

```bash
# Install dependencies if requirements.txt is available
pip install -r requirements.txt

# Install GT4SD
pip install gt4sd==1.5.0 --upgrade

# Compatibility fix
pip uninstall numpy -y
pip cache purge
pip install numpy==1.23.5 --no-cache-dir --force-reinstall

# Refresh certificates if downloads fail
apt-get update --allow-unauthenticated
apt-get install -y ca-certificates
update-ca-certificates --fresh

# Upgrade PyTorch with CUDA support for RTX 4090
pip install torch==2.5.1 torchvision==0.20.1 torchaudio==2.5.1 --index-url https://download.pytorch.org/whl/cu118

# Check PyTorch and GPU visibility
python -c "import torch; print(f'Version: {torch.__version__} | CUDA: {torch.cuda.is_available()} | CUDA version: {torch.version.cuda}')"
```

The forced PyTorch upgrade enables GPU execution even when GT4SD brings older dependency constraints.

## Run Predictions

Make the batch script executable and run it:

```bash
chmod +x run_all_predictions.sh
./run_all_predictions.sh
```

## Example `run_all_predictions.sh`

```bash
#!/bin/bash

# =========================
# GENERAL SETTINGS
# =========================
FASTA="all_mapping_active_sites.fasta"
SCRIPT="run_pred.py"
PYTHON="python"

OUTDIR="results"
mkdir -p "${OUTDIR}"

# =========================
# EXECUTION FUNCTION
# =========================
run_prediction () {
    local TAG=$1
    local SUBSTRATE=$2
    local PRODUCT=$3

    local OUT_FEAS="${OUTDIR}/${TAG}_feasibility.csv"
    local OUT_KCAT="${OUTDIR}/${TAG}_kcat.csv"

    echo "===================================================="
    echo "   Running optimization: ${TAG}"
    echo "   Substrate: ${SUBSTRATE}"
    echo "   Product:   ${PRODUCT}"
    echo "===================================================="

    ${PYTHON} ${SCRIPT} \
        ${FASTA} \
        ${SUBSTRATE} \
        ${PRODUCT} \
        ${OUT_FEAS} \
        ${OUT_KCAT}

    if [ $? -ne 0 ]; then
        echo "Execution failed: ${TAG}"
    else
        echo "Finished successfully: ${TAG}"
        echo "   -> ${OUT_FEAS}"
        echo "   -> ${OUT_KCAT}"
    fi

    echo
}

# =========================
# RUNS
# =========================

# Inulin, exo reaction
run_prediction \
    "inulin_exo" \
    "inulin_substrato.smi" \
    "produto_exo_inulin.smi"

# Levan, exo reaction
run_prediction \
    "levan_exo" \
    "substrato-levan.smi" \
    "produto_levan_exo.smi"

# Invertase
run_prediction \
    "invertase" \
    "invertase_substrato.smi" \
    "invertase_produto.smi"

# Inulin to F4 + F4
run_prediction \
    "inulin_endo" \
    "inulin_substrato.smi" \
    "inulin_f4_produto.smi"

echo "All predictions finished."
```

## Outputs

The script runs predictions for four reactions and writes CSV files to `./results/`:

- `${TAG}_feasibility.csv`: catalytic feasibility score for each enzyme
- `${TAG}_kcat.csv`: predicted `kcat` values and evolutionary optimization results

## Runtime Notes

Runtime can range from hours to days depending on the number of sequences and evolutionary iterations. In the reference workflow, about 3,177 sequences were processed.

Use `nvidia-smi` to confirm GPU usage while the workflow is running.
