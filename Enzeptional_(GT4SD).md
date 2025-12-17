# Predição e Otimização de Enzimas para Reações Específicas com Enzeptional (GT4SD)

Este README descreve a solução para predizer **feasibility** (viabilidade catalítica) e **kcat** (turnover number) de enzimas usando o framework **Enzeptional** do **Generative Toolkit for Scientific Discovery (GT4SD)**.

O foco está em reações de hidrólise de polissacarídeos:
- Exo-inulinase (inulin → frutose)
- Exo-levanase (levan → frutose)
- Invertase (sacarose → glicose + frutose)
- Endo-inulinase (inulin → inulooligossacarídeos F4)

## Tecnologias Utilizadas
- **GT4SD v1.5.0** com framework Enzeptional para otimização de enzimas.
- Modelos principais:
  - ESM2 650M (facebook/esm2_t33_650M_UR50D) para embeddings de proteínas.
  - ChemBERTa (seyonec/ChemBERTa-zinc-base-v1) para representação de substratos/produtos.
- PyTorch (2.5.1) com suporte a **GPU (RTX 4090)** para aceleração.
- Execução em container Docker com NVIDIA GPU support.

## Requisitos
- Docker
- Imagem base: `drugilsberg/gt4sd-base:v1.4.2-gpu`.
- Arquivos necessários na pasta de trabalho:
  - `all_mapping_active_sites.fasta` (sequências de enzimas com sítios ativos mapeados).
  - Arquivos SMILES: `inulin_substrato.smi`, `produto_exo_inulin.smi`, etc.
  - `run_pred.py` (script principal do Enzeptional, adaptado).
  - `run_all_predictions.sh` (script bash para executar múltiplas predições).

## Instalação e Configuração (dentro do Container)

Para criar o container, execute:

```bash
docker run --gpus all -dit --name gt4sd-gpu -v /home/node03/Alunos/lucas/gt4sd-core-gpu:/workspace drugilsberg/gt4sd-base:v1.4.2-gpu bash
```

Acesse o container:

```bash
docker exec -it gt4sd-gpu bash
```

Execute esses comandos no container Docker:

```bash

# Instale dependências (se houver requirements.txt)
pip install -r requirements.txt

# Instale GT4SD
pip install gt4sd==1.5.0 --upgrade

# Correções para compatibilidade
pip uninstall numpy -y
pip cache purge
pip install numpy==1.23.5 --no-cache-dir --force-reinstall

# Atualize certificados (opcional, para downloads)
apt-get update --allow-unauthenticated
apt-get install -y ca-certificates
update-ca-certificates --fresh

# Upgrade para PyTorch com suporte a CUDA (RTX 4090)
pip install torch==2.5.1 torchvision==0.20.1 torchaudio==2.5.1 --index-url https://download.pytorch.org/whl/cu118

# Verifique PyTorch + GPU
python -c "import torch; print(f'Versão: {torch.__version__} | CUDA: {torch.cuda.is_available()} | Versão CUDA: {torch.version.cuda}')"


**Nota:** O upgrade forçado do PyTorch permite uso da GPU, mesmo com dependências antigas do GT4SD.

## Execução
Torne o script executável e rode:

```bash
chmod +x run_all_predictions.sh
./run_all_predictions.sh
```

**run_all_predictions.sh**
```bash
#!/bin/bash

# =========================
# CONFIGURAÇÕES GERAIS
# =========================
FASTA="all_mapping_active_sites.fasta"
SCRIPT="run_pred.py"
PYTHON="python"

# Pasta para resultados
OUTDIR="results"
mkdir -p "${OUTDIR}"

# =========================
# FUNÇÃO DE EXECUÇÃO
# =========================
run_prediction () {
    local TAG=$1
    local SUBSTRATE=$2
    local PRODUCT=$3

    local OUT_FEAS="${OUTDIR}/${TAG}_feasibility.csv"
    local OUT_KCAT="${OUTDIR}/${TAG}_kcat.csv"

    echo "===================================================="
    echo "   Executando otimização: ${TAG}"
    echo "   Substrato: ${SUBSTRATE}"
    echo "   Produto:   ${PRODUCT}"
    echo "===================================================="

    ${PYTHON} ${SCRIPT} \
        ${FASTA} \
        ${SUBSTRATE} \
        ${PRODUCT} \
        ${OUT_FEAS} \
        ${OUT_KCAT}

    if [ $? -ne 0 ]; then
        echo "Erro na execução: ${TAG}"
    else
        echo "Finalizado com sucesso: ${TAG}"
        echo "   → ${OUT_FEAS}"
        echo "   → ${OUT_KCAT}"
    fi

    echo
}

# =========================
# EXECUÇÕES
# =========================

# Inulin (exo)
run_prediction \
    "inulin_exo" \
    "inulin_substrato.smi" \
    "produto_exo_inulin.smi"

# Levan (exo)
run_prediction \
    "levan_exo" \
    "substrato-levan.smi" \
    "produto_levan_exo.smi"

# Invertase
run_prediction \
    "invertase" \
    "invertase_substrato.smi" \
    "invertase_produto.smi"

# Inulin → F4+F4
run_prediction \
    "inulin_endo" \
    "inulin_substrato.smi" \
    "inulin_f4_produto.smi"

echo "Todas as predições foram executadas!"#!/bin/bash

# =========================
# CONFIGURAÇÕES GERAIS
# =========================
FASTA="all_mapping_active_sites.fasta"
SCRIPT="run_pred.py"
PYTHON="python"

# Pasta para resultados
OUTDIR="results"
mkdir -p "${OUTDIR}"

# =========================
# FUNÇÃO DE EXECUÇÃO
# =========================
run_prediction () {
    local TAG=$1
    local SUBSTRATE=$2
    local PRODUCT=$3

    local OUT_FEAS="${OUTDIR}/${TAG}_feasibility.csv"
    local OUT_KCAT="${OUTDIR}/${TAG}_kcat.csv"

    echo "===================================================="
    echo "   Executando otimização: ${TAG}"
    echo "   Substrato: ${SUBSTRATE}"
    echo "   Produto:   ${PRODUCT}"
    echo "===================================================="

    ${PYTHON} ${SCRIPT} \
        ${FASTA} \
        ${SUBSTRATE} \
        ${PRODUCT} \
        ${OUT_FEAS} \
        ${OUT_KCAT}

    if [ $? -ne 0 ]; then
        echo "Erro na execução: ${TAG}"
    else
        echo "Finalizado com sucesso: ${TAG}"
        echo "   → ${OUT_FEAS}"
        echo "   → ${OUT_KCAT}"
    fi

    echo
}

# =========================
# EXECUÇÕES
# =========================

# Inulin (exo)
run_prediction \
    "inulin_exo" \
    "inulin_substrato.smi" \
    "produto_exo_inulin.smi"

# Levan (exo)
run_prediction \
    "levan_exo" \
    "substrato-levan.smi" \
    "produto_levan_exo.smi"

# Invertase
run_prediction \
    "invertase" \
    "invertase_substrato.smi" \
    "invertase_produto.smi"

# Inulin → F4+F4
run_prediction \
    "inulin_endo" \
    "inulin_substrato.smi" \
    "inulin_f4_produto.smi"

echo "Todas as predições foram executadas!"
```

Oript executa predições para as 4 reações, gerando arquivos CSV em `./results/`:
- `${TAG}_feasibility.csv`: Viabilidade catalítica por enzima.
- `${TAG}_kcat.csv`: Predições de kcat e otimização evolutiva de variantes.

**Tempo de execução:** Pode levar horas/dias dependendo do número de sequências (~3177 no exemplo) e iterações evolutivas. Monitore com `nvidia-smi` para confirmar uso da GPU.

## Resultados
- Os CSVs contêm scores de feasibility, kcat preditos e sequências otimizadas.
