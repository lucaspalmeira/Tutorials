# Execução de Dinâmica Molecular no AMBER (CHARMM-GUI)

Este README descreve **como executar minimização, equilíbrio e produção** utilizando os arquivos gerados pelo **CHARMM-GUI** para **AMBER**, incluindo:

* Explicação do workflow
* Comandos completos, sem variáveis
* Execução com CPU (sander) ou GPU (pmemd.cuda)
* Produção em **um único comando**, sem divisão em etapas

Os nomes de arquivos assumidos são:

* `step3_input.parm7`
* `step3_input.rst7`
* `step4.0_minimization.mdin`
* `step4.1_equilibration.mdin`
* `step5_production.mdin`

Se seus arquivos possuírem outros nomes, substitua conforme necessário.

---

# 1. MINIMIZAÇÃO

Este comando executa a minimização utilizando **sander** (CPU). Caso exista o arquivo `dihe.restraint`, é necessário atualizar o valor de força antes.

```bash
grep -q "FC" dihe.restraint && sed -e "s/FC/1.0/g" dihe.restraint > step4.0_minimization.rest
```

### **Rodar a minimização:**

```bash
sander -O -i step4.0_minimization.mdin -p step3_input.parm7 -c step3_input.rst7 -o step4.0_minimization.mdout -r step4.0_minimization.rst7 -inf step4.0_minimization.mdinfo -ref step3_input.rst7
```

### Se quiser rodar em GPU:

```bash
pmemd.cuda -O ... (mesmos argumentos)
```

---

# 2. EQUILIBRAÇÃO

Se `dihe.restraint` existir:

```bash
sed -e "s/FC/1.0/g" dihe.restraint > step4.1_equilibration.rest
```

### **Rodar a equilibração:**

```bash
sander -O -i step4.1_equilibration.mdin -p step3_input.parm7 -c step4.0_minimization.rst7 -o step4.1_equilibration.mdout -r step4.1_equilibration.rst7 -inf step4.1_equilibration.mdinfo -ref step3_input.rst7 -x step4.1_equilibration.nc
```

### GPU:

```bash
pmemd.cuda -O ...
```

---

# 3. PRODUÇÃO

```bash
sander -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5.mdout -r step5.rst7 -inf step5.mdinfo -x step5.nc
```

### Para executar em GPU (recomendado):

```bash
mpirun -np 18 sander.MPI -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5.mdout -r step5.rst7 -inf step5.mdinfo -x step5.nc
```

---

**Nota:** Para a um dinâmica QM/MM utilizando a teoria semiempírica AM1, não é possível executar utilizando GPU através de pmemd.cuda. Para executar uma QM/MM onde a região quântica será calculada com GPU deverá ser feito uso do QUICK. O QUICK é um software para cálculo quântico que pode ser utilizando através do AmberTools.

## Dinâmica moléculas QM/MM com o QUICK

O **SANDER** pode acessar diferentes tipos de instalação do QUICK para simulações QM/MM:  
- serial,  
- paralelo com MPI,  
- acelerado por CUDA,  
- acelerado por HIP,  
- paralelo CUDA + MPI,  
- paralelo HIP + MPI.

Se estiver usando a **interface baseada em arquivos (FBI - file-based interface)**, o executável **sander** é capaz de chamar qualquer um dos diferentes executáveis do QUICK:  
`quick`, `quick.MPI`, `quick.cuda`, `quick.hip`, `quick.cuda.MPI` ou `quick.hip.MPI`.

Se estiver usando a **interface de programação de aplicações (API)**, é necessário usar um executável diferente do SANDER para cada tipo de QUICK:  
- As versões serial e paralela com MPI do QUICK são acessadas a partir dos executáveis `sander` e `sander.MPI`, respectivamente;  
- As versões aceleradas por single-GPU e multi-GPU do QUICK são acessadas a partir dos executáveis  
  `sander.quick.cuda` / `sander.quick.hip` (single-GPU)  
  e  
  `sander.quick.cuda.MPI` / `sander.quick.hip.MPI` (multi-GPU).  

Esses executáveis são idênticos ao `sander` e `sander.MPI` em todas as funcionalidades do SANDER, exceto que realizam os cálculos QM/MM usando o código acelerado por GPU do QUICK por meio da API.

Exemplos de como usar tanto a API quanto a funcionalidade FBI podem ser encontrados nos conjuntos de testes nos seguintes locais dentro do código-fonte do AMBER:  
- `$AMBERHOME/test/qmmm_QUICK` → exemplos com API  
- `$AMBERHOME/test/qmmm_EXTERN/*Quick` → exemplos com FBI

**Nota importante:** Antes de executar qualquer simulação QM/MM com QUICK, os usuários devem obrigatoriamente executar o comando  
`source $AMBERHOME/amber.sh`  
(ou `source $AMBERHOME/amber.csh`, dependendo do seu shell).  
Esse passo garante que a localização dos executáveis e bibliotecas necessários seja corretamente configurada nas variáveis de ambiente do sistema operacional.

## Interface baseada em arquivos (FBI)

Abaixo está um exemplo das modificações necessárias no arquivo de entrada do SANDER para realizar uma simulação QM/MM com **embedding mecânico**. Neste exemplo, os dois primeiros resíduos do sistema são colocados na região QM e a simulação é executada no nível **B3LYP/def2-SVP** usando o `quick.cuda.MPI` com 2 GPUs:

```bash
&cntrl
  ...
  ifqnt = 1,
/
&qmmm
  qmmask = ':1-2',
  qm_theory = 'extern',     ! indica que será usada a interface FBI
  qmmm_int = 5,             ! 5 = embedding mecânico
/
&quick
  method = 'B3LYP',
  basis = 'def2-svp',
  executable = 'quick.cuda.MPI',      ! executável QUICK desejado
  do_parallel = 'mpirun -np 2',       ! só necessário quando usar versão MPI do QUICK
/
```

**Explicação das flags principais:**
- `ifqnt = 1` → ativa QM/MM
- `qmmask` → define quais átomos entram na região quântica
- `qm_theory = 'extern'` → diz ao SANDER que vai usar um programa externo (FBI)
- `qmmm_int = 5` → embedding mecânico (cargas MM não entram no Hamiltoniano QM)
- `executable` → pode ser qualquer um: `quick`, `quick.MPI`, `quick.cuda`, `quick.hip`, `quick.cuda.MPI`, `quick.hip.MPI`
- `do_parallel` → só é usado quando o executável é uma versão MPI

**Atenção importante:**
- Em geral, **deve-se usar o sander serial** (não o `sander.MPI`) porque há limitações em chamadas de sistema dentro de programas MPI.
- Isso significa que **a parte MM roda em apenas um núcleo** (serial), mesmo usando várias GPUs na parte QM.

## Interface de programação (API)

Exemplo para simulação QM/MM com **embedding eletrostático** (o mais comum). O mesmo arquivo de entrada funciona para todas as versões do QUICK (serial, MPI, single-GPU ou multi-GPU):

```bash
&cntrl
  ...
  ifqnt = 1,
/
&qmmm
  qmmask = ':1-2',
  qm_theory = 'quick',      ! indica que será usada a API do QUICK
  qmmm_int = 1,             ! 1 = embedding eletrostático
  qm_ewald = 0,             ! atualmente é obrigatório (cutoff normal)
/
&quick
  method = 'B3LYP',
  basis = 'def2-svp',
/
```

**Pontos importantes:**
- `qm_theory = 'quick'` → ativa a API (comunicação direta em memória, muito mais rápida)
- `qmmm_int = 1` → embedding eletrostático (cargas MM entram no Hamiltoniano QM)
- O mesmo input funciona com `sander`, `sander.MPI`, `sander.quick.cuda` ou `sander.quick.cuda.MPI`
- Será gerado um arquivo `quick.out` com toda a saída detalhada do QUICK a cada passo de MD

### 11.3.1.3. Variáveis do namelist &quick (resumo completo)

| Variável          | Tipo     | Apenas API/FBI | Descrição e valores comuns                                                                                  |
|-------------------|----------|----------------|-------------------------------------------------------------------------------------------------------------|
| `method`          | String   | ambas          | Método: HF, B3LYP, ωB97X-D, PBE, etc. (padrão: BLYP)                                                        |
| `basis`           | String   | ambas          | Conjunto de base: 6-31G*, def2-SVP, cc-pVDZ, etc. (padrão: 6-31G)                                           |
| `executable`      | String   | **FBI apenas** | Qual executável QUICK usar: `quick.cuda.MPI`, `quick.hip`, etc. (padrão: quick)                             |
| `do_parallel`     | String   | **FBI apenas** | Comando antes do executável: `'mpirun -np 4'` ou `'srun --gpus=4'` (sem padrão)                             |
| `scf_cyc`         | Integer  | ambas          | Máximo de ciclos SCF (padrão: 200)                                                                          |
| `reuse_dmx`       | Integer  | **API apenas** | Reusar matriz densidade do passo anterior: 0 = não, **1 = sim (padrão)** – acelera muito MD                 |
| `denserms`        | Float    | **API apenas** | Critério de convergência da matriz densidade (padrão: 1.0E6)                                                |
| `intcutoff`       | Float    | **API apenas** | Cutoff de integrais (padrão: 1.0E-8)                                                                        |
| `xccutoff`        | Float    | **API apenas** | Threshold para pruning na XC (padrão: 1.0E-8)                                                               |
| `basiscutoff`     | Float    | **API apenas** | Cutoff para funções de base insignificantes (padrão: 1.0E-6)                                                |
| `gradcutoff`      | Float    | **API apenas** | Critério de convergência do gradiente (padrão: 1.0E-7)                                                      | 
| `export`          | String   | **API apenas** | Exportar orbitais: `'molden'` (padrão: none)                                                                |
| `keywords`        | String   | **API apenas** | Em vez de usar as flags acima, pode colocar toda a linha de keywords do QUICK de uma vez                    |
|                   |          |                | Exemplo: `'B3LYP BASIS=cc-pVDZ CHARGE=0 MULT=1 GRADIENT EXTCHARGES'`                                        |
| `outfprefix`      | String   | **API apenas** | Prefixo do arquivo de saída do QUICK (padrão: `quick` → gera `quick.out`)                                   |
| `debug`           | Integer  | ambas          | 0 = sem debug (padrão), 1 = debug normal, 2 = debug extra (FBI)                                             |
| `use_template`    | Integer  | **FBI apenas** | 0 = não usar template, 1 = usar arquivo template do QUICK                                                   |

**Resumo prático:**
- Rapidez máxima e embedding eletrostático → use **API** (`qm_theory='quick'`) + `sander.quick.cuda` ou `sander.quick.cuda.MPI`
- Precisa de algo muito específico que ainda não está na API → use **FBI** (`qm_theory='extern'`), mas aceita que a parte MM será serial

## Diferença entre **FBI** (File-Based Interface) e **API** (Application Programming Interface)


| Característica                          | FBI (qm_theory = 'extern')                              | API (qm_theory = 'quick')                                      |
|-----------------------------------------|------------------------------------------------------------------|-----------------------------------------------------------------|
| Comunicação QM ↔ MM                     | Via arquivos (lento)                                            | Direto na memória (muito rápido)                                                                                                                                                              |
| Velocidade real em MD                   | 2× a 10× **mais lento** que a API                               | **Máxima velocidade possível**                                                                                                                                                           |
| Parte MM (força clássica)               | Sempre **serial** (1 núcleo) – mesmo se usar 8 GPUs no QM       | Pode ser **paralela com MPI/OpenMP** (dezenas de núcleos)                                                                                                                                                             |
| Parte QM pode usar multi-GPU?           | Sim (quick.cuda.MPI)                                            | Sim (sander.quick.cuda.MPI)                                                                                                                                              |
| Executável do SANDER que você usa       | Só o `sander` normal (serial)                                   | `sander.quick.cuda` ou `sander.quick.cuda.MPI`                                                                                                                                              |
| Embedding eletrostático                 | Sim                                                             | Sim (recomendado)                                                                                                                                                        |
| Embedding mecânico                      | Sim                                                             | Sim                                                            |
| Pode usar qualquer método/basis do QUICK? | Sim (mais flexível para coisas experimentais)                  | Sim (quase tudo já está disponível)                                                                                                                                                          |
| Arquivo de saída do QUICK               | Arquivos separados a cada passo (muitos GB)                     | Um único `quick.out` (ou outro nome que você escolher)         |
| Facilidade de uso                       | Mais complicado (tem que acertar `do_parallel`, etc.)           | Muito simples – input quase idêntico ao sander normal                                                                                                                                                               |
| Quando você é obrigado a usar FBI       | Só em casos **muito raros** hoje (ex.: alguma opção do QUICK ainda não exposta na API) | Quase nunca mais necessário                                    |

3. Use:

```bash
# Single GPU (mais simples e suficiente na maioria dos casos)
sander.quick.cuda -O -i md.in -o md.out -p topo.prmtop -c coord.rst

# ou multi-GPU + MM paralelo (máxima velocidade)
mpirun -np 8 sander.quick.cuda.MPI -O -i md.in -o md.out -p topo.prmtop -c coord.rst
```

E no input:

```input
&qmmm
  qm_theory = 'quick',     ! ← API (obrigatório para velocidade)
  qmmask = ':LIG | :120,130,145',   ! seus resíduos + ligante
/
&quick
  method = 'PM6',          ! ou B3LYP, DFTB3, etc.
  basis = '6-31G*',        ! se for ab initio
/
```
