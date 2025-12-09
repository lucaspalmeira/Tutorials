# QM/MM no AMBER: do CHARMM-GUI à Produção com QUICK (B3LYP) e Métodos Semiempíricos

Este README descreve **como executar minimização, equilíbrio e produção** utilizando os arquivos gerados pelo **CHARMM-GUI** para **AMBER**, incluindo:

* Explicação do workflow
* Comandos completos, sem variáveis
* Execução com CPU (sander) ou GPU (pmemd.cuda)
* Produção realizada com o método B2LYP (DFT) ou semiempírico

Os nomes de arquivos assumidos são:

* `step3_input.parm7`
* `step3_input.rst7`
* `step4.0_minimization.mdin`
* `step4.1_equilibration.mdin`
* `step5_production.mdin`

Se seus arquivos possuírem outros nomes, substitua conforme necessário.

---

## 1. MINIMIZAÇÃO

Este comando executa a minimização utilizando **sander** (CPU). Caso exista o arquivo `dihe.restraint`, é necessário atualizar o valor de força antes.

```bash
grep -q "FC" dihe.restraint && sed -e "s/FC/1.0/g" dihe.restraint > step4.0_minimization.rest
```

### **Executar a minimização:**

```bash
sander -O -i step4.0_minimization.mdin -p step3_input.parm7 -c step3_input.rst7 -o step4.0_minimization.mdout -r step4.0_minimization.rst7 -inf step4.0_minimization.mdinfo -ref step3_input.rst7
```

## Se quiser executar em GPU:

```bash
pmemd.cuda -O ... (mesmos argumentos)
```

---

## 2. EQUILIBRAÇÃO

Se `dihe.restraint` existir:

```bash
sed -e "s/FC/1.0/g" dihe.restraint > step4.1_equilibration.rest
```

### **Executar a equilibração:**

```bash
sander -O -i step4.1_equilibration.mdin -p step3_input.parm7 -c step4.0_minimization.rst7 -o step4.1_equilibration.mdout -r step4.1_equilibration.rst7 -inf step4.1_equilibration.mdinfo -ref step3_input.rst7 -x step4.1_equilibration.nc
```

### GPU:

```bash
pmemd.cuda -O ...
```

---

## 3. PRODUÇÃO

```bash
sander -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5.mdout -r step5.rst7 -inf step5.mdinfo -x step5.nc
```

### Para executar em GPU (recomendado):

```bash
mpirun -np 2 sander.MPI -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5.mdout -r step5.rst7 -inf step5.mdinfo -x step5.nc
```

**Nota:** Para realizar uma dinâmica QM/MM utilizando a teoria semiempírica AM1, não é possível executar com aceleração por GPU através do `pmemd.cuda`. Para executar uma simulação QM/MM na qual a região quântica será calculada em GPU, deverá ser feito uso do QUICK. O QUICK é um software de cálculo quântico que pode ser utilizado diretamente através do AmberTools.

---

# Dinâmica moléculas QM/MM com o QUICK

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

O arquivo de entrada (.mdin) deverá conter:

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

Exemplo para simulação QM/MM com **embedding eletrostático** (o mais comum). O mesmo arquivo de entrada (.mdin) funciona para todas as versões do QUICK (serial, MPI, single-GPU ou multi-GPU):

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

### Variáveis do namelist &quick (resumo completo)

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

E no input (.mdin):

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

---

## Neste estudo, utilizamos as seguintes flags no arquivo `.mdin` para a simulação de produção executada com o programa **QUICK**:

```input
 &qmmm
  iqmatoms=593, 594, 595, 596, 597, 598, 3722, 3723, 3724, 3725, 3726, 3727, 3728, 3729, 3730, 3731, 45573, 45574, 45575, 45624, 45625, 45626, 56283, 56284, 56285, 58281, 58282, 58283, 58680, 58681, 58682, 58608, 58609, 58610, 58863, 58864, 58865, 9066, 9067, 9068, 9069, 9070, 9071, 9072, 9073, 9074, 9075, 9076, 9077, 9078, 9079, 9080, 9081, 9082, 9083, 9084, 9085, 9086, 9087, 9088, 
  qmcharge=-1,
  qm_theory='PM3',
  qmcut=12.0,
  qmshake=0,
  adjust_q=1
 /
 &quick
  method = 'B3LYP',
  basis = 'def2-svp',
 /
```

A simulação foi executada com o comando:

```bash
mpirun -np 2 sander.quick.cuda.MPI -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5.mdout -r step5.rst7 -inf step5.mdinfo -x step5.nc
```

Nota: O parâmetro -np foi limitado a 2 porque o uso de mais núcleos de CPU provoca overflow instantâneo de memória nas GPUs (2× RTX 4090).
Caso ocorra erro MPI_ABORT e o arquivo quick.out contenha a mensagem:

```output
QMMM: System specified with odd number of electrons (107)
QMMM: but odd spin (1). You most likely have the charge of
QMMM: QM region (qmcharge) set incorrectly.
```

Isso acontece porque a região QM selecionada possui 107 elétrons (número ímpar). Um sistema com número ímpar de elétrons não pode ser um singlet (spin = 1 é o valor padrão do Amber para singlet). Ele deve ser, no mínimo, um dubleto (multiplicidade 2 → 1 elétron desemparelhado).
Para corrigir, adicionamos explicitamente as flags spin=1 (dubleto) e qmcharge=-1 ao final do arquivo .mdin:

```input
 &qmmm
  iqmatoms=3729, 45573, 45574, 45575, 3728, 3730, 3731, 598, 596, 597, 45624, 45625, 45626, 56283, 56284, 56285, 58281, 58282, 58283, 58590, 58591, 58592, 9067, 9068, 9066, 9089, 9106, 9109, 9110, 9107, 9108, 
  qmcharge=-1,
  spin=1,
  qm_theory='quick',
  qmcut=12.0,
  qmshake=0,
  adjust_q=1
 /
 &quick
  method = 'B3LYP',
  basis = 'def2-svp',
 /
```

---

Infelizmente, mesmo utilizando `mpirun -np 2 sander.quick.cuda.MPI`, a estimativa de tempo foi de 6123,8 horas (~0,22 ns/dia), o que implicaria cerca de 255 dias para completar a simulação mesmo com aceleração por GPU.

Devido ao elevado custo computacional, a dinâmica molecular QM/MM foi realizada utilizando um método semiempírico com `qm_theory='AM1'`. 

As flags utilizadas no arquivo .mdin da produção foram:

```ìnput
 &qmmm
  iqmatoms=3729, 45573, 45574, 45575, 3728, 3730, 3731, 598, 596, 597, 45624, 45625, 45626, 56283, 56284, 56285, 58281, 58282, 58283, 58590, 58591, 58592, 9067, 9068, 9066, 9089, 9106, 9109, 9110, 9107, 9108, 
  qmcharge=-1,
  spin=1,
  qm_theory='AM1',
  qmcut=12.0,
  qmshake=1,
  adjust_q=1,
  qm_ewald=1, qm_pme=1,
 /
```

Mantivemos `spin=1` e `qmcharge=-1`. A simulação foi então executada utilizando apenas CPU (28 núcleos):

```bash
mpirun -np 28 sander.MPI -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5.mdout -r step5.rst7 -inf step5.mdinfo -x step5.nc
```

Embora este protocolo não utilize GPU, o tempo estimado caiu significativamente para 112,3 horas (mantendo ~0,22 ns/dia), tornando o cálculo viável.

---

## Possíveis erros:

**Erro *“QM region + cutoff larger than box”* no AMBER**

Durante uma simulação QM/MM com o AMBER (sander/pmemd), pode surgir o erro:

```
****************************************************
ERROR: QM region + cutoff larger than box dimension:
QM-MM Cutoff = 35.0000
 Coord   Lower     Upper    Size    Radius of largest sphere inside unit cell
   X   -40.743    39.788    80.531    46.761
   Y   -40.855    52.676    93.531    46.761
   Z   -42.214    43.047    85.261    46.761
****************************************************
SANDER BOMB in subroutine QM_CHECK_PERIODIC<qm_mm.f>
QM region + cutoff larger than box
cannot continue, need larger box.
```

---

### Causa do erro

O AMBER verifica se a região QM “cabe” dentro da caixa de simulação, considerando:

* O **raio da região QM** (distância máxima entre o centróide QM e seus átomos)
* O **cutoff QM/MM** (`qmcut`)
* A **menor dimensão da caixa** (por exemplo, `80.531 Å` no eixo X)

O parâmetro qmcut define o raio de corte eletrostático, em angstrons, usado para determinar quais átomos da região clássica (MM) interagem eletrostaticamente com a região quântica (QM) em simulações QM/MM no AMBER. Em outras palavras, qualquer átomo MM que esteja dentro dessa distância em relação a qualquer átomo QM será incluído na lista de pares que participam das interações eletrostáticas QM–MM. Esse valor, por padrão, é igual ao cutoff clássico (cut) usado para as interações MM–MM, e geralmente não precisa ser alterado. Em cálculos que usam Ewald ou PME, o qmcut afeta apenas a parte direta (real-space) da decomposição eletrostática, enquanto a parte recíproca continua independente do cutoff. É importante destacar que esse parâmetro não altera as interações internas da região QM, já que todos os átomos QM interagem entre si independentemente da distância, e também não influencia as interações de van der Waals entre QM e MM, que continuam sendo tratadas classicamente com o cutoff especificado por cut.

A condição necessária para **não gerar erro** é:

```
QM_radius + qmcut  <  half_box_dimension
```

Onde:

```
half_box_dimension = (menor dimensão da caixa) / 2
```

Se essa condição não for satisfeita, o AMBER aborta com o erro acima.

---

### Exemplo do problema encontrado

Com:

* `qmcut = 35 Å`
* `QM_radius = 8.147 Å`
* menor dimensão da caixa = `80.531 Å`

Temos:

```
QM_radius + qmcut = 43.147 Å
half_box_dimension = 80.531 / 2 = 40.265 Å
```

Comparando:

```
43.147 Å  >  40.265 Å   ← NÃO CABE NA CAIXA
```

Portanto, o AMBER interrompe a execução.

---

### Solução

**Reduzir o `qmcut`**

Se você usar, por exemplo, `qmcut = 30 Å`:

```
QM_radius + qmcut = 38.147 Å
half_box_dimension = 40.265 Å
```

Agora:

```
38.147 < 40.265  ← ✔️ OK
```

A simulação QM/MM roda sem erros.

---

### Como verificar automaticamente se sua caixa comporta a região QM

Use o script abaixo para calcular:

* centróide da região QM
* raio da região QM
* tamanho mínimo necessário da caixa
* se o `qmcut` escolhido é possível

---

### Script: **calc_qm_radius.py**

salve como `calc_qm_radius.py` e execute: 
```bash
python calc_qm_radius.py step3_input.pdb <qmcut> <size>
```

```python
import sys, math, statistics
fn=sys.argv[1]
# lista do seu &qmmm
iqm_list=[3729,45573,45574,45575,3728,3730,3731,598,596,597,45624,45625,45626,56283,56284,56285,58281,58282,58283,58590,58591,58592,9067,9068,9066,9089,9106,9109,9110,9107,9108]

atoms=[]
with open(fn) as f:
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            serial=int(line[6:11].strip())
            x=float(line[30:38].strip()); y=float(line[38:46].strip()); z=float(line[46:54].strip())
            resseq=line[22:26].strip()
            atoms.append((serial,int(resseq) if resseq.isdigit() else None,x,y,z))

# try match by atom serials first, else by residue numbers
coords=[(x,y,z) for (s,rs,x,y,z) in atoms if s in iqm_list]
if not coords:
    coords=[(x,y,z) for (s,rs,x,y,z) in atoms if rs in iqm_list]

if not coords:
    print("Nenhum átomo IQM encontrado: verifique se os números em iqm_list são ATOM serials ou residue IDs.")
    print("Total de átomos no PDB:", len(atoms))
    sys.exit(1)

cx=sum(x for x,y,z in coords)/len(coords)
cy=sum(y for x,y,z in coords)/len(coords)
cz=sum(z for x,y,z in coords)/len(coords)
dists=[math.sqrt((x-cx)**2+(y-cy)**2+(z-cz)**2) for x,y,z in coords]
qm_radius=max(dists)
print("IQM atom count:", len(coords))
print("QM centroid: {:.3f} {:.3f} {:.3f}".format(cx,cy,cz))
print("QM radius (max distance to centroid): {:.3f}Å".format(qm_radius))
qmcut = int(sys.argv[2])
required_half_box = qm_radius + qmcut
print(f'required half box: {required_half_box}Å')
required_box_size = 2 * required_half_box
print(f'required_box_size: {required_box_size}Å')

menor_dimensao_box = float(sys.argv[3]) #80.531
metade_da_menor_dimensao_box = menor_dimensao_box / 2

if required_half_box > metade_da_menor_dimensao_box:
    print(f"ATENÇÃO: QM region + QMcut ({required_half_box:.2f} Å) > metade da caixa ({metade_da_menor_dimensao_box:.2f} Å)")
    print(f"         Solução: aumente o box ou reduza QMcut (atual = {qmcut} Å).")
else:
    print(f"OK → QM region + QMcut = {required_half_box:.2f} Å (cabe na caixa de {menor_dimensao_box:.2f} Å)")
```

---

# Como usar

```bash
python calc_qm_radius.py step3_input.pdb 35 80.531
```

Exemplo de saída indicando problema:

```
IQM atom count: 31
QM centroid: -8.948 8.242 8.976
QM radius (max distance to centroid): 8.147Å
required half box: 38.14681655042923Å
required_box_size: 76.29363310085846Å
 QM_radius (qm_radius) + QMcut (30) = 38.14681655042923 < 40.2655(Metade da menor dimensao caixa)
```

```bash
python calc_qm_radius.py step3_input.pdb 30 80.531
```

Exemplo indicando que está OK:

```
IQM atom count: 31
QM centroid: -8.948 8.242 8.976
QM radius (max distance to centroid): 8.147Å
required half box: 38.14681655042923Å
required_box_size: 76.29363310085846Å
OK → QM region + QMcut = 38.15 Å (cabe na caixa de 80.53 Å)
```

> As instruções acima foram extraídas e interpretadas do manual do Amber 2025 (https://ambermd.org/doc12/Amber25.pdf)
