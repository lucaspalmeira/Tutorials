# Dinâmica molecular clássica (proteína-ligante) utilizando o AMBER

> Este tutorial assume que o complexo proteína-ligante foi preparado através do servidor CHARMM-GUI em `Solution Builder`.

## Geração de `dihe.restraint` para glicanos (inulina) no AMBER

Esta etapa descreve **passo a passo** como identificar um ligante glicano preparado pelo **CHARMM-GUI (Solution Builder)** e gerar manualmente um arquivo **`dihe.restraint`** para uso em simulações de dinâmica molecular no **AMBER**.

O procedimento é especialmente útil quando o CHARMM-GUI **não cria automaticamente** restrições diedrais para carboidratos/polímeros (caso comum para frutanos/inulina). Caso o arquivo tenha sido criado, ignore esta etapa e pule para a etapa de minimização.

Assim, o workflow descrito aqui foi aplicado a um oligômero de frutose β-(2→1) contendo os resíduos 571–576.

---

### Contexto

* Sistema preparado no **CHARMM-GUI**
* Dinâmica será executada no **AMBER**
* Ligante: **inulina (β-frutano)**
* Resíduos do ligante: `1CU` e `0CU`
* Arquivos principais:

  * `amber/step3_input.parm7`
  * `amber/step3_input.rst7`
  * `residues.dat`
  * `calc_inulin_dihedrals.in`
  * `write_dihe_restraint.py`

---

### Identificar os resíduos do ligante

Dentro do diretório do sistema, gere a lista de resíduos:

```bash
cpptraj step3_input.parm7 << EOF > residues.dat
resinfo :*
EOF
```

A partir do arquivo `residues.dat`, o ligante foi identificado como:

```
  571 1CU    8794   8814     21   571     2    
  572 1CU    8815   8835     21   572     2    
  573 1CU    8836   8856     21   573     2    
  574 1CU    8857   8877     21   574     2    
  575 1CU    8878   8898     21   575     2    
  576 0CU    8899   8920     22   576     2 
```

Isso indica uma cadeia de **6 unidades de frutose**, sendo a última terminal (`0CU`).

---

### Inspecionar nomes e índices dos átomos

Entre no `cpptraj`:

```bash
cpptraj step3_input.parm7
```

Liste os átomos de um resíduo do ligante:

```cpptraj
atominfo :571
```

Átomos relevantes para diedros glicosídicos:

* `O5`
* `C2`
* `O1`
* `C1`

Repita para o próximo resíduo:

```cpptraj
atominfo :572
```

---

### Definição correta dos diedros para β(2→1)-frutano

Para cada ligação glicosídica entre os resíduos *i* e *i+1*:

#### Diedro φ (phi)

```
O5(i) – C2(i) – O1(i+1) – C1(i+1)
```

#### Diedro ψ (psi)

```
C2(i) – O1(i+1) – C1(i+1) – C2(i+1)
```

Esses são os **únicos diedros que devem ser restringidos**.

**Nota:** Nunca restrinja diedros internos do anel.

---

### Cálculo dos diedros com cpptraj

Criar o arquivo `calc_inulin_dihedrals.in`:

```cpptraj
parm step3_input.parm7
trajin step3_input.rst7 1 1

# Ligação 571–572
dihedral phi_571_572 :571@O5 :571@C2 :572@O1 :572@C1 out phi_571_572.dat
dihedral psi_571_572 :571@C2 :572@O1 :572@C1 :572@O5 out psi_571_572.dat

# Ligação 572–573
dihedral phi_572_573 :572@O5 :572@C2 :573@O1 :573@C1 out phi_572_573.dat
dihedral psi_572_573 :572@C2 :573@O1 :573@C1 :573@O5 out psi_572_573.dat

# Ligação 573–574
dihedral phi_573_574 :573@O5 :573@C2 :574@O1 :574@C1 out phi_573_574.dat
dihedral psi_573_574 :573@C2 :574@O1 :574@C1 :574@O5 out psi_573_574.dat

# Ligação 574–575
dihedral phi_574_575 :574@O5 :574@C2 :575@O1 :575@C1 out phi_574_575.dat
dihedral psi_574_575 :574@C2 :575@O1 :575@C1 :575@O5 out psi_574_575.dat

# Ligação 575–576
dihedral phi_575_576 :575@O5 :575@C2 :576@O1 :576@C1 out phi_575_576.dat
dihedral psi_575_576 :575@C2 :576@O1 :576@C1 :576@O5 out psi_575_576.dat

run
```

Execute:

```bash
cpptraj -i calc_inulin_dihedrals.in
```
* Serão gerados arquivos .dat com os valores dos diedros na estrutura inicial.

Verifique os valores:

```bash
head phi_571_572.dat
head psi_571_572.dat
```

---

### Criação automática do arquivo `dihe.restraint` (Python)

Criar o script `write_dihe_restraint.py`:

```python
import glob

RK = 20.0        # força da restrição (kcal/mol·rad²)
DELTA = 10.0     # largura do poço (± graus)
OUTFILE = "dihe.restraint"

# Diedros φ — β(2→1)-frutano obit
# Mapeamento diedro → índices de átomos (iat)
# (obtidos via atominfo)
DIHEDRALS = {
    "phi_571_572": [8795, 8794, 8835, 8832],
    "phi_572_573": [8816, 8815, 8856, 8853],
    "phi_573_574": [8837, 8836, 8877, 8874],
    "phi_574_575": [8858, 8857, 8898, 8895],
    "phi_575_576": [8879, 8878, 8920, 8916],
}

def read_dihedral_value(filename):

    with open(filename) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                return float(line.split()[1])
    raise RuntimeError(f"Valor não encontrado em {filename}")

def build_restraint_window(angle, delta):
    """
    Constrói (R1, R2, R3, R4) garantindo:
    R1 <= R2 < R3 <= R4  e  limites dentro de [-180, 180]
    """
    r1 = -180.0
    r4 =  180.0

    r2 = angle - delta
    r3 = angle + delta

    if r2 < -180.0:
        r2 = -180.0
    if r3 > 180.0:
        r3 = 180.0

    # caso delta seja grande demais
    if r2 >= r3:
        raise ValueError(
            f"Janela inválida para diedro {angle:.2f}° "
            f"(r2={r2:.2f}, r3={r3:.2f})"
        )

    return r1, r2, r3, r4

with open(OUTFILE, "w") as out:
    for dat in sorted(glob.glob("phi_*.dat")):
        key = dat.replace(".dat", "")
        if key not in DIHEDRALS:
            continue

        angle = read_dihedral_value(dat)
        r1, r2, r3, r4 = build_restraint_window(angle, DELTA)
        iat = DIHEDRALS[key]

        out.write("&rst\n")
        out.write(f" iat={iat[0]},{iat[1]},{iat[2]},{iat[3]},\n")
        out.write(
            f" r1={r1:.1f}, r2={r2:.1f}, r3={r3:.1f}, r4={r4:.1f},\n"
        )
        out.write(f" rk2={RK:.1f}, rk3={RK:.1f},\n")
        out.write("/\n\n")

print("Arquivo dihe.restraint gerado.")
```

Execute:

```bash
python write_dihe_restraint.py
```

A saída será:

```
dihe.restraint
```

---

### Passar para a flag DISANG (em .mdin) o arquivo `dihe.restraint`

Nos arquivos `step4.0_minimization.mdin` e `step4.1_equilibration.mdin`:

```ini
&cntrl
  nmropt=1,
/
&wt type='END' /
DISANG=dihe.restraint
```

---

### Boas práticas recomendadas

* Use `dihe.restraint` **somente até o fim da equilibration**
* Produção → **remova completamente**
* `rk2 = rk3 = 10–20` é ideal
* Nunca restrinja diedros do anel
* Método compatível com literatura de MD de glicanos

---

Este procedimento garante:

* Estabilidade conformacional inicial do glicano
* Evita colapsos não físicos

---

**Execução da dinâmica molecular**

Para especificar a GPU device=1 (ou seja, a segunda GPU, já que a contagem começa em 0) ao rodar pmemd.cuda no AMBER, a forma padrão e recomendada é usar a variável de ambiente CUDA_VISIBLE_DEVICES. Atualmente, a seleção da GPU a ser usada em execuções com uma única GPU é automática se as GPUs estiverem configuradas no modo exclusivo de processo (nvidia-smi -c 3), mas a abordagem recomendada é usar a variável de ambiente CUDA_VISIBLE_DEVICES para selecionar qual GPU deve ser usada.

Desta forma, execute:

```bash
export CUDA_VISIBLE_DEVICES=1
```

ou

```bash
CUDA_VISIBLE_DEVICES="1" pmemd.cuda ...
```

---

### 1. Minimização

```bash
sed -e "s/FC/1.0/g" dihe.restraint > step4.0_minimization.rest
```

```bash
pmemd.cuda -O -i step4.0_minimization.mdin -p step3_input.parm7 -c step3_input.rst7 -o step4.0_minimization.mdout -r step4.0_minimization.rst7 -inf step4.0_minimization.mdinfo -ref step3_input.rst7
```

### 2. Equilibração

```bash
sed -e "s/FC/1.0/g" dihe.restraint > step4.1_equilibration.rest
```

```bash
pmemd.cuda -O -i step4.1_equilibration.mdin -p step3_input.parm7 -c step4.0_minimization.rst7 -o step4.1_equilibration.mdout -r step4.1_equilibration.rst7 -inf step4.1_equilibration.mdinfo -ref step3_input.rst7 -x step4.1_equilibration.nc
```

### Produção

```bash
pmemd.cuda -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5_production.mdout -r step5_production.rst7 -inf step5_production.mdinfo -x step5_production.nc
```

---

## Executando Simulações Aceleradas por GPU
Para executar uma simulação de dinâmica molecular acelerada por GPU, a única alteração necessária é usar o executável
pmemd.cuda em vez de pmemd. Exemplo:

```bash
pmemd.cuda -O -i mdin -o mdout -p prmtop -c inpcrd -r restrt -x mdcrd
```
Isso executará automaticamente o cálculo na GPU com mais memória, mesmo que essa GPU já esteja em uso.

Se você tiver apenas uma GPU compatível com CUDA em sua máquina, isso é suficiente; no entanto, se você quiser controlar
qual GPU será usada ou se quiser executar várias simulações independentes usando GPUs diferentes, você precisará especificar manualmente a GPU a ser usada com a variável de ambiente CUDA VISIBLE DEVICES.

CUDA_VISIBLE_DEVICES Especifica qual GPU deve ser usada para executar um cálculo PMEMD acelerado por GPU. Isso se baseia no ID de hardware da placa de vídeo, que pode ser obtido desativando a variável (unset CUDA_VISIBLE_DEVICES) e executando o comando deviceQuery do SDK CUDA da NVIDIA. Os valores válidos são uma lista de números inteiros de 0 a 32. Várias GPUs podem ser
listadas com vírgulas entre elas, e aquela com mais memória será selecionada. Por
exemplo:

```bash
export CUDA VISIBLE DEVICES=1.3
pmemd.cuda -O -i mdin -o mdout -p prmtop -c inpcrd -r restrt -x mdcrd
```

---

## Pós-processamento de MD no AMBER (Produção)

Arquivos principais:

* Topologia: `step3_input.parm7`
* Trajetória de produção: `step5_production.nc`
* Coordenadas de referência: `step5_production.rst7`

O ligante corresponde aos resíduos **1CU / 0CU** conforme o `step3_input.pdb`.

---

### 1. Centralização da trajetória (cpptraj)

Centraliza o sistema na proteína, removendo PBC e alinhando ao primeiro frame.

> Substitua `XXX` pelo último resíduo da proteína (excluindo solvente e ligante) em todas as etapas seguintes.


```bash
cpptraj -p step3_input.parm7 -y step5_production.nc << EOF
center :1-XXX mass origin
image origin center
rms first :1-XXX@CA
trajout step5_centered.nc
EOF
```

---

### 2. RMSD do ligante

Cálculo do RMSD do ligante ao longo do tempo (ajustando pela proteína).

```bash
cpptraj -p step3_input.parm7 -y step5_centered.nc << EOF
rms Protein first :1-XXX@CA
rms Ligand first :1CU,0CU out rmsd_ligand.dat
EOF
```

Arquivos gerados:

* `rmsd_ligand.dat`

---

### 3. RMSF por resíduo (proteína)

Flutuação por resíduo baseada nos átomos Cα.

```bash
cpptraj -p step3_input.parm7 -y step5_centered.nc << EOF
rms first :1-XXX@CA
atomicfluct out rmsf_ca.dat :1-XXX@CA byres
EOF
```

Arquivos gerados:

* `rmsf_ca.dat`

---

### 4. Raio de giro (Rg)

Cálculo do raio de giro da proteína.

```bash
cpptraj -p step3_input.parm7 -y step5_centered.nc << EOF
radgyr :1-XXX out rg_protein.dat
EOF
```

Arquivos gerados:

* `rg_protein.dat`

---

### 5. Ligações de hidrogênio ligante–proteína

Detecção de H-bonds entre ligante e proteína.

```bash
cpptraj -p step3_input.parm7 -y step5_centered.nc << EOF
hbond HB out hbond_lig_prot.dat \
  solventdonor :WAT \
  donormask :1CU,0CU \
  acceptormask :1-XXX
EOF
```

Arquivos gerados:

* `hbond_lig_prot.dat`

---

### 6. MM-PBSA e MM-GBSA

#### 6.1 Arquivo de entrada (`mmpbsa.in`)

```
&general
  startframe=1,
  endframe=LAST,
  interval=10,
  verbose=1,
/
&gb
  igb=5,
/
&pb
  istrng=0.150,
/
```

#### 6.2 Execução

```bash
MMPBSA.py -O -i mmpbsa.in -cp step3_input.parm7 -rp receptor.parm7 -lp ligand.parm7 -y step5_centered.nc
```

**Nota:** `receptor.parm7` e `ligand.parm7` devem ser previamente gerados com `ante-MMPBSA.py` ou `cpptraj`.

Arquivos gerados:

* `FINAL_RESULTS_MMPBSA.dat`

---

### 7. Clustering (conformação mais representativa)

Clusterização baseada no RMSD do ligante.

```bash
cpptraj -p step3_input.parm7 -y step5_centered.nc << EOF
rms first :1CU,0CU
cluster hieragglo clusters 3 linkage average \
  summary cluster_summary.dat \
  repout cluster_rep repfmt pdb
EOF
```

Arquivos gerados:

* `cluster.dat` – estatísticas dos clusters
* `cluster_rep.c0.pdb` – estrutura representativa do cluster dominante

---

> As instruções acima foram extraídas e interpretadas do manual do Amber 2025 (https://ambermd.org/doc12/Amber25.pdf)
