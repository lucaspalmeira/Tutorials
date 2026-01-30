# Tutorial de Dinâmica Molecular com GROMACS (proteína-ligante)

Este tutorial descreve, passo a passo, a execução de uma simulação de Dinâmica Molecular (MD) utilizando o **GROMACS**, com foco nos estágios clássicos de **minimização de energia**, **equilibração** e **produção**. Além disso, explica o papel de cada arquivo envolvido e o significado das principais *flags* utilizadas nos comandos.

O fluxo apresentado aqui é baseado nos seguintes comandos:

* Minimização de energia
* Equilibração
* Produção

---

## 1. Conceitos básicos

A **Dinâmica Molecular** é uma técnica computacional que utiliza mecânica newtoniana para simular um sistema de átomos, permitindo estudar sua evolução temporal, estabilidade estrutural, interações e propriedades termodinâmicas.

Uma simulação típica em GROMACS segue este fluxo:

1. Preparação do sistema (topologia, caixa, solvatação, íons) que pode ser realizada no CHARMM-GUI
2. Minimização de energia
3. Equilibração (geralmente em NVT e/ou NPT)
4. Produção (trajetória final para análise)

Este tutorial cobre as etapas 2, 3 e 4.

---

## 2. Arquivos do GROMACS: o que é cada um?

### `.gro`

Arquivo de **estrutura** do sistema. Contém:

* Coordenadas atômicas
* Tipo de átomo
* Caixa de simulação

É usado como entrada (`-c`) e saída em praticamente todas as etapas.

---

### `.top`

Arquivo de **topologia** do sistema. Define:

* Tipos de átomos
* Parâmetros de força
* Conectividade (ligações, ângulos, diedros)
* Inclusão de arquivos `.itp`

É essencial para que o GROMACS saiba como calcular as forças.

---

### `.mdp`

Arquivo de **parâmetros da simulação** (*Molecular Dynamics Parameters*). Controla:

* Tipo de simulação (minimização, MD)
* Integrador
* Temperatura e pressão
* Restrições
* Tempo de simulação

Cada etapa (minimização, equilibração, produção) possui seu próprio `.mdp`.

---

### `.ndx`

Arquivo de **índices**. Contém grupos de átomos personalizados, por exemplo:

* Proteína
* Ligante
* Solvente
* Íons

Esses grupos são usados para:

* Restrições
* Análises
* Acoplamento térmico

---

### `.tpr`

Arquivo **binário portátil** (*run input file*). É gerado pelo `gmx grompp` e contém:

* Estrutura (`.gro`)
* Topologia (`.top`)
* Parâmetros (`.mdp`)
* Índices (`.ndx`)

É o arquivo que o `gmx mdrun` realmente executa.

---

### `.edr`

Arquivo de **energia**. Armazena dados como:

* Energia potencial
* Energia cinética
* Temperatura
* Pressão

Usado em análises com `gmx energy`.

---

### `.xtc`

Arquivo de **trajetória comprimida**. Contém:

* Coordenadas ao longo do tempo
* Menor precisão (mas menor tamanho)

Muito usado para análises estruturais.

---

### `.trr`

Arquivo de **trajetória completa**. Pode conter:

* Coordenadas
* Velocidades
* Forças

---

### `.cpt`

Arquivo de **checkpoint**. Permite:

* Reiniciar simulações interrompidas
* Continuação da produção

---

### `.psf`

Arquivo de **Protein Structure File**. Mais comum em CHARMM/NAMD, mas pode aparecer em fluxos híbridos. Contém a conectividade do sistema.

---

## 3. Etapa 1 – Minimização de Energia

### Objetivo

Remover contatos ruins, tensões geométricas e colisões atômicas antes da MD propriamente dita.

### Comando

```
gmx grompp -f step4.0_minimization.mdp -o step4.0_minimization.tpr -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx
```

#### Flags explicadas

* `-f` : arquivo `.mdp` da minimização
* `-o` : arquivo `.tpr` de saída
* `-c` : estrutura inicial
* `-r` : estrutura de referência (usada para restrições)
* `-p` : topologia
* `-n` : arquivo de índices

Execução:

```
gmx mdrun -v -deffnm step4.0_minimization
```

* `-v` : modo verboso
* `-deffnm` : nome base para todos os arquivos de saída
* Restrições posicionais são comuns, especialmente para o ligante.

---

## 4. Etapa 2 – Equilibração

### Objetivo

Adaptar o sistema às condições termodinâmicas desejadas (temperatura e/ou pressão), mantendo o sistema estável.

Normalmente ocorre em duas fases:

* **NVT** (volume constante)
* **NPT** (pressão constante)

Neste exemplo, ambas estão condensadas em uma única etapa.

### Comando

```
gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_input.gro -p topol.top -n index.ndx
```

Execução:

```
gmx mdrun -v -deffnm step4.1_equilibration
```

Durante a equilibração:

* Restrições posicionais são comuns, especialmente para o ligante.

---

## 5. Etapa 3 – Produção

### Objetivo

Gerar a trajetória final usada para análises científicas. Nesta etapa:

* Não há (ou há poucas) restrições

### Comando

```
gmx grompp -f step5_production.mdp -o step5_production.tpr -c step4.1_equilibration.gro -p topol.top -n index.ndx
```

Execução:

```
gmx mdrun -v -ntomp 12 -deffnm step5_production -nb gpu -gpu_id 0
```

Se o cálculo cair, execute para retomar:

```
gmx mdrun -v -ntomp 12 -deffnm step5_production -cpi -nb gpu -gpu_id 0
```

#### Flags explicadas

* `-v` : modo verboso
* `-ntomp 12` : número de threads OpenMP
* `-deffnm` : nome base para todos os arquivos de saída
* `-cpi` : continua a partir de um checkpoint
* `-nb gpu` : cálculo de interações não ligadas na GPU
* `-gpu_id 0` : seleciona a GPU 0

---

## 6. Arquivos finais gerados na produção

* `.xtc` : trajetória para análise
* `.edr` : dados energéticos
* `.log` : log da simulação
* `.cpt` : checkpoint
* `.gro` : estrutura final

---

## 7. Pré-processamento da Trajetória

Antes de qualquer análise pós-Dinâmica Molecular ou cálculo de energia livre, é essencial realizar o **pré-processamento da trajetória**. Isso garante que artefatos de condições periódicas de contorno (PBC) não afetem os resultados.

Os principais arquivos utilizados nesta etapa são:

* `traj_comp.xtc`
* `step5_production.tpr`
* `index.ndx`

---

### 7.1 Remoção de PBC e reconstrução da molécula inteira

Remove quebras artificiais causadas pela caixa periódica:

```
gmx trjconv -s step5_production.tpr -f traj_comp.xtc -o traj_noPBC.xtc -pbc mol
```

Selecionar o grupo 0 (System)

---

### 7.2 Centralização do sistema

Centraliza a proteína (ou complexo) na caixa de simulação:

```
gmx trjconv -s step5_production.tpr -f traj_noPBC.xtc -o traj_center.xtc -center
```

Selecionar 1 (Protein) e depois 0 (System)

---

### 7.3 Ajuste rotacional/translacional (fitting)

Remove movimentos globais para análises estruturais:

```
gmx trjconv -s step5_production.tpr -f traj_center.xtc -o traj_fit.xtc -fit rot+trans
```

Selecionar 4 (Backbone) e depois 0 (System)

---

Nesta seção, são apresentados os comandos mais comuns usando **gmx** e considerando os arquivos disponíveis no diretório.

Os principais arquivos usados aqui são:

* `traj_fit.xtc` : trajetória da produção
* `step5_production.tpr` : entrada binária da produção
* `index.ndx` : grupos de átomos

---

## 8. Análises Pós-Dinâmica Molecular

As análises abaixo devem ser realizadas preferencialmente usando a trajetória pré-processada (`traj_fit.xtc`).

### 8.1 RMSD (Root Mean Square Deviation)

Avalia a estabilidade estrutural ao longo do tempo em relação a uma estrutura de referência.

```
gmx rms -s step5_production.tpr -f traj_fit.xtc -n index.ndx -o rmsd.xvg
```

Se o RMSD for da proteína: selecionar 1 (Protein) e depois 1 (Protein). Se o RMSD for do ligante: selecionar 13 (LIG) e depois 13 (LIG)

Para proteína, normalmente usa-se a estrutura inicial ou média como referência e grupos como **Protein** ou **Backbone**.

---

### 7.2 RMSF (Root Mean Square Fluctuation)

Mede a flexibilidade média de cada resíduo ao longo da simulação.

```
gmx rmsf -s step5_production.tpr -f traj_fit.xtc -n index.ndx -o rmsf_residue.xvg -res
```

Selecionar 1 (Protein)

---

### 7.3 Raio de Giro (Radius of Gyration – Rg)

Avalia o grau de compactação da estrutura ao longo do tempo.

```
gmx gyrate -s step5_production.tpr -f traj_fit.xtc -n index.ndx -o gyrate.xvg
```

Selecionar 1 (Protein)

---

### 7.4 Ligações de Hidrogênio (Hydrogen Bonds)

```
gmx hbond -s step5_production.tpr -f traj_fit.xtc -n index.ndx -num hbonds.xvg
```

Selecionar 1 (Protein) e depois 13 (LIG)

---

## 9. Cálculo de Energia Livre com GMX_MMPBSA (MM/PBSA e MM/GBSA)

O **GMX_MMPBSA** permite estimar energias livres de ligação via **MM/PBSA** ou **MM/GBSA**.

### 9.1 MM/PBSA

```
gmx_MMPBSA -O -i mmpbsa.in -cs step5_production.tpr -ct traj_fit.xtc -ci index.ndx -cp topol.top -o FINAL_RESULTS_MMPBSA.dat
```

---

### 9.2 MM/GBSA

Para MM/GBSA, utilize a seção &gb no arquivo de entrada:

```
&general
   startframe=1,
   endframe=5000,
   interval=10,
/

&gb
   igb=5,
   saltcon=0.15,
/
```

Execução:

```
gmx_MMPBSA -O -i mmpbsa_gb.in -cs step5_production.tpr -ct traj_fit.xtc -ci index.ndx -cp topol.top -o FINAL_RESULTS_MMGBSA.dat
```

---

### 9.3 Decomposição de Energia por Resíduo

Permite identificar quais resíduos mais contribuem para a energia de ligação.

No arquivo de entrada:

```
&decomp
   idecomp=1,
   print_res="all",
/
```

Execução:

```
gmx_MMPBSA -O -i mmpbsa_decomp.in -cs step5_production.tpr -ct traj_fit.xtc -ci index.ndx -cp topol.top
```

As saídas incluem contribuições eletrostáticas, de van der Waals e solvatação por resíduo.

---

```
gmx_MMPBSA -O -i mmpbsa.in -cs step5_production.tpr -ct traj_comp.xtc -ci index.ndx -cp topol.top -o FINAL_RESULTS_MMPBSA.dat
```
