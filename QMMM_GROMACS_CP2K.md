# Simulação QM/MM com GROMACS 2021.3 + CP2K

Este protocolo descreve como configurar e executar uma simulação QM/MM usando GROMACS (2021.3) acoplado ao CP2K, com a região quântica (QM) composta por:

**Região QM:**
- Ligantes: resíduos 603 (glucose, AGLC) e 604 (frutose, BFRU) – sacarose
- Resíduos da proteína próximos ao sítio ativo: 64 (ASP), 82 (LEU), 121 (PHE), 122 (ASP), 193 (ARG), 194 (HIS), 271 (GLU), 272 (THR)
- Moléculas de água num raio de 4 Å em volta do ligante (serão restritas posicionalmente durante minimização e equilíbrio).

**Nota**: Para este tutorial, o complexo **proteína + ligante** deverá ser preparado no **CHARMM-GUI**, onde deverá usar o módulo `Solution Builder` para checar o $pKa$, aplicar o pH, utilizar uma caixa retangular (octaédrica tem causado o erro de sobreposição de átomos - `atom overlap`), adicionar íons (**NaCl**), utilizar o campo de força `CHARMM36` e por fim verificar a **temperatura** desejada para as etapas de equilíbrio e produção.

## 1. Criação do grupo de índice QM (QMatoms)

### 1.1 Criar grupos básicos no index.ndx

```bash
gmx make_ndx -f step3_input.gro -n index.ndx
```

Dentro do make_ndx, execute os seguintes comandos:

```bash
# Átomos selecionados do ligante que fazem parte da reação (sacarose)
a 9068 | a 9090 | a 9066 | a 9098 | a 9099 | a 9100 | a 9101 | a 9102 | a 9067
name 4 LIG

# Átomos da proteína do sítio ativo
a 596 | a 597 | a 598 | a 3728 | a 3729 | a 3730 | a 3731
name 5 atomos_sitio

q
```

### 1.2 Selecionar oxigênios de água dentro de 4 Å do ligante

```bash
gmx select -f step3_input.gro -s step4.0_minimization.tpr -select "name OH2 and (within 0.4 of group AGLC or within 0.4 of group BFRU)" -on wat_oxygens_near_lig.ndx

gmx select -f step3_input.gro -s step3_input.gro -select "resname TIP3 and same residue as (atomname OW and within 0.5 of resnr 246)" -on near246

gmx select -f structure.gro -s structure.gro -select " name OH2  and ( within 0.4 of resnr 246)" -on near246.ndx
```

O arquivo `wat_oxygens_near_lig.ndx` conterá algo como:

```
[ name_OH2_and_(within_0.4_of_group_AGLC_or_within_0.4_of_group_BFRU)_f0_t0.000 ]
44331 45630 45666 45681 56370 58383 58677 58692 58710 58767 58947 58965
```

### 1.3 Identificar os números dos resíduos de água (usando PyMOL)

Abra o arquivo `step3_input.gro` no PyMOL e execute no console:

```python
residuos = []
iterate (index 44346 or index 45525 or index 45690 or index 45726 or index 45741 or index 56379 or index 58365 or index 58665 or index 58680 or index 58698 or index 58755 or index 58926 or index 58929 or index 58944), residuos.append(resi)

residuos_unicos = sorted(list(set(residuos)))
print(residuos_unicos)
cmd.select("waters_near_lig", "resi " + "+".join(map(str, residuos_unicos)))
```

Anote os números dos resíduos de água (exemplo fictício abaixo).

### 1.4 Adicionar águas próximas ao index.ndx

```bash
gmx make_ndx -f step3_input.gro -n index.ndx
```

Dentro do make_ndx:

```bash
# Exemplo com os resíduos reais obtidos no PyMOL
r 12439 | r 12872 | r 12884 | r 12889 | r 16452 | r 17123 | r 17221 | r 17226 | r 17232 | r 17251 | r 17311 | r 17317
name 6 WAT_4A
```

# Grupo final QM: ligante + resíduos do sítio + águas próximas

```bash
4 | 5 | 6
name 7 QMatoms

q
```

## 2. Restrição posicional das águas próximas ao ligante (apenas minimização e equilíbrio)

Nos arquivos `.mdp` de **minimização** e **equilíbrio** (`step4.0_minimization.mdp` e `step4.1_equilibration.mdp`), adicione:

```mdp
; Restringir águas próximas ao ligante (evita que saiam da região QM)
pull                     = yes
pull-ngroups             = 2
pull-ncoords             = 1
pull-group1-name         = WAT_4A      ; grupo criado acima
pull-group2-name         = LIG
pull-coord1-type         = constraint
pull-coord1-geometry     = direction
pull-coord1-groups       = 1 2
pull-coord1-dim          = Y Y Y
pull-coord1-start        = yes
pull-coord1-rate         = 0
pull-coord1-k            = 1000        ; kJ mol⁻¹ nm⁻²
```

**Importante:** Essas linhas não devem contar no arquivo `step5_production.mdp` (não queremos restrição durante a produção).

## 3. Parâmetros QM/MM (adicionar em todos os .mdp: minimização, equilíbrio e produção)

No final de cada arquivo `.mdp`:

```mdp
; Parâmetros QM/MM com CP2K
qmmm-cp2k-active         = yes
qmmm-cp2k-qmgroup        = QMatoms
qmmm-cp2k-qmmethod       = PBE
qmmm-cp2k-qmcharge       = 0
qmmm-cp2k-qmmultiplicity = 1
qmmm-cp2k-core           = auto
```

## 4. Execução com Docker (GROMACS + CP2K + GPU)

```bash
# Baixar imagem
sudo docker pull kimjoochan/gromacs-cp2k:2022.2-9.1-cuda

# Iniciar container (ajuste o caminho local e a GPU conforme necessário)
sudo docker run -v "$(pwd):/home" -v "/etc/localtime:/etc/localtime:ro" -e TZ=America/Sao_Paulo -v "/tmp/.X11-unix:/tmp/.X11-unix" -e QT_X11_NO_MITSHM=1 --shm-size=1g --env="DISPLAY" --net=host --gpus '"device=1"' --privileged -itd --name gmx_cp2k kimjoochan/gromacs-cp2k:2022.2-9.1-cuda /bin/bash

# Verificar GPU dentro do container
docker exec gmx_cp2k nvidia-smi
```

**Nota:** altere '"device=1"' de acordo com a GPU disponível.

### 4.1 Minimização

Para acessar o conteiner e realizar todas as etapas seguinte da dinâmica molecular, execute:

```bash
sudo docker exec -it gmx_cp2k bash
```

```bash
cd home
```

A seguir, execute:

```bash
gmx_mpi_d grompp -f step4.0_minimization.mdp -o step4.0_minimization.tpr -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx

gmx_mpi_d mdrun -v -deffnm step4.0_minimization
```

### 4.2 Equilíbrio

```bash
gmx_mpi_d grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_input.gro -p topol.top -n index.ndx

gmx_mpi_d mdrun -v -deffnm step4.1_equilibration
```

### 4.3 Produção

```bash
gmx_mpi_d grompp -f step5_production.mdp -c step4.1_equilibration.gro -p topol.top -n index.ndx -o step5_production.tpr

gmx_mpi_d mdrun -deffnm step5_production -nb gpu -gpu_id 0
```

### 4.4 Continuar simulação interrompida

```bash
gmx_mpi_d mdrun -v -deffnm step5_production -cpi step5_production.cpt -nb gpu -gpu_id 0
```
