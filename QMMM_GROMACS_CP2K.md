# Simulação QM/MM (GROMACS 2021.3 + CP2K)

**Região QM (224 átomos):**
- Ligante: resíduos 603 (AGLC) + 604 (BFRU)
- Resíduos da proteína: 39 (ASP), 57 (LEU), 96 (PHE), 97 (ASP), 166 (ALA), 167 (PHE), 246 (GLU), 247 (THR)
- 8 moléculas de água: 12630, 12642, 12647, 12657, 12659, 12668, 12688, 12689

### 1. Criar o grupo de índice QM (QMatoms) – FAZER SÓ UMA VEZ

```bash
gmx make_ndx -f step3_input.gro -o index.ndx
```

No prompt interativo >, cole exatamente estas linhas, uma por uma (Enter após cada uma):

```
# Ligante completo (cadeia do ligante, resíduos 603 e 604)
r 603 | r 604

# Resíduos da proteína
r 39 | r 57 | r 96 | r 97 | r 166 | r 167 | r 246 | r 247

# 8 águas selecionadas
r 12893 | r 12910 | r 16456 | r 17118 | r 17223 | r 17229 | r 17305 | r 17311

# Junta tudo num único grupo
0 | 1 | 2

# Nomeia o grupo QM
name 3 QMatoms

# Sai e salva
q
```

2. Minimização + Equilibração clássica (igual ao CHARMM-GUI)

```bash
gmx grompp -f step4.0_minimization.mdp -o step4.0_minimization.tpr -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx -maxwarn 1
gmx mdrun -v -deffnm step4.0_minimization -gpu_id 0

gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_input.gro -p topol.top -n index.ndx
gmx mdrun -v -deffnm step4.1_equilibration -gpu_id 0
```

3. Arquivo .mdp da produção QM/MM (step5_production.mdp)

Adicione as seguintes linhas ao final do arquivo 'step5_production.mdp'

```
... Linhas anteriores presentes no step5_production.mdp
; CP2K QMMM parameters
qmmm-cp2k-active        = yes
qmmm-cp2k-qmgroup       = QMatoms ; Index group of QM atoms
qmmm-cp2k-qmmethod      = PBE
qmmm-cp2k-qmcharge      = 0
qmmm-cp2k-qmmultiplicity = 1
;
```

4. Executar a produção QM/MM

```bash
# Escolher GPU 0 ou 1
export CUDA_VISIBLE_DEVICES=0   # ou 1

gmx grompp -f step5_production.mdp -o step5_production.tpr -c step4.1_equilibration.gro -t step4.1_equilibration.cpt -p topol.top -n index.ndx

# Execução (escolha uma das opções)

# Opção A – uma única RTX 4090 (mais simples e estável)
gmx mdrun -v -deffnm step5_production -gpu_id 0 -ntomp 14

# Opção B – duas RTX 4090 (mais rápido se CP2K foi compilado com CUDA+DBCSR_ACC)
mpirun -np 2 gmx_mpi mdrun -deffnm step5_production -ntomp 14 -gpu_id 01
```

Continuar uma simulação interrompida

```bash
gmx mdrun -v -deffnm step6_qmmm -cpi step6_qmmm.cpt -gpu_id 0
```
