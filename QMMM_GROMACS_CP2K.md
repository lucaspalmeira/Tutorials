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
resnr 12630 | resnr 12642 | resnr 12647 | resnr 12657 | resnr 12659 | resnr 12668 | resnr 12688 | resnr 12689

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

3. Arquivo .mdp da produção QM/MM (step6_qmmm.mdp)

```
; Copie TODAS as linhas de tcoupl, pcoupl, vdw, PME, etc. do seu step5_production.mdp
; (são exatamente as mesmas do CHARMM-GUI)

; ==================== QM/MM com CP2K ====================
qmmm-cp2k-active        = yes
qmmm-cp2k-qmgroup       = QMatoms
qmmm-cp2k-qmcharge      = 0
qmmm-cp2k-qmmultiplicity = 1
qmmm-cp2k-qmmethod      = PBE               ; padrão rápido e robusto
; qmmm-cp2k-qmmethod    = INPUT             ; descomente + use -qmi se quiser input CP2K customizado
```

4. Executar a produção QM/MM

```bash
# Escolher GPU 0 ou 1
export CUDA_VISIBLE_DEVICES=0   # ou 1

gmx grompp -f step6_qmmm.mdp -o step6_qmmm.tpr -c step4.1_equilibration.gro -t step4.1_equilibration.cpt -p topol.top -n index.ndx

# Execução (escolha uma das opções)

# Opção A – uma única RTX 4090 (mais simples e estável)
gmx mdrun -v -deffnm step6_qmmm -gpu_id 0 -ntomp 28

# Opção B – duas RTX 4090 (mais rápido se CP2K foi compilado com CUDA+DBCSR_ACC)
# mpirun -np 2 gmx_mpi mdrun -deffnm step6_qmmm -ntomp 14 -gpu_id 01
```

Continuar uma simulação interrompida

```bash
gmx mdrun -v -deffnm step6_qmmm -cpi step6_qmmm.cpt -gpu_id 0
```
