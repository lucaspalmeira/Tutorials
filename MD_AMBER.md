# Dinâmica molecular clássica (proteína-ligante) utilizando o AMBER

Para especificar a GPU device=1 (ou seja, a segunda GPU, já que a contagem começa em 0) ao rodar pmemd.cuda no AMBER, a forma padrão e recomendada é usar a variável de ambiente CUDA_VISIBLE_DEVICES. Atualmente, a seleção da GPU a ser usada em execuções com uma única GPU é automática se as GPUs estiverem configuradas no modo exclusivo de processo (nvidia-smi -c 3), mas a abordagem recomendada é usar a variável de ambiente CUDA_VISIBLE_DEVICES para selecionar qual GPU deve ser usada.

Desta forma, execute:

```bash
export CUDA_VISIBLE_DEVICES=1
```

ou

```bash
CUDA_VISIBLE_DEVICES="1" pmemd.cuda ...
```

### 1. Minimização

```bash
pmemd.cuda -O -i step4.0_minimization.mdin -p step3_input.parm7 -c step3_input.rst7 -o step4.0_minimization.mdout -r step4.0_minimization.rst7 -inf step4.0_minimization.mdinfo -ref step3_input.rst7
```

### 2. Equilibração

```bash
pmemd.cuda -O -i step4.1_equilibration.mdin -p step3_input.parm7 -c step4.0_minimization.rst7 -o step4.1_equilibration.mdout -r step4.1_equilibration.rst7 -inf step4.1_equilibration.mdinfo -ref step3_input.rst7 -x step4.1_equilibration.nc
```

### Produção

```bash
pmemd.cuda -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5_production.mdout -r step5_production.rst7 -inf step5_production.mdinfo -x step5_production.nc
```


