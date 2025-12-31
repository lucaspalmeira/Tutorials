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

> As instruções acima foram extraídas e interpretadas do manual do Amber 2025 (https://ambermd.org/doc12/Amber25.pdf)
