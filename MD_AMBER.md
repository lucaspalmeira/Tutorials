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

## Geração de `dihe.restraint` para glicanos (inulina) no AMBER

Esta etapa descreve **passo a passo** como identificar um ligante glicano preparado pelo **CHARMM-GUI (Solution Builder)** e gerar manualmente um arquivo **`dihe.restraint`** para uso em simulações de dinâmica molecular no **AMBER**.

O procedimento é especialmente útil quando o CHARMM-GUI **não cria automaticamente** restrições diedrais para carboidratos/polímeros (caso comum para frutanos/inulina). Caso o arquivo tenha sido criado, pule ignore esta etapa e pule para a etapa de minimização.

---

### Contexto

* Sistema preparado no **CHARMM-GUI**
* Dinâmica será executada no **AMBER**
* Ligante: **inulina (β-frutano)**
* Resíduos do ligante: `1CU` e `0CU`
* Arquivos principais:

  * `amber/step3_input.parm7`
  * `amber/step3_input.rst7`

---

### Identificar os resíduos do ligante

Dentro do diretório do sistema, gere a lista de resíduos:

```bash
cpptraj step3_input.parm7 << EOF > residues.dat
resinfo :*
EOF
```

No arquivo `residues.dat`, identifique o ligante. Exemplo:

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

Átomos relevantes para diedros glicosídicos (exemplo real):

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

#### 🔹 Diedro φ (phi)

```
O5(i) – C2(i) – O1(i+1) – C1(i+1)
```

#### 🔹 Diedro ψ (psi)

```
C2(i) – O1(i+1) – C1(i+1) – C2(i+1)
```

Esses são os **únicos diedros que devem ser restringidos**.

**Nota:** Nunca restrinja diedros internos do anel.

---

### Medir os valores iniciais dos diedros

Crie o arquivo `get_glycan_dihes.cpptraj`:

```cpptraj
parm step3_input.parm7
trajin step3_input.rst7 1 1

dihedral phi_571_572 :571@O5 :571@C2 :572@O1 :572@C1 out phi_571_572.dat
dihedral psi_571_572 :571@C2 :572@O1 :572@C1 :572@C2 out psi_571_572.dat

run
```

Execute:

```bash
cpptraj -i get_glycan_dihes.cpptraj
```

Verifique os valores:

```bash
head phi_571_572.dat
head psi_571_572.dat
```

---

### Obter os índices absolutos dos átomos (ParmEd)

Entre no ParmEd:

```bash
parmed step3_input.parm7
```

Liste os átomos envolvidos:

```parmed
printAtoms :571@O5,C2
printAtoms :572@O1,C1,C2
```

Exemplo de saída:

```
8795 O5
8794 C2
8835 O1
8832 C1
8815 C2
```

Esses números serão usados no `dihe.restraint`.

---

### Criar o arquivo `dihe.restraint`

Exemplo **correto e funcional**:

```text
&rst
 iat=8795,8794,8835,8832,
 r1=-180.0, r2=-75.0, r3=-55.0, r4=180.0,
 rk2=20.0, rk3=20.0,
/

&rst
 iat=8794,8835,8832,8815,
 r1=-180.0, r2=100.0, r3=130.0, r4=180.0,
 rk2=20.0, rk3=20.0,
/
```

Ajuste `r2` e `r3` com base nos valores medidos (±10° é o ideal).

Repita para todas as ligações:

* 571–572
* 572–573
* 573–574
* 574–575
* 575–576

Para o resíduo terminal (`0CU`), aplique apenas os diedros possíveis.

---

### Ativar o `dihe.restraint` no AMBER

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

* ✔ Use `dihe.restraint` **somente até o fim da equilibration**
* ✔ Produção → **remova completamente**
* ✔ `rk2 = rk3 = 10–20` é ideal
* ✔ Nunca restrinja diedros do anel
* ✔ Método compatível com literatura de MD de glicanos

---

Este procedimento garante:

* Estabilidade conformacional inicial do glicano
* Evita colapsos não físicos

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

> As instruções acima foram extraídas e interpretadas do manual do Amber 2025 (https://ambermd.org/doc12/Amber25.pdf)
