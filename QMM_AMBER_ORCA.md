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
pmemd.cuda -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5.mdout -r step5.rst7 -inf step5.mdinfo -x step5.nc
```

---

# 4. RESUMO DOS COMANDOS

### Minimização

```bash
sander -O -i step4.0_minimization.mdin -p step3_input.parm7 -c step3_input.rst7 \
 -o step4.0_minimization.mdout -r step4.0_minimization.rst7 -inf step4.0_minimization.mdinfo \
 -ref step3_input.rst7
```

### Equilibração

```bash
sander -O -i step4.1_equilibration.mdin -p step3_input.parm7 -c step4.0_minimization.rst7 \
 -o step4.1_equilibration.mdout -r step4.1_equilibration.rst7 -inf step4.1_equilibration.mdinfo \
 -ref step3_input.rst7 -x step4.1_equilibration.nc
```

### Produção (única etapa)

```bash
sander -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 \
 -o step5.mdout -r step5.rst7 -inf step5.mdinfo -x step5.nc
```

---

# 5. EXECUTAR TUDO EM GPU

Basta trocar **sander** por:

```
pmemd.cuda
```

ou, para paralelismo multi-GPU/NVLink:

```
pmemd.cuda.MPI
```
