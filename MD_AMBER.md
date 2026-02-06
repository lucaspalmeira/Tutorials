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

### Definição correta dos diedros para β(2→1)-frutano (inulina)

Para cada ligação glicosídica entre os resíduos *i* (doador) e *i+1* (aceptor):

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

#### Em caso da ligação ser do tipo β(2→6)-frutano (Levan):

O carbono anomérico C2 do resíduo terminal (0CU) liga-se ao O6 do resíduo anterior (6CU).

O eixo glicosídico é definido por:

```
O5(i) – C6(i) – O6(i) – C2(i+1)
```

Somente dois diedros devem ser controlados:

#### Diedro φ (phi)
```
O5(i) – C6(i) – O6(i) – C2(i+1)
```

#### Diedro ψ (psi)
```
C6(i) – O6(i) – C2(i+1) – O5(i+1)
```

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
autoimage anchor :1-XXX
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
donormask :1CU,0CU \
acceptormask :1-XXX
EOF
```

Se quiser considerar ligações de hidrogênio entre solvente e ligante bem como solvente e proteína, use:

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

#### Cálculo de Energia Livre de Ligação por MMPBSA (AMBER)
Visão Geral do Método

A técnica MM/PBSA (ou MM/GBSA) permite estimar a energia livre de ligação de um complexo receptor-ligante a partir de snapshots de uma simulação de dinâmica molecular. Usaremos o script MMPBSA.py do AmberTools para calcular a energia livre de binding (ΔG_bind) e também decompor as contribuições energéticas por resíduo. Esse procedimento requer três topologias separadas (complexo, receptor e ligante) derivadas da topologia do sistema completo, além de um arquivo de trajetória e um arquivo de controle de parâmetros para o cálculo.


##### Preparação das Topologias “Secas” (sem solvente)

Antes de executar o MMPBSA, precisamos gerar arquivos de topologia (.prmtop) para: (1) o complexo proteína-ligante sem solvente, (2) a proteína (receptor) isolada, e (3) o ligante isolado. Como você já possui a topologia completa do sistema solvado (step3_input.parm7), utilizaremos o utilitário ante-MMPBSA.py para criar essas topologias filtradas, removendo água e íons e separando receptor e ligante.

*Passo 1* – Identificar o ligante no sistema: Abra o PDB ou a topologia do complexo e identifique quais resíduos correspondem ao ligante. Por exemplo, se o seu ligante for uma molécula distinta das cadeias proteicas (como uma pequena molécula ou polímero de açúcar), descubra o intervalo de resíduos ou o nome do residuário do ligante. Suponha, por exemplo, que a proteína contenha resíduos 1-630 e o ligante seja o resíduo 631 (ou um pequeno intervalo no final da numeração).

*Passo 2* – Executar o ante-MMPBSA.py: Com a informação acima, use o comando ante-MMPBSA.py para gerar três topologias secas. Por exemplo:

```bash
ante-MMPBSA.py -p step3_input.parm7 \
               -c complex.prmtop \
               -r receptor.prmtop \
               -l ligand.prmtop \
               -m ':1-630' \
               -s ':WAT,Na+,Cl-' \
               --radii=mbondi2
```

No exemplo acima, assumimos que :1-630 define a máscara Amber dos resíduos do receptor (proteína), de forma que o restante dos átomos (após remover solvente) será tratado como ligante. Ajuste a máscara -m de acordo com o seu sistema (pode ser um intervalo de resíduos ou identificador de cadeia do receptor). A opção -s ':WAT,Na+,Cl-' indica os átomos a remover (strip) da topologia original – aqui removemos moléculas de água (WAT) e íons sódio e cloreto. A opção --radii=mbondi2 define os raios atômicos (modelo Bondi modificado) adequados para cálculos de PB. Com este comando único, serão gerados os arquivos complex.prmtop (complexo sem solvente), receptor.prmtop (proteína) e ligand.prmtop (ligante). Dica: Verifique se o script identificou corretamente o receptor e ligante; ele normalmente imprime a suposição de máscara, mas usando a opção -m garantimos a seleção correta.


##### Configurando o Arquivo de Input do MMPBSA

Crie um arquivo de texto (por exemplo, mmpbsa.in) que define os parâmetros do cálculo MM/PBSA. O formato segue o estilo dos arquivos de entrada do sander do AMBER, contendo vários namelists começando com & e terminando com /. Para nosso objetivo, incluiremos seções para executar tanto o cálculo GB quanto PB em um único passo, e ativaremos a decomposição por resíduo. Abaixo um exemplo de input:

```
&general
   interval=1,           # usar todos os frames; ajuste se quiser pular frames
   verbose=1, 
   keep_files=0,         # 0 para descartar arquivos temporários (_MMPBSA_*), ou 2 para mantê-los
   strip_mask=":WAT,Na+,Cl-"
/
&gb
   igb=5,                # modelo Generalized Born (OBC II):contentReference[oaicite:5]{index=5}
   saltcon=0.100         # concentração salina 0.1 M
/
&pb
   istrng=0.100          # força iônica 0.1 M no cálculo PB:contentReference[oaicite:6]{index=6}
/
&decomp
   idecomp=1,            # decomposição por resíduo (1 ou 2 são per-residue):contentReference[oaicite:7]{index=7}
   dec_verbose=1         # informações detalhadas de decomposição
/
```

Explicação dos parâmetros principais:

* *&general*: Aqui definimos configurações gerais. interval=1 indica que cada frame da trajetória será usado (poderia ajustar, e.g. interval=5 para usar um a cada 5 frames). endframe poderia ser usado para limitar até um frame específico, mas caso não seja definido, o padrão é usar todos os frames disponíveis. A opção strip_mask=":WAT,Na+,Cl-" é uma camada extra de segurança para ignorar solvente/íons da trajetória (o script já os removerá usando as topologias secas, mas incluímos por garantia). keep_files=0 evita salvar arquivos temporários extensos (MMPBSA*), deixando apenas os resultados finais. verbose=1 aumenta a verbosidade das saídas.

* *&gb*: Configura o cálculo MM-GBSA (solvente implícito Generalized Born). Aqui escolhemos igb=5 (modelo GB OBC(II)) com saltcon=0.1 (equivalente a solução de 0,1 M de sal). (Nota: igb=2 (OBC(I)) também é comum; ambos são modelos GB disponíveis.)

* *&pb*: Configura o cálculo MM-PBSA (solvente implícito Poisson-Boltzmann). Usamos istrng=0.100 para definir força iônica de 0,1 M no solver PB. (Outros parâmetros PB, como tamanho de malha, usam padrões internos do Amber PBSA.)

* *&decomp*: Ativa a decomposição de energia por resíduo. Definimos idecomp=1 para decomposição por resíduo padrão (que inclui interações 1-4 de forma separada). A configuração dec_verbose=1 faz o script imprimir detalhes de cada componente de energia por resíduo no arquivo de decomposição. (Valores idecomp=1 ou 2 são usados para decomposição por resíduo; 2 agrupa as interações 1-4 dentro dos termos eletrostático e van der Waals. Ambos fornecem resultados de contribuição por resíduo semelhantes.) Observação: Devido a limitações do Amber, o termo de solvação não polar no PB não é decomposto por resíduo – esse componente aparecerá zero ou ausente na lista por resíduo para PB (enquanto para GB todos os termos são decompostos normalmente).

Com o arquivo de input pronto, revise para garantir que as máscaras (se usadas) correspondem ao seu sistema. No nosso exemplo, fornecemos strip_mask global para remover solvente. O script MMPBSA.py também tenta adivinhar máscaras de receptor e ligante a partir das topologias fornecidas, mas como geramos topologias separadas, essa etapa interna serve apenas para referência.


##### Executando o Cálculo MMPBSA

Agora podemos rodar o cálculo de duas formas: em modo serial (um processo) ou em modo paralelo (MPI), já que você possui a versão MPI (MMPBSA.py.MPI) instalada. A escolha depende do tamanho da trajetória e dos recursos disponíveis – para 100 ns de simulação, o modo paralelo pode acelerar significativamente o processamento dividindo os frames entre vários núcleos.

Modo Serial (único processo)

Para executar em modo serial, use o script MMPBSA.py diretamente. Exemplo de comando no terminal dentro do diretório de trabalho:

```bash
MMPBSA.py -O -i mmpbsa.in \
          -o FINAL_RESULTS_MMPBSA.dat \
          -do FINAL_DECOMP_MMPBSA.dat \
          -sp step3_input.parm7 \
          -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop \
          -y step5_production.nc
```


Descrição dos argumentos: -O permite sobrescrever arquivos de saída existentes. -i especifica o arquivo de input que criamos (mmpbsa.in). -o define o nome do arquivo de resultados de energia média (nesse caso, escolhemos FINAL_RESULTS_MMPBSA.dat). -do define o nome do arquivo de saída da decomposição por resíduo (listaremos as contribuições energéticas de cada resíduo nesse arquivo). Em seguida, -sp aponta para a topologia do sistema solvado original (usamos step3_input.parm7 como solvated complex topology para que o script saiba lidar com a trajetória contendo solvente). As opções -cp, -rp e -lp fornecem as topologias secas do complexo, receptor e ligante, respectivamente, que geramos no passo anterior. Finalmente, -y aponta para o arquivo de trajetória da simulação (step5_production.nc).

*Dica*: O MMPBSA.py aceita vários arquivos de trajetória ou uso de curingas (*). Se sua simulação estiver dividida em múltiplos arquivos (por exemplo, prod1.nc, prod2.nc, ...), você pode usar -y "*.nc" para ler todos. No nosso caso, temos um único step5_production.nc. Certifique-se de que o AmberTools reconhece o formato .nc (NetCDF); versões atuais do AmberTools geralmente suportam NetCDF nativamente.

Ao rodar o comando acima, o script executará iterativamente cálculos de energia para cada frame (usando sander internamente). Mensagens de progresso serão exibidas no terminal (stdout) e eventuais avisos/erros em stderr. Ao finalizar, você terá os arquivos de saída especificados.

Modo Paralelo (MPI)

Para aproveitar seu CPU multi-core (Ryzen 9 7950X3D) e acelerar o cálculo, você pode rodar o MMPBSA em paralelo. Certifique-se de ter um ambiente MPI configurado (por exemplo, OpenMPI ou MPICH) e use o executável MMPBSA.py.MPI. O comando é similar, mas prefixado por mpirun -np <N> onde <N> é o número de processos desejado. Por exemplo, para rodar com 16 núcleos em paralelo:

```bash
mpirun -np 16 MMPBSA.py.MPI -O -i mmpbsa.in \
          -o FINAL_RESULTS_MMPBSA.dat \
          -do FINAL_DECOMP_MMPBSA.dat \
          -sp step3_input.parm7 \
          -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop \
          -y step5_production.nc > mmpbsa_progress.log 2>&1
```

No exemplo acima, redirecionamos a saída para mmpbsa_progress.log para guardar o log (recomendado, pois em execução paralela as mensagens podem ser extensas). Você pode acompanhar esse log para ver o andamento. O resultado final (arquivos .dat) será o mesmo produzido no modo serial.

Observações sobre paralelização: O script MMPBSA.py.MPI divide os frames da trajetória entre os processos MPI para calcular as energias em paralelo. Ele funciona de forma mais eficiente quando o número de frames é múltiplo do número de processos, distribuindo carga igual entre eles. Porém, isso não é estritamente necessário – caso não seja múltiplo, o programa distribui os frames restantes entre os processos iniciais. Importante: não use mais processos do que o número de frames, pois isso causaria erro. Em geral, com ~100 ns de simulação, você terá milhares de frames, então utilizar 16 ou até 32 processos deve ser viável dado que você possui 32 threads de CPU. Ajuste -np conforme seus recursos e note que o uso de muitos processos aumentará a memória usada (cada processo carrega uma cópia parcial dos dados). Seu hardware de 128 GB RAM deve suportar confortavelmente esse cálculo.

##### Análise dos Resultados

Após a execução, o arquivo FINAL_RESULTS_MMPBSA.dat conterá um resumo das componentes de energia média do complexo, receptor e ligante, bem como as diferenças calculadas (Complexo - Receptor - Ligante) que correspondem à energia livre de ligação estimada. Por exemplo, você verá seções para GB e PB separadamente (pois calculamos ambos). Cada seção lista termos como energia de van der Waals (VDWAALS), eletrostática de Coulomb (EEL), solvatção polar (EGB para GB, EPB para PB) e solvatção não polar (ESURF), com seus valores médios e desvios padrão. A última parte de cada seção mostrará as Diferenças (Complex - Receptor - Ligand) para cada termo e a soma total (TOTAL), que é a ΔG_binding estimada. Uma ΔG negativa indica uma ligação favorável (espontânea), enquanto positiva indica desfavorável. Lembre-se de que não incluímos entropia conformacional nesse cálculo, então o valor obtido corresponde principalmente à diferença de entalpia livre; para estimar a energia livre absoluta de ligação, normalmente seria necessário incluir o termo entrópico (via análise de modos normais, etc.), mas isso pode ser omitido para comparar compostos ou mutantes.

No arquivo FINAL_DECOMP_MMPBSA.dat (gerado porque usamos -do e ativamos &decomp), você encontrará as contribuições energéticas por resíduo. Cada linha tipicamente corresponde a um resíduo do receptor ou do ligante, indicando, por exemplo, energia van der Waals, eletrostática, solvatção polar e não polar associadas àquele resíduo para a energia de interação. Esses valores representam quanto cada resíduo contribui para a energia de ligação (valores negativos significam que o resíduo favorece a ligação, positivos desfavorecem). Como mencionado, a contribuição não polar de PB pode não aparecer por resíduo devido a limitações do Amber, mas o efeito desse termo geralmente é pequeno e pode ser considerado de forma global. Use esse arquivo para identificar quais aminoácidos do receptor (ou átomos do ligante) são mais importantes energeticamente para a interação.

Por fim, você terá realizado com sucesso o cálculo de energia livre de binding tanto via MM-GBSA quanto MM-PBSA, em modo serial ou paralelo, e obtido a decomposição por resíduo para uma análise detalhada. Com esses resultados, é possível inferir a estabilidade relativa do complexo proteína-ligante e os hot-spots de interação ao nível de resíduo.

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
