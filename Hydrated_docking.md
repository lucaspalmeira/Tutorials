# Docking Molecular com Águas Explícitas (Hydrated Docking)  
Usando AutoDock Vina 1.2.x + Meeko + AutoGrid4

**Objetivo:** Realizar docking molecular considerando explicitamente moléculas de água no sítio de ligação  
**Método recomendado:** AutoDock4 scoring (`--scoring ad4`) + águas explícitas no ligante

## 1. Pré-requisitos e Instalação

### 1.1 Instalação do AutoDock Vina 1.2.x

```bash
# 1. Instalar compilador e dependências
# Ubuntu/Debian
sudo apt update
sudo apt install build-essential libboost-all-dev swig

# 2. Baixar e compilar Vina
git clone https://github.com/ccsb-scripps/AutoDock-Vina.git
cd AutoDock-Vina/build/linux/release
make

# 3. Colocar o executável no PATH (opcional mas recomendado)
sudo cp vina /usr/local/bin/
# ou adicione ao PATH:
# export PATH=$PATH:/caminho/para/AutoDock-Vina/build/linux/release
```

### 1.2 Scripts necessários

Baixe os scripts da Meeko e do mapwater:

```bash
https://github.com/forlilab/Meeko/release/meeko/cli/mk_prepare_receptor.py
https://github.com/forlilab/Meeko/release/meeko/cli/mk_prepare_ligand.py
https://github.com/forlilab/molscrub/blob/develop/scripts/scrub.py
https://github.com/ccsb-scripps/AutoDock-Vina/develop/example/autodock_scripts/mapwater.py
https://github.com/ccsb-scripps/AutoDock-Vina/blob/develop/example/autodock_scripts/dry.py
```

### 1.3 Outros programas necessários

- PDB2PQR (pH e protonação)
- PyMOL
- AutoGrid4 (vem com MGLTools ou AutoDockTools)

## 2. Preparação do Receptor

```bash
# 1. Corrigir e protonar o receptor em pH desejado (exemplo: pH 5.0)
pdb2pqr --ff=AMBER --ph-calc-method=propka --with-ph=5.0 --keep-chain --pdb-output receptor_pH5.0.pdb receptor_original.pdb receptor_pqr.pdb

# 2. Preparar o receptor para AutoDock (Meeko)
python3 mk_prepare_receptor.py -i receptor_pH5.0.pdb -o receptor --box_center 12.0 9.0 -1.0 --box_size 25.5 31.5 25.5 -p -g

# Arquivos gerados:
# - receptor.pdbqt
# - receptor.gpf
# - receptor.box.pdb (para visualização)
```

**Dica:** Abra `receptor.box.pdb` no PyMOL para conferir se a caixa está bem posicionada.

## 3. Preparação do Ligante

No Pymol:

```pymol
load ligante.pdb
save ligante.sdf
```

Em seguinda, no terminal:

```bash
scrub.py ligante.sdf -o liganteH.sdf
```

```bash
python3 mk_prepare_ligand.py -i liganteH.sdf -o ligante_hydrated.pdbqt -w
```

## 4. Geração dos Mapas de Afinidade (AutoGrid4)

```bash
# Verificar tipos de átomos do ligante
grep "^ATOM\|^HETATM" ligante_hydrated.pdbqt | cut -c 78-79 | sort | uniq

# Editar receptor.gpf → manter apenas os tipos que existem no ligante!
# Exemplo comum para moléculas orgânicas:
# ligand_types C OA HD H

# Rodar AutoGrid4
autogrid4 -p receptor.gpf -l receptor.glg
```

## 5. Criar o Mapa Especial de Água

```bash
python3 mapwater.py -r receptor.pdbqt -s receptor.W.map
```

## 6. Realizar o Docking (OBRIGATÓRIO: scoring ad4)

```bash
vina --ligand ligante.pdbqt --maps receptor --scoring ad4 --exhaustiveness 32 --out ligante_out_ad4.pdbqt
```

## 7. Pós-processamento – Classificar as águas

```bash

python3 dry.py -r receptor.pdbqt -m receptor.W.map -i ligante_out_ad4.pdbqt
```

Resultado: arquivo com sufixo `_DRY_SCORED.pdbqt` com as águas classificadas:

- STRONG   → água bem mantida
- WEAK     → água marginal
- DISPLC   → água provavelmente deslocada

## 8. Visualização Final (recomendado)

```bash
pymol receptor_pH5.0.pdb receptor.box.pdb ligante_out_ad4_DRY_SCORED.pdbqt
```

## Resumo dos Arquivos Mais Importantes

| Etapa                       | Arquivo principal gerado                     | Finalidade                                  |
|-----------------------------|----------------------------------------------|---------------------------------------------|
| Protonação                  | receptor_pH5.0.pdb, receptor_pqr.pdb         | Receptor corrigido e protonado              |
| Preparo receptor            | receptor.pdbqt, receptor.gpf                 | Arquivos AutoDock                           |
| Preparo ligante             | ligante_hydrated.pdbqt                       | Ligante + águas explícitas                  |
| Mapas                       | receptor.*.map + receptor.W.map              | Mapas de afinidade + mapa de água           |
| Docking                     | ligante_out_ad4.pdbqt                        | Resultados do docking                       |
| Pós-processamento           | ligante_out_ad4_DRY_SCORED.pdbqt             | Resultados finais com águas classificadas   |

## Considerações Finais

- Use **sempre** `--scoring ad4` para docking hidratado
- Aumente `--exhaustiveness` (32–64) para fragmentos ou ligantes flexíveis
- O método **não é recomendado** para virtual screening de milhares de compostos
