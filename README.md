# Tutorials

Este repositório reúne tutoriais e scripts práticos para bioinformática estrutural, docking molecular, dinâmica molecular, QM/MM e rotinas auxiliares de uso em Linux/HPC.

## Conteúdo do repositório

### `Enzeptional_(GT4SD).md`
Tutorial para predição e otimização de enzimas com o framework **Enzeptional (GT4SD)**, com foco em viabilidade catalítica (*feasibility*) e `kcat` para reações específicas.

### `Hydrated_docking.md`
Guia de **docking molecular com águas explícitas**, usando **AutoDock Vina**, **Meeko** e **AutoGrid4**, com foco em *hydrated docking* no sítio de ligação.

### `MD_AMBER.md`
Tutorial completo de **dinâmica molecular clássica com AMBER** para sistemas proteína-ligante, incluindo preparação, minimização, equilibração, produção e análises com foco em carbohidratos.

### `MD_GROMACS.md`
Tutorial de **dinâmica molecular com GROMACS**, cobrindo minimização, equilibração, produção, pré-processamento de trajetória, análises estruturais e cálculo de energia livre com **gmx_MMPBSA**.

### `MD_SLURM.md`
Guia prático para submissão e monitoramento de **jobs no SLURM**, com foco em fluxos de trabalho com **GROMACS**, **AMBER** e **MMPBSA** em ambiente HPC.

### `QMMM_AMBER_QUICK.md`
Tutorial de **QM/MM no AMBER** usando **QUICK**, incluindo etapas de minimização, equilíbrio e produção com métodos semiempíricos e DFT.

### `QMMM_GROMACS_CP2K.md`
Protocolo para **QM/MM com GROMACS + CP2K**, incluindo definição da região quântica, criação de grupos de índice e execução em ambiente com GPU.

## Scripts

### `run_md_amber_three_replicates.sh`
Script para executar **minimização, equilibração e três replicatas de produção** no AMBER com `pmemd.cuda`.

### `run_cpptraj_replicatas.sh`
Script para análises pós-produção com **cpptraj** em múltiplas replicatas, incluindo:
- centralização da trajetória
- RMSD do ligante
- RMSF por resíduo
- raio de giro
- ligações de hidrogênio
- clustering conformacional

## Guias rápidos

### `commands-docker.md`
Resumo de **comandos úteis do Docker**, incluindo imagens, contêineres, volumes, limpeza de ambiente e backup/restauração.

### `git_commands.md`
Lista de **comandos essenciais do Git** para configuração, versionamento, branches, commits, inspeção e recuperação de alterações.

## Objetivo

O objetivo deste repositório é servir como uma coleção prática de protocolos e comandos para apoiar estudos em:

- docking molecular
- dinâmica molecular
- QM/MM
- uso de containers
- execução em ambientes HPC
- automação de análises

## Observação

Os tutoriais foram escritos com foco prático e podem exigir adaptação conforme:
- o sistema estudado
- os arquivos gerados pelo CHARMM-GUI
- o cluster/HPC utilizado
- as versões dos programas instalados

---
**Autor:** Lucas Palmeira
