# Simulation Scripts

This page documents the project scripts used to automate AMBER production runs and cpptraj analysis across replicas.

## `run_md_amber_three_replicates.sh`

This script is intended to run minimization, equilibration, and three production replicas in AMBER with `pmemd.cuda`.

### Typical Responsibilities

- prepare or reuse restart files from previous stages
- run staged minimization
- run equilibration
- launch three independent production replicas
- keep output names consistent across replicas
- preserve logs for debugging

### Expected Outputs

- minimization restart files
- equilibration restart files
- production trajectory files
- production restart files
- AMBER output logs for each stage and replica

## `run_cpptraj_replicatas.sh`

This script is intended to run post-production analysis with **cpptraj** across multiple replicas.

### Typical Analyses

- trajectory centering
- ligand RMSD
- residue RMSF
- radius of gyration
- hydrogen-bond analysis
- conformational clustering

### Expected Outputs

- processed trajectory files
- tabular analysis data
- clustering outputs
- analysis logs

## Adaptation Notes

Before running either script, confirm that file names, topology paths, replica labels, GPU settings, and analysis masks match the current project. Keep these assumptions visible near the top of each script.
