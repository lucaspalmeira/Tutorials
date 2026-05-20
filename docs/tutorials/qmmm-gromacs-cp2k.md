# QM/MM with GROMACS 2021.3 and CP2K

This protocol describes how to configure and run a QM/MM simulation with **GROMACS 2021.3** coupled to **CP2K**.

The example QM region contains:

- ligands: residues 603, glucose `AGLC`, and 604, fructose `BFRU`, forming sucrose
- active-site protein residues: 64 `ASP`, 82 `LEU`, 121 `PHE`, 122 `ASP`, 193 `ARG`, 194 `HIS`, 271 `GLU`, and 272 `THR`
- water molecules within 4 angstrom of the ligand, positionally restrained during minimization and equilibration

Prepare the protein-ligand complex with CHARMM-GUI using `Solution Builder`. Check pKa/protonation, apply the desired pH, use a rectangular box, add NaCl ions, select the CHARMM36 force field, and confirm the target temperature for equilibration and production. In this workflow, rectangular boxes are preferred because octahedral boxes have caused atom-overlap errors in some systems.

## 1. Create the QM Index Group

### 1.1 Create Basic Groups in `index.ndx`

```bash
gmx make_ndx -f step3_input.gro -n index.ndx
```

Inside `make_ndx`, run:

```bash
# Selected ligand atoms involved in the sucrose reaction
a 9068 | a 9090 | a 9066 | a 9098 | a 9099 | a 9100 | a 9101 | a 9102 | a 9067
name 4 LIG

# Selected active-site protein atoms
a 596 | a 597 | a 598 | a 3728 | a 3729 | a 3730 | a 3731
name 5 active_site_atoms

q
```

### 1.2 Select Water Oxygens within 4 Angstrom of the Ligand

```bash
gmx select -f step3_input.gro -s step4.0_minimization.tpr -select "name OH2 and (within 0.4 of group AGLC or within 0.4 of group BFRU)" -on wat_oxygens_near_lig.ndx

gmx select -f step3_input.gro -s step3_input.gro -select "resname TIP3 and same residue as (atomname OW and within 0.5 of resnr 246)" -on near246

gmx select -f structure.gro -s structure.gro -select " name OH2  and ( within 0.4 of resnr 246)" -on near246.ndx
```

The file `wat_oxygens_near_lig.ndx` should contain a group similar to:

```text
[ name_OH2_and_(within_0.4_of_group_AGLC_or_within_0.4_of_group_BFRU)_f0_t0.000 ]
44331 45630 45666 45681 56370 58383 58677 58692 58710 58767 58947 58965
```

### 1.3 Identify Water Residue Numbers with PyMOL

Open `step3_input.gro` in PyMOL and run:

```python
residues = []
iterate (index 44346 or index 45525 or index 45690 or index 45726 or index 45741 or index 56379 or index 58365 or index 58665 or index 58680 or index 58698 or index 58755 or index 58926 or index 58929 or index 58944), residues.append(resi)

unique_residues = sorted(list(set(residues)))
print(unique_residues)
cmd.select("waters_near_lig", "resi " + "+".join(map(str, unique_residues)))
```

Record the residue numbers of the nearby waters.

### 1.4 Add Nearby Waters to `index.ndx`

```bash
gmx make_ndx -f step3_input.gro -n index.ndx
```

Inside `make_ndx`, use the real water residue numbers obtained in PyMOL:

```bash
r 12439 | r 12872 | r 12884 | r 12889 | r 16452 | r 17123 | r 17221 | r 17226 | r 17232 | r 17251 | r 17311 | r 17317
name 6 WAT_4A
```

Create the final QM group:

```bash
4 | 5 | 6
name 7 QMatoms

q
```

## 2. Positional Restraints for Nearby Waters

Add the following only to the minimization and equilibration `.mdp` files, such as `step4.0_minimization.mdp` and `step4.1_equilibration.mdp`:

```mdp
; Restrain waters near the ligand to keep them in the QM region
pull                     = yes
pull-ngroups             = 2
pull-ncoords             = 1
pull-group1-name         = WAT_4A
pull-group2-name         = LIG
pull-coord1-type         = constraint
pull-coord1-geometry     = direction
pull-coord1-groups       = 1 2
pull-coord1-dim          = Y Y Y
pull-coord1-start        = yes
pull-coord1-rate         = 0
pull-coord1-k            = 1000
```

Do not include these lines in `step5_production.mdp`, because production should not restrain these waters unless the scientific design requires it.

## 3. QM/MM Parameters

Add the following to the end of each `.mdp` file used for minimization, equilibration, and production:

```mdp
; QM/MM parameters with CP2K
qmmm-cp2k-active         = yes
qmmm-cp2k-qmgroup        = QMatoms
qmmm-cp2k-qmmethod       = PBE
qmmm-cp2k-qmcharge       = 0
qmmm-cp2k-qmmultiplicity = 1
qmmm-cp2k-core           = auto
```

## 4. Run with Docker

```bash
# Pull the image
sudo docker pull kimjoochan/gromacs-cp2k:2022.2-9.1-cuda

# Start the container. Adjust the local path and GPU as needed.
sudo docker run --user $(id -u):$(id -g) -v "$(pwd):/home" -v "/etc/localtime:/etc/localtime:ro" -e TZ=America/Sao_Paulo -v "/tmp/.X11-unix:/tmp/.X11-unix" -e QT_X11_NO_MITSHM=1 -e OMPI_ALLOW_RUN_AS_ROOT=1 -e OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 --shm-size=1g --env="DISPLAY" --net=host --gpus '"device=1"' --privileged -itd --name gmx_cp2k kimjoochan/gromacs-cp2k:2022.2-9.1-cuda /bin/bash

# Check the GPU inside the container
docker exec gmx_cp2k nvidia-smi
```

Change `'"device=1"'` according to the available GPU.

### 4.1 Minimization

Enter the container:

```bash
sudo docker exec -it gmx_cp2k bash
cd home
```

Run:

```bash
gmx_mpi_d grompp -f step4.0_minimization.mdp -o step4.0_minimization.tpr -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx

gmx_mpi_d mdrun -v -deffnm step4.0_minimization
```

### 4.2 Equilibration

```bash
gmx_mpi_d grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_input.gro -p topol.top -n index.ndx

gmx_mpi_d mdrun -v -deffnm step4.1_equilibration
```

### 4.3 Production

```bash
gmx_mpi_d grompp -f step5_production.mdp -c step4.1_equilibration.gro -p topol.top -n index.ndx -o step5_production.tpr

gmx_mpi_d mdrun -deffnm step5_production -nb gpu -gpu_id 0
```

### 4.4 Continue an Interrupted Simulation

```bash
gmx_mpi_d mdrun -v -deffnm step5_production -cpi step5_production.cpt -nb gpu -gpu_id 0
```
