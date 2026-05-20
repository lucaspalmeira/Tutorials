# Molecular Dynamics with GROMACS

This tutorial describes a protein-ligand molecular dynamics workflow with **GROMACS**, focusing on energy minimization, equilibration, production, trajectory preprocessing, structural analysis, and free-energy estimation with **gmx_MMPBSA**.

The workflow assumes that the system topology, box, solvation, and ion placement may have been prepared with CHARMM-GUI or an equivalent preparation pipeline.

## 1. Basic Concepts

Molecular dynamics uses Newtonian mechanics to simulate the time evolution of atoms. A typical GROMACS workflow contains:

1. system preparation: topology, box, solvation, and ions
2. energy minimization
3. equilibration, usually NVT and/or NPT
4. production dynamics
5. trajectory preprocessing and analysis

This tutorial focuses on steps 2 through 5.

## 2. GROMACS File Types

| File | Meaning | Use |
| --- | --- | --- |
| `.gro` | Structure file | Coordinates, atom types, and simulation box |
| `.top` | Topology file | Atom types, force-field parameters, connectivity, and included `.itp` files |
| `.mdp` | Molecular dynamics parameter file | Integrator, temperature, pressure, restraints, and simulation length |
| `.ndx` | Index file | Custom atom groups for restraints, analysis, and coupling |
| `.tpr` | Portable binary run-input file | Output of `gmx grompp`; executed by `gmx mdrun` |
| `.edr` | Energy file | Potential energy, kinetic energy, temperature, pressure, and related data |
| `.xtc` | Compressed trajectory | Coordinates over time, commonly used for structural analyses |
| `.trr` | Full trajectory | Coordinates, velocities, and forces |
| `.cpt` | Checkpoint file | Restart or continue interrupted simulations |
| `.psf` | Protein Structure File | Common in CHARMM/NAMD and sometimes present in hybrid workflows |

## 3. Energy Minimization

### Objective

Remove bad contacts, steric clashes, and geometric strain before molecular dynamics.

### Preprocess

```bash
gmx grompp -f step4.0_minimization.mdp -o step4.0_minimization.tpr -c step3_input.gro -r step3_input.gro -p topol.top -n index.ndx
```

Flags:

- `-f`: minimization `.mdp` file
- `-o`: output `.tpr`
- `-c`: initial structure
- `-r`: reference structure for restraints
- `-p`: topology
- `-n`: index file

### Run

```bash
gmx mdrun -v -deffnm step4.0_minimization
```

Flags:

- `-v`: verbose output
- `-deffnm`: base name for output files

Position restraints are common at this stage, especially for ligands.

## 4. Equilibration

### Objective

Adapt the system to the target thermodynamic conditions, such as temperature and pressure, while maintaining stability.

Equilibration is commonly split into NVT and NPT phases. In the example below, equilibration is represented as a single stage.

### Preprocess

```bash
gmx grompp -f step4.1_equilibration.mdp -o step4.1_equilibration.tpr -c step4.0_minimization.gro -r step3_input.gro -p topol.top -n index.ndx
```

### Run

```bash
gmx mdrun -v -deffnm step4.1_equilibration
```

During equilibration, position restraints are often used for sensitive parts of the system, such as the ligand or backbone.

## 5. Production

### Objective

Generate the final trajectory used for scientific analysis. Production should normally have no restraints, or only the restraints justified by the study design.

### Preprocess

```bash
gmx grompp -f step5_production.mdp -o step5_production.tpr -c step4.1_equilibration.gro -p topol.top -n index.ndx
```

### Run

```bash
gmx mdrun -v -ntomp 12 -deffnm step5_production -nb gpu -gpu_id 0
```

### Continue an Interrupted Run

```bash
gmx mdrun -v -ntomp 12 -deffnm step5_production -cpi -nb gpu -gpu_id 0
```

Flags:

- `-ntomp 12`: number of OpenMP threads
- `-cpi`: continue from checkpoint
- `-nb gpu`: run nonbonded interactions on GPU
- `-gpu_id 0`: select GPU 0

## 6. Production Output Files

- `.xtc`: trajectory for analysis
- `.edr`: energy data
- `.log`: simulation log
- `.cpt`: checkpoint
- `.gro`: final structure

## 7. Trajectory Preprocessing

Preprocess the trajectory before post-MD analysis or free-energy calculations. This avoids artifacts from periodic boundary conditions.

Main files:

- `traj_comp.xtc`
- `step5_production.tpr`
- `index.ndx`

### 7.1 Remove PBC and Rebuild Molecules

```bash
gmx trjconv -s step5_production.tpr -f traj_comp.xtc -o traj_noPBC.xtc -pbc mol
```

Select group `0` for `System`.

### 7.2 Center the System

```bash
gmx trjconv -s step5_production.tpr -f traj_noPBC.xtc -o traj_center.xtc -center
```

Select `Protein`, then `System`.

### 7.3 Rotational and Translational Fitting

```bash
gmx trjconv -s step5_production.tpr -f traj_center.xtc -o traj_fit.xtc -fit rot+trans
```

Select `Backbone`, then `System`.

## 8. Post-MD Analyses

Use the preprocessed trajectory, usually `traj_fit.xtc`.

### 8.1 RMSD

RMSD evaluates structural stability over time relative to a reference.

```bash
gmx rms -s step5_production.tpr -f traj_fit.xtc -n index.ndx -o rmsd.xvg
```

For protein RMSD, select `Protein` twice. For ligand RMSD, select the ligand group twice.

### 8.2 RMSF

RMSF measures average residue flexibility.

```bash
gmx rmsf -s step5_production.tpr -f traj_fit.xtc -n index.ndx -o rmsf_residue.xvg -res
```

Select `Protein`.

### 8.3 Radius of Gyration

Radius of gyration evaluates structural compactness over time.

```bash
gmx gyrate -s step5_production.tpr -f traj_fit.xtc -n index.ndx -o gyrate.xvg
```

Select `Protein`.

### 8.4 Hydrogen Bonds

```bash
gmx hbond -s step5_production.tpr -f traj_fit.xtc -n index.ndx -num hbonds.xvg
```

Select `Protein`, then the ligand group.

## 9. Free-Energy Estimation with gmx_MMPBSA

`gmx_MMPBSA` estimates binding free energies with MM/PBSA or MM/GBSA.

### Enthalpy

Create `mmpbsa_enthalpy.in`:

```text
Sample input file for PB calculation
This input file is meant to show only that gmx_MMPBSA works. Although, we tried to use the input files as recommended
in the Amber manual, some parameters have been changed to perform more expensive calculations in a reasonable amount of time.
Feel free to change the parameters according to what is better for your system.

&general
sys_name="Prot-Lig-CHARMM",
startframe=1,
endframe=9999999999,
interval=1,
# In gmx_MMPBSA v1.5.0 we have added a new PB radii set named charmm_radii.
# This radii set should be used only with systems prepared with CHARMM force fields.
# Uncomment the line below to use charmm_radii set
#PBRadii=7,
/
&pb
# radiopt=0 is recommended, which means using radii from the prmtop file for both
# the PB calculation and the NP calculation.
istrng=0.15, fillratio=4.0, radiopt=0
/
```

Run:

```bash
mpirun -np 12 gmx_MMPBSA -O -i mmpbsa_enthalpy.in -cs step5_production.tpr -ct traj_fit.xtc -ci index.ndx -cg 1 13 -cp topol.top -o FINAL_RESULTS_MMPBSA_enthalpy.dat -eo FINAL_RESULTS_MMPBSA_enthalpy.CSV
```

### Entropy

Create `mmpbsa_entropy.in`:

```text
Sample input file for entropy calculations (IE)
This input file is meant to show only that gmx_MMPBSA works. Although,
we tried to use the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters
according to what is better for your system.

&general
sys_name="IE",
startframe=1,
endframe=99999999999,
# Interaction Entropy approximation
interaction_entropy=1, ie_segment=50, temperature=323.15
/

&gb
igb=2, saltcon=0.150,
/
```

Run:

```bash
mpirun -np 12 gmx_MMPBSA -O -i mmpbsa_entropy.in -cs step5_production.tpr -ct traj_fit.xtc -ci index.ndx -cg 1 13 -cp topol.top -o FINAL_RESULTS_MMPBSA_entropy.dat -eo FINAL_RESULTS_MMPBSA_entropy.CSV
```

### Per-Residue Energy Decomposition

Create `decomposition.in`:

```text
Sample input file for decomposition analysis
Make sure to include at least one residue from both the receptor
and ligand in the print_res mask of the &decomp section.
This is automatically guaranteed when using the "within" keyword.

&general
startframe=1, endframe=9999999999999, interval=1,
/

&gb
igb=5, saltcon=0.150,
/

&decomp
idecomp=2, dec_verbose=3,
print_res="within 6"
/
```

Run:

```bash
mpirun -np 12 gmx_MMPBSA -O -i decomposition.in -cs step5_production.tpr -ct traj_fit.xtc -ci index.ndx -cg 1 13 -cp topol.top -o FINAL_RESULTS_MMPBSA_decomposition.dat -eo FINAL_RESULTS_MMPBSA_decomposition.CSV
```
