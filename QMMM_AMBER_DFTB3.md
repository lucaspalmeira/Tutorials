# DFTB3-Based QM/MM MD in AMBER: A Tutorial on Sucrose Hydrolysis by Invertase
## First step of the sucrose hydrolysis reaction by invertase

## Introduction

<p align="justify">
DFTB3 is an extension of SCC-DFTB, also known as DFTB2, based on the expansion of the total DFT energy up to third order. This improvement allows a more accurate description of charge variations, electrostatic interactions, and hydrogen bonds in biomolecular systems. Unlike SCC-DFTB, in which the Hubbard parameter is fixed for each element, DFTB3 allows chemical hardness to depend on the atomic charge state, making the method more robust for reactions involving electronic redistribution.
</p>

<p align="justify">
In this tutorial, the system considered is an invertase complexed with sucrose, previously prepared using CHARMM-GUI at pH 4.0. The simulation is performed in AMBER using the FF19SB force field for the classical region, at a temperature of 323.15 K, with a total simulation time of 1 ns. This preparation provides a biologically consistent starting structure, with protonation states adjusted to the acidic environment and simulation conditions compatible with the catalytic activity of the enzyme.
</p>

<p align="justify">
The reaction of interest is sucrose hydrolysis, characterized by the cleavage of the beta-glycosidic beta(2->1) bond between the glucose and fructose units. Since this process involves bond breaking and bond formation, a purely classical description is not sufficient. Therefore, the use of DFTB3 within a QM/MM approach is suitable for treating the reactive region, while the remaining protein, solvent, and ions can be described by molecular mechanics.
</p>

<p align="justify">
Despite its advantages, DFTB3 still has limitations, especially in systems containing sp3 nitrogen atoms, where proton affinities may be underestimated. Even so, for biomolecular systems involving charge transfer, hydrogen bonding, and enzymatic reactions, the method provides a good balance between computational cost and quantum mechanical description. Thus, it is particularly useful for investigating catalytic mechanisms in enzymes such as invertases, allowing the structural and electronic dynamics associated with sucrose hydrolysis to be monitored.
</p>

---

This document describes the full workflow used to analyze the reaction mechanism involving:

- nucleophilic attack of Asp38 on fructose C2
- proton donation from GLH245
- transition state stabilization by Asp168

The reactive fructose residue is **0CU**, and the reactive atom is:

C2 → 0CU@C2

Target bond formation:
Asp38@OD2 → 0CU@C2

---

### 1. Run Minimization

```bash
sed -e "s/FC/1.0/g" dihe.restraint > step4.1_equilibration.rest
```

```bash
mpirun -np 30 sander.MPI -O -i step4.0_minimization.mdin -p step3_input.parm7 -c step3_input.rst7 -o step4.0_minimization.mdout -r step4.0_minimization.rst7 -inf step4.0_minimization.mdinfo -ref step3_input.rst7
```

### 2. Run Equilibration

```bash
sed -e "s/FC/1.0/g" dihe.restraint > step4.0_minimization.rest
```

```bash
mpirun -np 30 sander.MPI -O -i step4.1_equilibration.mdin -p step3_input.parm7 -c step4.0_minimization.rst7 -o step4.1_equilibration.mdout -r step4.1_equilibration.rst7 -inf step4.1_equilibration.mdinfo -ref step3_input.rst7 -x step4.1_equilibration.nc
```

### 3. Run Production

```bash
mpirun -np 30 sander.MPI -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5_production.mdout -r step5_production.rst7 -inf step5_production.mdinfo -x step5_production.nc
```

Restart production
```bash
mpirun -np 30 sander.MPI -O -i step5_production.mdin -p step3_input.parm7 -c step5_production.rst7 -o step5_production_2.mdout -r step5_production_2.rst7 -inf step5_production_2.mdinfo -x step5_production_2.nc
```

---

### 4. Center trajctory

```bash
cpptraj -p step3_input.parm7 -y step5_production.nc << EOF
autoimage anchor :1-XXX
center :1-XXX mass origin
image origin center
rms first :1-XXX@CA
trajout step5_centered.nc

EOF
```

### 5. Extract representative cluster frame

File: extract_cluster_rep.in

```cpptraj
parm step3_input.parm7
trajin step5_centered.nc

trajout cluster2_rep_frame19.pdb pdb onlyframes 19

run
quit
```

Run:

```bash
cpptraj -i extract_cluster_rep.in
```

---

### 6. Measure mechanism distances

File: mechanism_distances.in

```cpptraj
parm step3_input.parm7
trajin step5_centered.nc

distance d_Asp38_OD2_C2_0CU :38@OD2 :0CU@C2 out d_Asp38_OD2_C2_0CU.dat
distance d_Asp38_CG_C2_0CU  :38@CG  :0CU@C2 out d_Asp38_CG_C2_0CU.dat
distance d_Asp38_OD1_C2_0CU :38@OD1 :0CU@C2 out d_Asp38_OD1_C2_0CU.dat

distance d_C2_0CU_O1_1GA :0CU@C2 :1GA@O1 out d_C2_0CU_O1_1GA.dat

distance d_GLH245_HE2_O1_1GA :245@HE2 :1GA@O1 out d_GLH245_HE2_O1_1GA.dat
distance d_GLH245_OE2_HE2 :245@OE2 :245@HE2 out d_GLH245_OE2_HE2.dat
distance d_GLH245_OE2_O1_1GA :245@OE2 :1GA@O1 out d_GLH245_OE2_O1_1GA.dat

distance d_Asp168_OD2_C2_0CU :168@OD2 :0CU@C2 out d_Asp168_OD2_C2_0CU.dat
distance d_Asp168_OD2_O1_1GA :168@OD2 :1GA@O1 out d_Asp168_OD2_O1_1GA.dat
distance d_Asp168_CG_C2_0CU  :168@CG  :0CU@C2 out d_Asp168_CG_C2_0CU.dat

run

writedata mechanism_distances.agr   d_Asp38_OD2_C2_0CU   d_Asp38_CG_C2_0CU   d_Asp38_OD1_C2_0CU   d_C2_0CU_O1_1GA   d_GLH245_HE2_O1_1GA   d_GLH245_OE2_HE2   d_GLH245_OE2_O1_1GA   d_Asp168_OD2_C2_0CU   d_Asp168_OD2_O1_1GA   d_Asp168_CG_C2_0CU

quit
```

Run:

```bash
cpptraj -i mechanism_distances.in
```

---

### 7. Measure distances only in representative frame

File: rep_frame_measure.in

```cpptraj
parm step3_input.parm7
trajin step5_centered.nc 19 19

distance d_Asp38_OD2_C2_0CU :38@OD2 :0CU@C2 out rep_d_Asp38_OD2_C2_0CU.dat
distance d_Asp38_CG_C2_0CU :38@CG :0CU@C2 out rep_d_Asp38_CG_C2_0CU.dat
distance d_C2_0CU_O1_1GA :0CU@C2 :1GA@O1 out rep_d_C2_0CU_O1_1GA.dat

distance d_GLH245_HE2_O1_1GA :245@HE2 :1GA@O1 out rep_d_GLH245_HE2_O1_1GA.dat
distance d_GLH245_OE2_HE2 :245@OE2 :245@HE2 out rep_d_GLH245_OE2_HE2.dat

distance d_Asp168_OD2_C2_0CU :168@OD2 :0CU@C2 out rep_d_Asp168_OD2_C2_0CU.dat

trajout cluster2_rep_frame10_measured.pdb pdb

run
quit
```

Run:

```bash
cpptraj -i rep_frame_measure.in
```

---

### 8. Python plotting script

File: plot_mechanism_distances.py

```python
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR = Path(".")
REP_FRAME = 10
OUTPUT_FIG = "mechanism_distances.png"

FILES_AND_LABELS = {
    "d_Asp38_OD2_C2_0CU.dat": "Asp38 OD2 - 0CU C2",
    "d_Asp38_CG_C2_0CU.dat": "Asp38 CG - 0CU C2",
    "d_Asp38_OD1_C2_0CU.dat": "Asp38 OD1 - 0CU C2",
    "d_C2_0CU_O1_1GA.dat": "0CU C2 - 1GA O1",
    "d_GLH245_HE2_O1_1GA.dat": "GLH245 HE2 - 1GA O1",
    "d_GLH245_OE2_HE2.dat": "GLH245 OE2 - HE2",
    "d_Asp168_OD2_C2_0CU.dat": "Asp168 OD2 - 0CU C2",
}

def read_cpptraj_dat(filepath):
    x,y=[],[]
    for line in open(filepath):
        if line.startswith("#") or not line.strip():
            continue
        f,v=line.split()[:2]
        x.append(float(f))
        y.append(float(v))
    return np.array(x),np.array(y)

plt.figure(figsize=(10,6))

for f,label in FILES_AND_LABELS.items():
    if Path(f).exists():
        x,y=read_cpptraj_dat(f)
        plt.plot(x,y,label=label)

plt.axvline(REP_FRAME,linestyle="--")

plt.xlabel("Frame")
plt.ylabel("Distance (Å)")
plt.legend()
plt.tight_layout()

plt.savefig(OUTPUT_FIG,dpi=600)
plt.show()
```

Run:

```bash
python plot_mechanism_distances.py
```

---

### 9. Expected outputs

cluster2_rep_frame10.pdb
mechanism_distances.agr
mechanism_distances.png

---

### 10. Interpretation guidelines

Covalent intermediate formation:

Asp38 OD2 – C2 distance ≈ 1.4–1.7 Å

Bond breaking:

C2 – O1 distance increases

Proton transfer:

GLH245 HE2 approaches O1

Transition-state stabilization:

Asp168 interacts with C2 or O1

---

### 11. Full workflow

```bash
cpptraj -i extract_cluster_rep.in
cpptraj -i mechanism_distancesd.in
cpptraj -i rep_frame_measure.in
python plot_mechanism_distances.py
```
# Second step of the sucrose hydrolysis reaction by invertase using QM/MM-DFTB3 in AMBER

## Formation of the ASJ/LIG complex and simulation of water attack on the fructosyl–Asp38 covalent intermediate

---

## 1. Purpose of this tutorial

<p align="justify">
This tutorial describes the complete workflow used to prepare and simulate the second step of the reaction catalyzed by invertase using QM/MM-DFTB3 in AMBER. The first reaction step consisted of the nucleophilic attack of Asp38 on the fructose C2 atom, forming a covalent intermediate between fructose and the catalytic Asp38 residue. The step described here starts from this covalent intermediate and prepares the system to investigate the attack of a water molecule on the covalent bond between Asp38 and fructose, allowing recovery of Asp38 and release of fructose in the active site. <strong>Frame 14</strong> from the trajectory obtained in the first reaction step was selected because, in this frame, the covalent bond <code>Asp38@OD2–fructose@C2</code> had a length of ~1.4 Å, which is consistent with the literature (JITONNOM, KETUDAT-CAIRNS, HANNONGBUA, 2018).
</p>

<p align="justify">
The water molecule used for the nucleophilic attack was selected from <strong>frame 14</strong> of the dynamics obtained in the first reaction step. This frame was chosen because, among the analyzed frames, it contained a water molecule in a favorable position to attack the covalent bond between fructose C2 and Asp38. In addition, this same water molecule was oriented in a way compatible with subsequent proton donation to Glu245, the residue that acts as a base in this step of the mechanism.
</p>

<p align="justify">
Therefore, the goal is to build a reactive system containing: the Asp38–fructose covalent intermediate, represented as a modified residue named <strong>ASJ</strong>; the reactive water molecule, kept as an independent ligand named <strong>LIG</strong>; and Glu245, included in the quantum region because it participates in proton transfer.
</p>

---

## 2. Why this step was not prepared directly in CHARMM-GUI

<p align="justify">
The first reaction step could be prepared with the aid of CHARMM-GUI because the system still contained sucrose bound noncovalently to the active site. However, after formation of the covalent intermediate between Asp38 and fructose, the direct use of CHARMM-GUI was no longer suitable. The observed problem was the automatic addition of an <strong>ROH</strong> group in the region of the bond between fructose C2 and Asp38. This modification altered the chemical identity of the covalent intermediate and prevented the <code>Asp38@OD2–fructose@C2</code> bond from being represented exactly as intended.
</p>

<p align="justify">
For this reason, the second step was prepared outside CHARMM-GUI. The adopted solution was to transform Asp38 and the covalently bound fructose into a single modified residue named <strong>ASJ</strong>. Thus, the covalent bond between the <code>OD2</code> oxygen of Asp38 and the <code>C2</code> carbon of fructose becomes part of the system topology. The water molecule selected for the attack was kept separate, under the name <strong>LIG</strong>, so that it would not be mixed with regular solvent water molecules.
</p>

---

## 3. General description of the reaction step

The second modeled step involves:

- nucleophilic attack of the `LIG` water molecule on the fructose `C2` carbon;
- cleavage of the covalent bond between `ASJ38@OD2` and `ASJ38@C2`;
- proton transfer from water to `GLU245@OE2`;
- recovery of the catalytic Asp38 residue;
- release of fructose after hydrolysis of the covalent intermediate.

Central atoms of the reaction:

```text
ASJ38@OD2   oxygen of Asp38 covalently bound to fructose
ASJ38@C2    fructose C2 carbon in the covalent intermediate
LIG@O       oxygen of the water molecule selected for the attack
LIG@H1      water hydrogen oriented toward Glu245
GLU245@OE2  oxygen of Glu245 that receives the proton
```

Main distances and angle:

```text
ASJ38@OD2 -- ASJ38@C2
LIG@O     -- ASJ38@C2
LIG@H1    -- GLU245@OE2
ASJ38@OD2 -- ASJ38@C2 -- LIG@O
```

---

## 4. Generic repository organization

<p align="justify">
The repository was organized to contain only the scripts, control files, and templates required to reproduce the preparation. Files that directly expose the molecular system, such as PDBs, MOL2s, topologies, coordinates, trajectories, and execution logs, should not be included on GitHub. The user must provide their own input file in the <code>input/</code> directory.
</p>

Recommended structure:

```text
.
├── input/
│   └── reactive_intermediate_ASP38_fructose_LIG.pdb
├── scripts/
│   ├── merge_ASP38_0CU_to_ASJ_keep_LIG.py
│   ├── validate_ASJ_contiguity.py
│   ├── check_reactive_geometry.py
│   ├── check_reactive_geometry_flexible.py
│   ├── extract_LIG_as_ligand.py
│   ├── extract_ASJ_model_from_complex.py
│   ├── make_tleap_ready_pdb.py
│   ├── patch_ASJ_mol2_DU_types.py
│   ├── fix_mol2_counts.py
│   ├── fix_atoms_per_molecule_with_parmed.py
│   ├── generate_LIG_reactive_restraints.py
│   ├── check_ASJ_internal_geometry.py
│   ├── check_after_mm_qmready.sh
│   ├── diagnose_after_prepare.sh
│   └── diagnose_tleap_errors.sh
├── templates/
│   ├── mdin/
│   │   ├── 01_min_hydrogens.in
│   │   ├── 02_min_restrained_active_site.in
│   │   ├── 03_min_weak_restraints.in
│   │   ├── 04_heat_NVT_323K.in
│   │   ├── 05_eq_NPT_323K.in
│   │   ├── 06a_min_QMMM_DFTB3.in
│   │   ├── 06_test_QMMM_DFTB3_100steps.in
│   │   └── 06_prod_QMMM_DFTB3_ASJ_LIG.in
│   ├── run/
│   │   ├── run_mm_pmemd_cuda.sh
│   │   ├── run_qmmm_minimization.sh
│   │   ├── run_qmmm_test_100steps.sh
│   │   └── run_qmmm_sanderMPI.sh
│   └── tleap/
│       └── tleap_ASJ_LIG_template.in
├── prepare_ASJ_LIG_system.sh
├── README.md
└── .gitignore
```

---

## 5. Expected input file

The user must provide a PDB containing:

```text
1. the complete protein;
2. the covalent intermediate with Asp38 bound to fructose;
3. fructose originally named 0CU;
4. the reactive water molecule renamed as LIG;
5. Glu245 present as GLU;
6. water, ions, and all other system components.
```

Recommended generic name:

```text
input/reactive_intermediate_ASP38_fructose_LIG.pdb
```

The main script accepts the PDB as an argument:

```bash
bash prepare_ASJ_LIG_system.sh input/reactive_intermediate_ASP38_fructose_LIG.pdb
```

If no argument is provided, the script uses the default name internally defined in the `INPUT_PDB` variable.

---

## 6. General execution workflow

```bash
bash prepare_ASJ_LIG_system.sh input/reactive_intermediate_ASP38_fructose_LIG.pdb

bash templates/run/run_mm_pmemd_cuda.sh

bash scripts/check_after_mm_qmready.sh

grep -R "qmcharge" templates/mdin/06*.in
sed -i 's/qmcharge=0/qmcharge=-1/g' templates/mdin/06*.in
sed -i 's/qmcharge = 0/qmcharge=-1/g' templates/mdin/06*.in
grep -R "qmcharge" templates/mdin/06*.in

NP=8 bash templates/run/run_qmmm_minimization.sh

NP=8 bash templates/run/run_qmmm_test_100steps.sh

NP=8 bash templates/run/run_qmmm_sanderMPI.sh
```

---

## 7. Preparation of the ASJ/LIG system

Run:

```bash
bash prepare_ASJ_LIG_system.sh input/reactive_intermediate_ASP38_fructose_LIG.pdb
```

The script performs the following steps:

```text
1. Merges Asp38 and fructose into a single ASJ residue.
2. Keeps the reactive water molecule as LIG.
3. Validates whether ASJ appears as a single contiguous residue.
4. Checks the initial reactive geometry.
5. Extracts LIG as an independent ligand.
6. Parameterizes LIG with antechamber.
7. Generates missing LIG parameters with parmchk2.
8. Extracts a reduced ASJ model.
9. Prepares the PDB for tleap.
10. Parameterizes ASJ with antechamber.
11. Fixes DU types and inappropriate types in ASJ.mol2.
12. Fixes the MOL2 header counts.
13. Generates ASJ.frcmod.
14. Runs tleap.
15. Fixes ATOMS_PER_MOLECULE with ParmEd.
16. Checks ASJ internal integrity.
17. Generates reactive restraints for LIG.
18. Checks the initial ASJ/LIG/GLU245 geometry.
```

---

## 8. Internal preparation commands

### 8.1 Generation of the ASJ residue and preservation of LIG

```bash
python scripts/merge_ASP38_0CU_to_ASJ_keep_LIG.py \
  input/reactive_intermediate_ASP38_fructose_LIG.pdb \
  01_complex_ASJ_LIG.pdb
```

Function:

```text
Transforms Asp38 + 0CU into a single ASJ residue.
Keeps the reactive water molecule as LIG.
```

---

### 8.2 Validation of ASJ contiguity

```bash
python scripts/validate_ASJ_contiguity.py 01_complex_ASJ_LIG.pdb
```

Function:

```text
Confirms whether there is only one ASJ corresponding to residue 38.
Prevents the error of a duplicated ASJ or an ASJ split into two blocks.
```

---

### 8.3 Check of the initial reactive geometry

```bash
python scripts/check_reactive_geometry.py 01_complex_ASJ_LIG.pdb
```

Function:

```text
Measures distances between ASJ, LIG, and GLU245 in the initial PDB.
```

---

### 8.4 Extraction of LIG as an independent ligand

```bash
python scripts/extract_LIG_as_ligand.py \
  01_complex_ASJ_LIG.pdb \
  LIG_only.pdb
```

Function:

```text
Extracts only the LIG molecule for parameterization as a small ligand.
```

---

### 8.5 LIG parameterization

```bash
antechamber -i LIG_only.pdb -fi pdb \
  -o LIG.mol2 -fo mol2 \
  -rn LIG -at gaff2 -c bcc -nc 0 -s 2
```

Explanation of the flags:

```text
-i LIG_only.pdb  input file.
-fi pdb          input file format.
-o LIG.mol2      output file.
-fo mol2         output file format.
-rn LIG          residue name in the MOL2 file.
-at gaff2        uses GAFF2 atom types.
-c bcc           computes AM1-BCC charges.
-nc 0            total charge of the LIG molecule.
-s 2             more detailed message level.
```

Generation of the complementary parameter file:

```bash
parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod
```

Explanation of the flags:

```text
-i LIG.mol2    input MOL2 file.
-f mol2        specifies that the input format is MOL2.
-o LIG.frcmod  missing-parameter file.
```

---

### 8.6 Extraction of the ASJ model

```bash
python scripts/extract_ASJ_model_from_complex.py \
  01_complex_ASJ_LIG.pdb \
  ASJ_model.pdb
```

Function:

```text
Extracts a model containing only the modified ASJ residue for parameterization.
```

---

### 8.7 Preparation of the PDB for tleap

```bash
python scripts/make_tleap_ready_pdb.py \
  01_complex_ASJ_LIG.pdb \
  ASJ_model.pdb \
  01_complex_ASJ_LIG_for_tleap.pdb
```

Function:

```text
Generates a PDB suitable for tleap.
Removes redundant records.
Keeps ASJ as a single residue.
```

---

### 8.8 Initial parameterization of ASJ

```bash
antechamber -i ASJ_model.pdb -fi pdb \
  -o ASJ_raw.mol2 -fo mol2 \
  -rn ASJ -at gaff2 -c bcc -nc 0 -s 2
```

Explanation:

```text
ASJ_model.pdb  model of the modified residue.
ASJ_raw.mol2   initially generated MOL2 file.
-rn ASJ        forces the residue name to ASJ.
-at gaff2      uses GAFF2 as the starting point.
-c bcc         computes AM1-BCC charges.
-nc 0          total charge used for the ASJ model.
```

---

### 8.9 Correction of ASJ atom types

```bash
python scripts/patch_ASJ_mol2_DU_types.py \
  ASJ_raw.mol2 \
  ASJ.mol2
```

Function:

```text
Removes DU types.
Fixes inappropriate types.
Forces essential internal bonds in ASJ.
Adjusts the Asp38 portion of ASJ to protein-compatible types.
```

Main types used in the Asp38 portion:

```text
CA   CT
CB   CT
CG   C
OD1  O
OD2  OS
HA   H1
HB2  H1
HB3  H1
```

`OD2` remains as `OS` because it is covalently bound to fructose C2.

---

### 8.10 Correction of MOL2 counts

```bash
python scripts/fix_mol2_counts.py ASJ.mol2 ASJ.mol2.tmp
mv -f ASJ.mol2.tmp ASJ.mol2
```

Explanation:

```text
fix_mol2_counts.py  fixes the number of atoms and bonds in the MOL2 header.
mv -f               replaces the old file with the corrected file.
```

`mv` flag:

```text
-f  forces replacement of the destination file if it already exists.
```

---

### 8.11 Check for DU types

```bash
grep -n ' DU ' ASJ.mol2
```

Explanation:

```text
grep  searches for patterns in the file.
-n    shows the line number.
' DU ' searches for DU atom types.
```

If any `DU` appears, the preparation must be stopped.

---

### 8.12 Generation of ASJ.frcmod

```bash
parmchk2 -i ASJ.mol2 -f mol2 -o ASJ.frcmod
```

Function:

```text
Generates missing parameters for the modified ASJ residue.
```

---

### 8.13 Running tleap

```bash
tleap -f templates/tleap/tleap_ASJ_LIG_template.in
```

Explanation of the flags:

```text
tleap  AMBER program used to build topology and coordinates.
-f     specifies the tleap command file.
```

`tleap` generates:

```text
system_ASJ_LIG.parm7
system_ASJ_LIG.rst7
system_ASJ_LIG.pdb
leap.log
```

---

### 8.14 Fixing ATOMS_PER_MOLECULE with ParmEd

```bash
python scripts/fix_atoms_per_molecule_with_parmed.py \
  system_ASJ_LIG.parm7 \
  system_ASJ_LIG_fixed.parm7
```

Then:

```bash
mv -f system_ASJ_LIG.parm7 system_ASJ_LIG_raw_from_tleap.parm7
mv -f system_ASJ_LIG_fixed.parm7 system_ASJ_LIG.parm7
```

Function:

```text
Correctly recalculates the ATOMS_PER_MOLECULE block.
Prevents molecular inconsistencies in systems with a modified covalent residue.
```

---

## 9. Problems encountered during parameterization

### 9.1 Unwanted ROH addition

Problem:

```text
The automatic preparation added ROH to the bond between Asp38 and fructose C2.
```

Solution:

```text
Do not use this automatic preparation for this step.
Build the ASJ residue manually.
```

---

### 9.2 Fragmented or duplicated ASJ

Problem:

```text
Asp38 and fructose appeared as separate parts.
The ASJ residue was not contiguous.
```

Solution:

```text
Merge Asp38 and 0CU into a single ASJ38 block.
Validate with validate_ASJ_contiguity.py.
```

---

### 9.3 C-N-c3 and N-c3-C-N parameter errors

Problem observed in `leap.log`:

```text
Could not find angle parameter for atom types: C - N - c3
Could not find angle parameter for atom types: c3 - C - N
No torsion terms for atom types: N-c3-C-N
```

Cause:

```text
The ASJ CA/CB atoms were initially treated as GAFF types incompatible with the peptide environment.
```

Solution:

```text
Use protein-like types in the Asp38 portion of ASJ.
```

---

### 9.4 Internal ASJ disruption during MM

Problem:

```text
After MM, the OD2-C2 bond was correct, but the Asp38 side chain of ASJ separated.
```

Example of error:

```text
CB-CG     > 40 Å
CG-OD1    > 30 Å
CG-OD2    > 40 Å
```

Solution:

```text
Explicitly force the essential internal bonds in ASJ.mol2.
```

Essential bonds:

```text
N-CA
CA-C
C-O
CA-CB
CB-CG
CG-OD1
CG-OD2
OD2-C2
```

---

### 9.5 Escape of the LIG molecule during MM

Problem:

```text
The LIG molecule moved away from fructose C2 during NPT equilibration.
```

Example:

```text
O(LIG)-C2(fructose):     ~20 Å
H1(LIG)-OE2(GLU245):     ~19 Å
```

Solution:

```text
Apply distance and angle restraints during minimization, heating, and equilibration.
```

---

### 9.6 Odd-electron error in QM/MM

Problem:

```text
QMMM: System specified with odd number of electrons (125)
QMMM: but odd spin (1).
```

Cause:

```text
Using qmcharge=0 produced an odd number of electrons for the corrected QM region.
```

Solution:

```text
Use qmcharge=-1.
```

---

## 10. Checking ASJ internal integrity

Run:

```bash
python scripts/check_ASJ_internal_geometry.py system_ASJ_LIG.pdb
```

Expected output:

```text
N-CA       ~1.4–1.5 Å  OK
CA-C       ~1.5–1.6 Å  OK
C-O        ~1.2 Å      OK
CA-CB      ~1.5 Å      OK
CB-CG      ~1.5 Å      OK
CG-OD1     ~1.2 Å      OK
CG-OD2     ~1.3 Å      OK
OD2-C2     ~1.4 Å      OK
```

If any distance is much above 2.2 Å, the system must not be used.

---

## 11. Generation of restraints to keep LIG in the active site

Run:

```bash
python scripts/generate_LIG_reactive_restraints.py \
  system_ASJ_LIG.pdb \
  DISANG_LIG_reactive.RST
```

This script automatically identifies the indices of:

```text
ASJ38@OD2
ASJ38@C2
GLU245@OE2
LIG@O
LIG@H1
```

and generates the file:

```text
DISANG_LIG_reactive.RST
```

---

## 12. Restraints applied to LIG

### 12.1 Distance between fructose C2 and water oxygen

```text
&rst
  iat=ASJ_C2,LIG_O,
  r1=2.4, r2=2.8, r3=3.6, r4=4.3,
  rk2=50.0, rk3=50.0,
/
```

Purpose:

```text
Keep the water oxygen close to fructose C2.
```

---

### 12.2 Distance between water H1 and Glu245 OE2

```text
&rst
  iat=LIG_H1,GLU245_OE2,
  r1=1.3, r2=1.5, r3=2.6, r4=3.2,
  rk2=30.0, rk3=30.0,
/
```

Purpose:

```text
Keep the water hydrogen oriented toward Glu245.
```

---

### 12.3 OD2-C2-O(LIG) attack angle

```text
&rst
  iat=ASJ_OD2,ASJ_C2,LIG_O,
  r1=130.0, r2=150.0, r3=180.0, r4=181.0,
  rk2=10.0, rk3=10.0,
/
```

Purpose:

```text
Maintain a geometry compatible with nucleophilic attack on fructose C2.
```

---

## 13. Meaning of the restraint parameters

```text
iat   indices of the atoms involved in the restraint.
r1    extreme lower limit.
r2    beginning of the no-penalty region.
r3    end of the no-penalty region.
r4    extreme upper limit.
rk2   force constant applied below r2.
rk3   force constant applied above r3.
```

The region between `r2` and `r3` is flat-bottomed, meaning that no penalty is applied if the distance or angle remains within this interval.

---

## 14. Contents of the input files in templates/mdin/

### 14.1 `01_min_hydrogens.in`

```text
Initial minimization: hydrogens and solvent only, heavy solute restrained
&cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  ntb=1, cut=10.0,
  ntpr=100,
  restraint_wt=25.0,
  restraintmask='!:WAT,Na+,Cl- & !@H=',
  nmropt=1,
/

&wt
  type='END',
/
DISANG=DISANG_LIG_reactive.RST
LISTOUT=POUT_LIG_reactive_01_min_hydrogens
```

---

### 14.2 `02_min_restrained_active_site.in`

```text
Minimization with active site strongly restrained
&cntrl
  imin=1, maxcyc=10000, ncyc=5000,
  ntb=1, cut=10.0,
  ntpr=100,
  restraint_wt=10.0,
  restraintmask='(:38 | :245 | :LIG)',
  nmropt=1,
/

&wt
  type='END',
/
DISANG=DISANG_LIG_reactive.RST
LISTOUT=POUT_LIG_reactive_02_min_restrained_active_site
```

---

### 14.3 `03_min_weak_restraints.in`

```text
Minimization with weak restraints on reactive region
&cntrl
  imin=1, maxcyc=10000, ncyc=5000,
  ntb=1, cut=10.0,
  ntpr=100,
  restraint_wt=2.0,
  restraintmask='(:38 | :245 | :LIG)',
  nmropt=1,
/

&wt
  type='END',
/
DISANG=DISANG_LIG_reactive.RST
LISTOUT=POUT_LIG_reactive_03_min_weak_restraints
```

---

### 14.4 `04_heat_NVT_323K.in`

```text
Heating NVT to 323.15 K with active site restrained
&cntrl
  imin=0, irest=0, ntx=1,
  nstlim=100000, dt=0.001,
  ntc=2, ntf=2,
  cut=10.0, ntb=1,
  ntt=3, gamma_ln=2.0,
  tempi=0.0, temp0=323.15,
  ntpr=500, ntwx=500, ntwr=5000,
  ioutfm=1,
  restraint_wt=5.0,
  restraintmask='(:38 | :245 | :LIG)',
  nmropt=1,
/

&wt
  type='END',
/
DISANG=DISANG_LIG_reactive.RST
LISTOUT=POUT_LIG_reactive_04_heat_NVT_323K
```

---

### 14.5 `05_eq_NPT_323K.in`

```text
Equilibration NPT at 323.15 K with weaker active-site restraints
&cntrl
  imin=0, irest=1, ntx=5,
  nstlim=250000, dt=0.001,
  ntc=2, ntf=2,
  cut=10.0, ntb=2, ntp=1, taup=2.0,
  ntt=3, gamma_ln=2.0,
  temp0=323.15,
  ntpr=500, ntwx=500, ntwr=5000,
  ioutfm=1,
  restraint_wt=1.0,
  restraintmask='(:38 | :245 | :LIG)',
  nmropt=1,
/

&wt
  type='END',
/
DISANG=DISANG_LIG_reactive.RST
LISTOUT=POUT_LIG_reactive_05_eq_NPT_323K
```

---

### 14.6 `06a_min_QMMM_DFTB3.in`

```text
QM/MM minimization DFTB3 before production
&cntrl
  imin=1,
  maxcyc=5000,
  ncyc=2500,
  ntmin=1,
  cut=10.0,
  ntb=1,
  ntpr=50,
  ifqnt=1,
/
&qmmm
  qmmask='(:38@CB,HB2,HB3,CG,OD1,OD2,C1,O1,C2,C3,O3,C4,O4,C5,O5,C6,O6,H3,H4,H5,H11,H12,H1O,H3O,H61,H62,H4O,H6O | :245@CB,HB2,HB3,CG,HG2,HG3,CD,OE1,OE2 | :LIG)',
  qmcharge=-1,
  qm_theory='DFTB3',
  qmcut=12.0,
  qm_ewald=1,
  qmshake=0,
  scfconv=1.0d-6,
  itrmax=2000,
/
```

---

### 14.7 `06_test_QMMM_DFTB3_100steps.in`

```text
Short test QM/MM-DFTB3 with ASJ38 + GLU245 + LIG in QM region
&cntrl
  imin=0, irest=1, ntx=5,
  nstlim=100, dt=0.0005,
  ntc=1, ntf=1,
  cut=10.0,
  ntb=2, ntp=1, taup=2.0,
  ntt=3, gamma_ln=2.0,
  temp0=323.15,
  ntpr=1, ntwx=10, ntwr=100,
  ioutfm=1,
  ifqnt=1,
/
&qmmm
  qmmask='(:38@CB,HB2,HB3,CG,OD1,OD2,C1,O1,C2,C3,O3,C4,O4,C5,O5,C6,O6,H3,H4,H5,H11,H12,H1O,H3O,H61,H62,H4O,H6O | :245@CB,HB2,HB3,CG,HG2,HG3,CD,OE1,OE2 | :LIG)',
  qmcharge=-1,
  qm_theory='DFTB3',
  qmcut=12.0,
  writepdb=1,
  qm_ewald=1,
  qmshake=0,
/
```

---

### 14.8 `06_prod_QMMM_DFTB3_ASJ_LIG.in`

```text
Production QM/MM-DFTB3 with ASJ38 + GLU245 + LIG in QM region
&cntrl
  imin=0, irest=1, ntx=5,
  nstlim=2000000, dt=0.0005,
  ntc=1, ntf=1,
  cut=10.0,
  ntb=2, ntp=1, taup=2.0,
  ntt=3, gamma_ln=2.0,
  temp0=323.15,
  ntpr=100, ntwx=100, ntwr=5000,
  ioutfm=1,
  ifqnt=1,
/
&qmmm
  qmmask='(:38@CB,HB2,HB3,CG,OD1,OD2,C1,O1,C2,C3,O3,C4,O4,C5,O5,C6,O6,H3,H4,H5,H11,H12,H1O,H3O,H61,H62,H4O,H6O | :245@CB,HB2,HB3,CG,HG2,HG3,CD,OE1,OE2 | :LIG)',
  qmcharge=-1,
  qm_theory='DFTB3',
  qmcut=12.0,
  writepdb=1,
  qm_ewald=1,
  qmshake=0,
/
```

---

## 15. Explanation of the main parameters in the `.in` files

### Classical minimization

```text
imin=1      runs energy minimization.
maxcyc      total number of minimization cycles.
ncyc        number of cycles before switching the minimization algorithm.
ntb=1       uses a periodic box at constant volume.
cut=10.0    nonbonded cutoff in Å.
ntpr        energy printing frequency.
```

### Classical dynamics

```text
imin=0      runs molecular dynamics.
irest=0     starts a new dynamics run.
irest=1     restarts dynamics from saved velocities.
ntx=1       reads coordinates only.
ntx=5       reads coordinates and velocities.
nstlim      number of steps.
dt          integration timestep in ps.
ntc=2       applies SHAKE to bonds involving hydrogen.
ntf=2       removes forces for these constrained bonds.
ntt=3       Langevin thermostat.
gamma_ln    Langevin collision frequency.
temp0       target temperature.
ntb=2       constant-pressure dynamics.
ntp=1       isotropic pressure coupling.
taup        pressure relaxation time.
ntwx        trajectory writing frequency.
ntwr        restart writing frequency.
ioutfm=1    NetCDF trajectory output.
```

### Classical restraints

```text
restraint_wt      positional restraint force.
restraintmask     mask of restrained atoms.
nmropt=1          activates AMBER/NMR-style restraints.
DISANG            file containing distance/angle restraints.
LISTOUT           restraint output file.
```

### QM/MM

```text
ifqnt=1             activates QM/MM.
qmmask              defines the atoms in the QM region.
qmcharge=-1         total charge of the QM region.
qm_theory='DFTB3'   DFTB3 quantum method.
qmcut=12.0          QM/MM cutoff.
qm_ewald=1          uses Ewald electrostatic treatment.
qmshake=0           does not apply SHAKE to QM atoms.
scfconv             SCF convergence criterion.
itrmax              maximum number of SCF iterations.
writepdb=1          writes a PDB of the QM region at the beginning of the simulation.
```

---

## 16. Running classical MM

Run:

```bash
bash templates/run/run_mm_pmemd_cuda.sh
```

The script runs:

```bash
pmemd.cuda -O -i templates/mdin/01_min_hydrogens.in \
  -p system_ASJ_LIG.parm7 -c system_ASJ_LIG.rst7 \
  -o 01_min_hydrogens.out -r 01_min_hydrogens.rst7 \
  -inf 01_min_hydrogens.mdinfo -ref system_ASJ_LIG.rst7

pmemd.cuda -O -i templates/mdin/02_min_restrained_active_site.in \
  -p system_ASJ_LIG.parm7 -c 01_min_hydrogens.rst7 \
  -o 02_min_restrained_active_site.out -r 02_min_restrained_active_site.rst7 \
  -inf 02_min_restrained_active_site.mdinfo -ref 01_min_hydrogens.rst7

pmemd.cuda -O -i templates/mdin/03_min_weak_restraints.in \
  -p system_ASJ_LIG.parm7 -c 02_min_restrained_active_site.rst7 \
  -o 03_min_weak_restraints.out -r 03_min_weak_restraints.rst7 \
  -inf 03_min_weak_restraints.mdinfo -ref 02_min_restrained_active_site.rst7

pmemd.cuda -O -i templates/mdin/04_heat_NVT_323K.in \
  -p system_ASJ_LIG.parm7 -c 03_min_weak_restraints.rst7 \
  -o 04_heat_NVT_323K.out -r 04_heat_NVT_323K.rst7 \
  -x 04_heat_NVT_323K.nc -inf 04_heat_NVT_323K.mdinfo \
  -ref 03_min_weak_restraints.rst7

pmemd.cuda -O -i templates/mdin/05_eq_NPT_323K.in \
  -p system_ASJ_LIG.parm7 -c 04_heat_NVT_323K.rst7 \
  -o 05_eq_NPT_323K.out -r 05_eq_NPT_323K.rst7 \
  -x 05_eq_NPT_323K.nc -inf 05_eq_NPT_323K.mdinfo \
  -ref 04_heat_NVT_323K.rst7
```

`pmemd.cuda` flags:

```text
-O      overwrites existing output files.
-i      input file.
-p      AMBER topology.
-c      input coordinates/restart.
-o      main output file.
-r      output restart.
-inf    dynamics information file.
-x      NetCDF trajectory.
-ref    reference coordinates for positional restraints.
```

---

## 17. Post-MM check

Run:

```bash
bash scripts/check_after_mm_qmready.sh
```

This script creates:

```text
05_eq_NPT_323K_qmready.pdb
05_eq_NPT_323K_qmready.rst7
```

Internal commands:

```bash
cat > check_final_mm.cpptraj << 'EOF'
parm system_ASJ_LIG.parm7
trajin 05_eq_NPT_323K.rst7
autoimage anchor :38
trajout 05_eq_NPT_323K_qmready.pdb pdb
trajout 05_eq_NPT_323K_qmready.rst7 restart
run
EOF

cpptraj -i check_final_mm.cpptraj
```

Explanation:

```text
parm                 loads the topology.
trajin               reads the final NPT restart.
autoimage anchor :38 recenters the system using ASJ as the anchor.
trajout ... pdb      writes a PDB for inspection.
trajout ... restart  writes a restart for QM/MM.
run                  executes the commands.
```

Then the script runs:

```bash
python scripts/check_ASJ_internal_geometry.py 05_eq_NPT_323K_qmready.pdb
python scripts/check_reactive_geometry_flexible.py 05_eq_NPT_323K_qmready.pdb
```

Expected geometry:

```text
OD2(ASJ38)-C2(fructose): ~1.4 Å
O(LIG)-C2(fructose):     ~2.8–3.6 Å
H1(LIG)-OE2(GLU245):     ~1.5–2.6 Å
OD2-C2-O(LIG):           ~150–180°
```

---

## 18. Mandatory correction of the QM charge

Before QM/MM minimization, check:

```bash
grep -R "qmcharge" templates/mdin/06*.in
```

If any file still contains `qmcharge=0`, correct it:

```bash
sed -i 's/qmcharge=0/qmcharge=-1/g' templates/mdin/06*.in
sed -i 's/qmcharge = 0/qmcharge=-1/g' templates/mdin/06*.in
```

Check again:

```bash
grep -R "qmcharge" templates/mdin/06*.in
```

Expected result:

```text
qmcharge=-1,
```

Explanation of `sed` flags:

```text
-i      edits the file in place.
s/a/b/g replaces all occurrences of a with b.
```

---

## 19. QM/MM minimization

Run:

```bash
NP=8 bash templates/run/run_qmmm_minimization.sh
```

Executed command:

```bash
mpirun -np ${NP:-8} sander.MPI -O \
  -i templates/mdin/06a_min_QMMM_DFTB3.in \
  -p system_ASJ_LIG.parm7 \
  -c 05_eq_NPT_323K_qmready.rst7 \
  -o 06a_min_QMMM_DFTB3.out \
  -r 06a_min_QMMM_DFTB3.rst7 \
  -inf 06a_min_QMMM_DFTB3.mdinfo
```

Flags:

```text
NP=8        defines 8 MPI processes.
mpirun      starts parallel execution.
-np         number of processes.
sander.MPI  AMBER engine used for QM/MM.
-O          overwrites old output files.
-i          QM/MM input file.
-p          topology.
-c          input restart.
-o          main output.
-r          output restart.
-inf        execution information file.
```

---

## 20. Short QM/MM test

Before the long production run, execute 100 steps:

```bash
NP=8 bash templates/run/run_qmmm_test_100steps.sh
```

Executed command:

```bash
mpirun -np ${NP:-8} sander.MPI -O \
  -i templates/mdin/06_test_QMMM_DFTB3_100steps.in \
  -p system_ASJ_LIG.parm7 \
  -c 06a_min_QMMM_DFTB3.rst7 \
  -o 06_test_QMMM_DFTB3_100steps.out \
  -r 06_test_QMMM_DFTB3_100steps.rst7 \
  -x 06_test_QMMM_DFTB3_100steps.nc \
  -inf 06_test_QMMM_DFTB3_100steps.mdinfo
```

If the script is configured to start from `05_eq_NPT_323K.rst7`, change the `CRD` variable to:

```bash
CRD=06a_min_QMMM_DFTB3.rst7
```

---

## 21. Checking QM/MM errors

After minimization or the short test:

```bash
grep -c "odd number of electrons" 06a_min_QMMM_DFTB3.out
grep -c "Convergence could not be achieved" 06a_min_QMMM_DFTB3.out
grep "NSTEP" -A8 06a_min_QMMM_DFTB3.out | tail -30
```

For the short test:

```bash
grep -c "odd number of electrons" 06_test_QMMM_DFTB3_100steps.out
grep -c "Convergence could not be achieved" 06_test_QMMM_DFTB3_100steps.out
grep "NSTEP" -A8 06_test_QMMM_DFTB3_100steps.out | tail -30
```

Expected results:

```text
0
0
```

or only a few convergence failures at the beginning of minimization.

---

## 22. QM/MM production

Run:

```bash
NP=8 bash templates/run/run_qmmm_sanderMPI.sh
```

Recommended command:

```bash
mpirun -np 8 sander.MPI -O \
  -i templates/mdin/06_prod_QMMM_DFTB3_ASJ_LIG.in \
  -p system_ASJ_LIG.parm7 \
  -c 06a_min_QMMM_DFTB3.rst7 \
  -o 06_prod_QMMM_DFTB3_ASJ_LIG.out \
  -r 06_prod_QMMM_DFTB3_ASJ_LIG.rst7 \
  -x 06_prod_QMMM_DFTB3_ASJ_LIG.nc \
  -inf 06_prod_QMMM_DFTB3_ASJ_LIG.mdinfo
```

If the script is configured to start from `05_eq_NPT_323K.rst7`, change the `CRD` variable to:

```bash
CRD=06a_min_QMMM_DFTB3.rst7
```

---

## 23. Production monitoring

During production:

```bash
tail -f 06_prod_QMMM_DFTB3_ASJ_LIG.out
```

or:

```bash
tail -f 06_prod_QMMM_DFTB3_ASJ_LIG.mdinfo
```

Check errors:

```bash
grep -c "odd number of electrons" 06_prod_QMMM_DFTB3_ASJ_LIG.out
grep -c "Convergence could not be achieved" 06_prod_QMMM_DFTB3_ASJ_LIG.out
```

Check progress:

```bash
grep "NSTEP" -A8 06_prod_QMMM_DFTB3_ASJ_LIG.out | tail -30
```

---

## 24. Expected files

After preparation:

```text
system_ASJ_LIG.parm7
system_ASJ_LIG.rst7
system_ASJ_LIG.pdb
system_ASJ_LIG_raw_from_tleap.parm7
ASJ.mol2
ASJ.frcmod
LIG.mol2
LIG.frcmod
DISANG_LIG_reactive.RST
```

After MM:

```text
01_min_hydrogens.rst7
02_min_restrained_active_site.rst7
03_min_weak_restraints.rst7
04_heat_NVT_323K.rst7
04_heat_NVT_323K.nc
05_eq_NPT_323K.rst7
05_eq_NPT_323K.nc
05_eq_NPT_323K_qmready.pdb
05_eq_NPT_323K_qmready.rst7
```

After QM/MM:

```text
06a_min_QMMM_DFTB3.out
06a_min_QMMM_DFTB3.rst7
06_test_QMMM_DFTB3_100steps.out
06_test_QMMM_DFTB3_100steps.rst7
06_prod_QMMM_DFTB3_ASJ_LIG.out
06_prod_QMMM_DFTB3_ASJ_LIG.rst7
06_prod_QMMM_DFTB3_ASJ_LIG.nc
```

---

## 25. Suggested analyses for the second step

File:

```text
mechanism_distances_step2.in
```

Content:

```cpptraj
parm system_ASJ_LIG.parm7
trajin 06_prod_QMMM_DFTB3_ASJ_LIG.nc

distance d_ASJ_OD2_C2 :38@OD2 :38@C2 out d_ASJ_OD2_C2.dat
distance d_LIG_O_C2 :LIG@O :38@C2 out d_LIG_O_C2.dat
distance d_LIG_H1_GLU245_OE2 :LIG@H1 :245@OE2 out d_LIG_H1_GLU245_OE2.dat
angle a_OD2_C2_LIG_O :38@OD2 :38@C2 :LIG@O out a_OD2_C2_LIG_O.dat

run
quit
```

Run:

```bash
cpptraj -i mechanism_distances_step2.in
```

Interpretation:

```text
ASJ@OD2 -- ASJ@C2 increases if the covalent bond is broken.
LIG@O -- ASJ@C2 decreases if water attacks C2.
LIG@H1 -- GLU245@OE2 decreases during proton transfer.
OD2-C2-O(LIG) indicates the geometry of the nucleophilic attack.
```

---

## 26. Simple script to plot distances

File:

```text
plot_step2_mechanism.py
```

Content:

```python
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

FILES = {
    "d_ASJ_OD2_C2.dat": "ASJ OD2 - ASJ C2",
    "d_LIG_O_C2.dat": "LIG O - ASJ C2",
    "d_LIG_H1_GLU245_OE2.dat": "LIG H1 - GLU245 OE2",
}

def read_dat(path):
    x, y = [], []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split()
        x.append(float(f[0]))
        y.append(float(f[1]))
    return np.array(x), np.array(y)

plt.figure(figsize=(10, 6))

for file, label in FILES.items():
    if Path(file).exists():
        x, y = read_dat(file)
        plt.plot(x, y, label=label)

plt.xlabel("Frame")
plt.ylabel("Distance (Å)")
plt.legend()
plt.tight_layout()
plt.savefig("step2_mechanism_distances.png", dpi=600)
plt.show()
```

Run:

```bash
python plot_step2_mechanism.py
```

---

## 27. Summary of the final workflow

```bash
# 1. Prepare the ASJ/LIG system
bash prepare_ASJ_LIG_system.sh input/reactive_intermediate_ASP38_fructose_LIG.pdb

# 2. Run MM with LIG restraints
bash templates/run/run_mm_pmemd_cuda.sh

# 3. Generate a QM/MM-ready restart and check geometry
bash scripts/check_after_mm_qmready.sh

# 4. Ensure the correct charge of the QM region
grep -R "qmcharge" templates/mdin/06*.in
sed -i 's/qmcharge=0/qmcharge=-1/g' templates/mdin/06*.in
sed -i 's/qmcharge = 0/qmcharge=-1/g' templates/mdin/06*.in
grep -R "qmcharge" templates/mdin/06*.in

# 5. QM/MM minimization
NP=8 bash templates/run/run_qmmm_minimization.sh

# 6. Short QM/MM test
NP=8 bash templates/run/run_qmmm_test_100steps.sh

# 7. QM/MM production
NP=8 bash templates/run/run_qmmm_sanderMPI.sh
```

---

## 28. Final remarks

<p align="justify">
The central point of this preparation was to correctly represent the covalent intermediate between Asp38 and fructose while avoiding the undesired modification automatically generated during graphical-interface preparation. The creation of the ASJ residue allowed the <code>ASJ38@OD2–ASJ38@C2</code> covalent bond to be explicitly maintained in the topology. The water molecule selected from frame 14 of the first-step dynamics was preserved as LIG and kept in the active site by geometric restraints during the classical stage. The final correction of the QM-region charge to <code>qmcharge=-1</code> removed the odd-electron error and allowed the QM/MM-DFTB3 simulation of the second reaction step to begin.
</p>

## References

<p align="justify">
JITONNOM, J.; KETUDAT-CAIRNS, J. R.; HANNONGBUA, S. QM/MM Modeling of the Hydrolysis and Transfructosylation Reactions of Fructosyltransferase from Aspergillus Japonicas, an Enzyme That Produces Prebiotic Fructooligosaccharide. Journal of Molecular Graphics and Modelling, v. 79, p. 175–184, Jan. 2018.
</p>

<p align="justify">
CASE, D. A.; AKTULGA, H. M.; BELFON, K.; CERUTTI, D. S.; CISNEROS, G. A.; CRUZEIRO, V. W. D.; FOROUZESH, N.; GIESE, T. J.; GÖTZ, A. W.; GOHLKE, H.; IZADI, S.; KASAVAJHALA, K.; KAYMAK, M. C.; KING, E.; KURTZMAN, T.; LEE, T.-S.; LI, P.; LIU, J.; LUCHKO, T.; LUO, R.; MANATHUNGA, M.; MACHADO, M. R.; NGUYEN, H. M.; O’HEARN, K. A.; ONUFRIEV, A. V.; PAN, F.; PANTANO, S.; QI, R.; RAHNAMOUN, A.; RISHEH, A.; SCHOTT-VERDUGO, S.; SHAJAN, A.; SWAILS, J.; WANG, J.; WEI, H.; WU, X.; WU, Y.; ZHANG, S.; ZHAO, S.; ZHU, Q.; CHEATHAM, T. E.; ROE, D. R.; ROITBERG, A.; SIMMERLING, C.; YORK, D. M.; NAGAN, M. C.; MERZ, K. M. AmberTools. Journal of Chemical Information and Modeling, v. 63, n. 20, p. 6183–6191, 23 Oct. 2023. 
</p>

<p align="justify">
CASE, D. A.; CHEATHAM, T. E.; DARDEN, T.; GOHLKE, H.; LUO, R.; MERZ, K. M.; ONUFRIEV, A.; SIMMERLING, C.; WANG, B.; WOODS, R. J. The Amber Biomolecular Simulation Programs. Journal of Computational Chemistry, v. 26, n. 16, p. 1668–1688, Dec. 2005. 
</p>
