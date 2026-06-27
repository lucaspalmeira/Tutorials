# DFTB3-Based QM/MM MD in AMBER: A Tutorial on Sucrose Hydrolysis by Invertase

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

trajout cluster2_rep_frame10.pdb pdb onlyframes 10
trajout cluster2_rep_frame11.pdb pdb onlyframes 11

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

## 12. Preparation of the covalent fructose intermediate for the second reaction step

After the first QM/MM simulation, a representative structure was selected in which the fructose unit remained covalently bound to the nucleophilic residue Asp38. In this second stage, the glucose molecule was removed, leaving only the covalent fructose intermediate. A catalytic water molecule was then included near the active site to promote hydrolysis of the covalent Asp38–fructose intermediate.

The second reaction step is therefore defined as:

```text
Asp38@OD2 — 0CD@C2  +  H2O  →  Asp38@OD2 + fructose-OH
```

The catalytic water is expected to attack the C2 atom of fructose, while one of its protons is transferred to the general acid/base residue GLU/GLH245.

For this system, the relevant residues are:

```text
0CD/fructose covalent intermediate = residue 601
Catalytic water                    = residue 757
Nucleophile                         = ASP38
Transition-state stabilizer         = ASP168
General acid/base                   = GLU/GLH245
```

The catalytic water was selected based on geometric criteria:

```text
O_water — C2_0CD distance
O_water — C2_0CD — OD2_ASP38 attack angle
O_water/H_water proximity to OE1/OE2 of GLU/GLH245
```

Before minimization, the catalytic water 757 showed:

```text
O_water757 — C2_0CD      = 5.5745 Å
O_water757 — C2 — OD2    = 140.0239 degrees
O_water757 — GLU245 OE1  = 2.7006 Å
O_water757 — GLU245 OE2  = 4.1050 Å
```

This geometry indicates that the water is not yet close enough for direct nucleophilic attack, but it is well oriented relative to the Asp38–fructose covalent bond and well positioned relative to GLU/GLH245. Therefore, it was retained as the catalytic water candidate and restrained during minimization and equilibration.

---

## 13. Correction of the covalent intermediate PDB

The intermediate PDB must preserve the covalent bond:

```text
ASP38@OD2 — 0CD@C2
```

The `CONECT` records should not define a bond between fructose and an artificial `ROH` group. The corrected PDB should contain only the intended covalent connection between Asp38 and fructose, plus the internal connectivity of the fructose residue.

In PyMOL, the covalent bond can be checked with:

```python
load invertase_with_water_best_frame16.pdb

select c2_fru, resn 0CD and resi 601 and name C2
select asp38_od2, resi 38 and name OD2
select wat757, resi 757
select wat757_O, resi 757 and name O
select glu245_oe, resi 245 and name OE1+OE2

distance d_Asp38_C2, asp38_od2, c2_fru
distance d_WAT757_C2, wat757_O, c2_fru
distance d_WAT757_GLU245, wat757_O, glu245_oe
angle angle_attack, wat757_O, c2_fru, asp38_od2

show sticks, resn 0CD or resi 38+168+245+757
zoom resn 0CD or resi 38+168+245+757, 8
```

Expected covalent distance:

```text
ASP38@OD2 — 0CD@C2 ≈ 1.4–1.7 Å
```

---

## 14. Problem introduced by CHARMM-GUI: artificial ROH group

During the preparation of the covalent ligand system, CHARMM-GUI may introduce an artificial residue named `ROH` near the covalent bond. In this system, the `ROH` group was placed between the fructose C2 atom and the Asp38 OD2 atom.

This is chemically incorrect for the intended QM/MM reaction, because the desired intermediate is:

```text
ASP38@OD2 — 0CD@C2
```

not:

```text
0CD@C2 — ROH@O1
```

After CHARMM-GUI preparation, the topology may contain:

```text
residue 601 = ROH
residue 602 = 0CD
```

After removing `ROH`, the fructose residue becomes:

```text
residue 601 = 0CD
```

In the corrected system used here:

```text
0CD/fructose = residue 601
atoms 9029–9050
```

---

## 15. Removing the ROH residue from the AMBER topology and coordinates

The `ROH` residue was removed using `cpptraj`. This avoids the need for `ParmEd`, which may fail depending on the AmberTools installation.

Create the script:

```bash
nano remove_roh_cpptraj.sh
```

Paste:

```bash
#!/bin/bash

# Remove ROH from the topology
cpptraj << EOF
parm step3_input.parm7
parmstrip :ROH
parmwrite out step3_input_noROH.parm7
quit
EOF

# Remove ROH from the restart coordinates
cpptraj << EOF
parm step3_input.parm7
trajin step3_input.rst7
strip :ROH
trajout step3_input_noROH.rst7 restart
run
quit
EOF

# Generate a PDB without ROH for visual inspection
cpptraj << EOF
parm step3_input.parm7
trajin step3_input.rst7
strip :ROH
trajout step3_input_noROH.pdb pdb
run
quit
EOF
```

Run:

```bash
chmod +x remove_roh_cpptraj.sh
./remove_roh_cpptraj.sh
```

Expected output files:

```text
step3_input_noROH.parm7
step3_input_noROH.rst7
step3_input_noROH.pdb
```

Check that `ROH` was removed:

```bash
grep ROH step3_input_noROH.pdb
```

If nothing is printed, the `ROH` group was removed successfully.

The corrected residue numbering should now be:

```text
601 0CD
757 WAT
```

---

## 16. Renaming the catalytic water from WAT to CWT in the topology

After removing `ROH`, the catalytic water was residue 757. However, this water contains only three atoms:

```text
757 WAT = O H1 H2
```

while the bulk waters generated by CHARMM-GUI contain an additional extra point:

```text
WAT = O H1 H2 EPW
```

Therefore, AMBER may fail with the error:

```text
Error: Fast 3-point water residue, name and bond data incorrect!
```

This happens because the catalytic water and the bulk water residues have the same residue name `WAT`, but different atom/bond definitions.

To avoid this, the catalytic water residue 757 was renamed from `WAT` to `CWT` in the topology. This separates the catalytic water from the bulk water model.

Create the script:

```bash
nano rename_wat757_to_cwt.py
```

Paste:

```python
#!/usr/bin/env python3

from pathlib import Path

input_parm = Path("step3_input_noROH.parm7")
output_parm = Path("step3_input_noROH_CWT757.parm7")

target_resid = 757
new_label = "CWT"

lines = input_parm.read_text().splitlines()

flag_idx = None
for i, line in enumerate(lines):
    if line.startswith("%FLAG RESIDUE_LABEL"):
        flag_idx = i
        break

if flag_idx is None:
    raise RuntimeError("Could not find %FLAG RESIDUE_LABEL in parm7.")

data_start = flag_idx + 2

data_end = data_start
while data_end < len(lines) and not lines[data_end].startswith("%FLAG"):
    data_end += 1

raw = "".join(lines[data_start:data_end])
labels = [raw[i:i+4] for i in range(0, len(raw), 4)]

if target_resid < 1 or target_resid > len(labels):
    raise RuntimeError(f"Residue {target_resid} is outside the valid range. Total residues: {len(labels)}")

old_label = labels[target_resid - 1].strip()
print(f"Residue {target_resid}: {old_label} -> {new_label}")

labels[target_resid - 1] = f"{new_label:<4}"[:4]

new_label_lines = []
for i in range(0, len(labels), 20):
    new_label_lines.append("".join(labels[i:i+20]))

new_lines = lines[:data_start] + new_label_lines + lines[data_end:]

output_parm.write_text("\n".join(new_lines) + "\n")

print(f"Written file: {output_parm}")
```

Run:

```bash
python3 rename_wat757_to_cwt.py
```

This generates:

```text
step3_input_noROH_CWT757.parm7
```

Check the result:

```bash
grep -n "CWT" step3_input_noROH_CWT757.parm7
```

Also check with `cpptraj`:

```bash
cpptraj -p step3_input_noROH_CWT757.parm7 << EOF
resinfo :757
resinfo :758
quit
EOF
```

The expected result is:

```text
757 CWT
758 WAT
```

From this point forward, the corrected topology for the second reaction step is:

```text
step3_input_noROH_CWT757.parm7
```

and the corrected coordinates are:

```text
step3_input_noROH.rst7
```

---

## 17. Checking catalytic water geometry before minimization

Before minimization, the catalytic water geometry was checked with `cpptraj`.

Create:

```bash
nano check_wat757_before_min.in
```

Paste:

```cpptraj
parm step3_input_noROH_CWT757.parm7

# Initial structure before minimization
trajin step3_input_noROH.rst7

# Distance between catalytic water oxygen and fructose C2
distance d_wat757_C2 :757@O :601@C2 out before_min_dist_wat757_C2.dat

# Distance of the covalent Asp38-fructose bond
distance d_Asp38_C2 :38@OD2 :601@C2 out before_min_dist_Asp38_C2.dat

# Nucleophilic attack angle: O_water — C2_fructose — OD2_Asp38
angle ang_attack :757@O :601@C2 :38@OD2 out before_min_angle_wat757_C2_Asp38.dat

# Distance between catalytic water oxygen and GLU/GLH245
distance d_Owat_GLU245_OE1 :757@O :245@OE1 out before_min_dist_Owat_GLU245_OE1.dat
distance d_Owat_GLU245_OE2 :757@O :245@OE2 out before_min_dist_Owat_GLU245_OE2.dat

# Distances between water hydrogens and GLU/GLH245
distance d_H1_GLU245_OE1 :757@H1 :245@OE1 out before_min_dist_H1_GLU245_OE1.dat
distance d_H1_GLU245_OE2 :757@H1 :245@OE2 out before_min_dist_H1_GLU245_OE2.dat
distance d_H2_GLU245_OE1 :757@H2 :245@OE1 out before_min_dist_H2_GLU245_OE1.dat
distance d_H2_GLU245_OE2 :757@H2 :245@OE2 out before_min_dist_H2_GLU245_OE2.dat

run
quit
```

Run:

```bash
cpptraj -i check_wat757_before_min.in
```

Inspect the results:

```bash
cat before_min_dist_wat757_C2.dat
cat before_min_angle_wat757_C2_Asp38.dat
cat before_min_dist_Owat_GLU245_OE1.dat
cat before_min_dist_Owat_GLU245_OE2.dat
cat before_min_dist_H1_GLU245_OE1.dat
cat before_min_dist_H1_GLU245_OE2.dat
cat before_min_dist_H2_GLU245_OE1.dat
cat before_min_dist_H2_GLU245_OE2.dat
```

The observed values before minimization were:

```text
O_water757 — C2_0CD      = 5.5745 Å
O_water757 — C2 — OD2    = 140.0239 degrees
O_water757 — GLU245 OE1  = 2.7006 Å
O_water757 — GLU245 OE2  = 4.1050 Å
```

The water was therefore retained because it had a good attack angle and good proximity to GLU/GLH245.

---

## 18. Restraining the covalent fructose and catalytic water during minimization

The QM region is not used during minimization. The minimization is classical, but the covalent fructose intermediate and the catalytic water are restrained to preserve the reactive geometry.

The restrained residues are:

```text
0CD/fructose = residue 601
CWT/water    = residue 757
```

The relevant section in `step4.0_minimization.mdin` is:

```text
ntr=1,
nmropt=1,
```

At the end of `step4.0_minimization.mdin`, use:

```text
Ligand 0CD and catalytic water 757 positional restraints
10.0
RES 601 601
RES 757 757
END
END
```

A complete example of the final section is:

```text
 &wt
    type='END'
 /
LISTIN=POUT
LISTOUT=POUT
&end
Ligand 0CD and catalytic water 757 positional restraints
10.0
RES 601 601
RES 757 757
END
END
```

Run minimization with the corrected topology:

```bash
pmemd.cuda -O \
  -i step4.0_minimization.mdin \
  -p step3_input_noROH_CWT757.parm7 \
  -c step3_input_noROH.rst7 \
  -o step4.0_minimization.mdout \
  -r step4.0_minimization.rst7 \
  -inf step4.0_minimization.mdinfo \
  -ref step3_input_noROH.rst7
```

The message below may appear:

```text
Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
```

This is commonly observed in AMBER calculations and is not necessarily a fatal error if the calculation completes normally and the output/restart files are generated.

---

## 19. Checking water geometry after minimization

After minimization, measure the catalytic geometry again:

```bash
cpptraj -p step3_input_noROH_CWT757.parm7 << EOF
trajin step4.0_minimization.rst7

distance d_wat757_C2 :757@O :601@C2 out after_min_dist_wat757_C2.dat
distance d_Asp38_C2 :38@OD2 :601@C2 out after_min_dist_Asp38_C2.dat
angle ang_attack :757@O :601@C2 :38@OD2 out after_min_angle_wat757_C2_Asp38.dat
distance d_Owat_GLU245_OE1 :757@O :245@OE1 out after_min_dist_Owat_GLU245_OE1.dat
distance d_Owat_GLU245_OE2 :757@O :245@OE2 out after_min_dist_Owat_GLU245_OE2.dat

run
quit
EOF
```

Inspect:

```bash
cat after_min_dist_wat757_C2.dat
cat after_min_angle_wat757_C2_Asp38.dat
cat after_min_dist_Owat_GLU245_OE1.dat
cat after_min_dist_Owat_GLU245_OE2.dat
```

---

## 20. Restraining the covalent fructose and catalytic water during equilibration

The equilibration is also classical. The QM region is only activated during production.

In `step4.1_equilibration.mdin`, add:

```text
ntr=1,
```

inside the `&cntrl` block.

At the end of the file, after `&end`, add:

```text
Ligand 0CD and catalytic water 757 positional restraints
5.0
RES 601 601
RES 757 757
END
END
```

A complete final section should look like:

```text
 &wt
    type='END'
 /
DISANG=step4.1_equilibration.rest
LISTIN=POUT
LISTOUT=POUT
&end
Ligand 0CD and catalytic water 757 positional restraints
5.0
RES 601 601
RES 757 757
END
END
```

Run equilibration with the corrected topology:

```bash
pmemd.cuda -O \
  -i step4.1_equilibration.mdin \
  -p step3_input_noROH_CWT757.parm7 \
  -c step4.0_minimization.rst7 \
  -o step4.1_equilibration.mdout \
  -r step4.1_equilibration.rst7 \
  -inf step4.1_equilibration.mdinfo \
  -ref step3_input_noROH.rst7 \
  -x step4.1_equilibration.nc
```

Note that `-ref step3_input_noROH.rst7` is used as the positional restraint reference.

---

## 21. Checking catalytic water geometry during equilibration

After equilibration, analyze whether the catalytic water remained near the reactive region.

Create:

```bash
nano check_wat757_equilibration.in
```

Paste:

```cpptraj
parm step3_input_noROH_CWT757.parm7
trajin step4.1_equilibration.nc

distance d_wat757_C2 :757@O :601@C2 out eq_dist_wat757_C2.dat
distance d_Asp38_C2 :38@OD2 :601@C2 out eq_dist_Asp38_C2.dat
angle ang_attack :757@O :601@C2 :38@OD2 out eq_angle_wat757_C2_Asp38.dat

distance d_Owat_GLU245_OE1 :757@O :245@OE1 out eq_dist_Owat_GLU245_OE1.dat
distance d_Owat_GLU245_OE2 :757@O :245@OE2 out eq_dist_Owat_GLU245_OE2.dat

distance d_H1_GLU245_OE1 :757@H1 :245@OE1 out eq_dist_H1_GLU245_OE1.dat
distance d_H1_GLU245_OE2 :757@H1 :245@OE2 out eq_dist_H1_GLU245_OE2.dat
distance d_H2_GLU245_OE1 :757@H2 :245@OE1 out eq_dist_H2_GLU245_OE1.dat
distance d_H2_GLU245_OE2 :757@H2 :245@OE2 out eq_dist_H2_GLU245_OE2.dat

run
quit
```

Run:

```bash
cpptraj -i check_wat757_equilibration.in
```

Useful criteria:

```text
O_water — C2_0CD              ≈ 2.8–3.5 Å is ideal for direct attack
O_water — C2_0CD              > 5.0 Å indicates the water is still far
O_water — C2 — OD2_ASP38      > 130 degrees is a good attack orientation
H_water — OE1/OE2_GLU245      ≈ 1.6–2.3 Å suggests proton transfer compatibility
O_water — OE1/OE2_GLU245      ≈ 2.6–3.2 Å suggests hydrogen-bond compatibility
```

---

## 22. QM/MM production for the second reaction step

The QM/MM region is used only during production.

The reactive QM region should include at least:

```text
0CD/fructose residue 601
ASP38 side chain
ASP168 side chain
GLU/GLH245 side chain
CWT757 catalytic water
```

The key mechanistic process is:

```text
CWT757@O attacks 0CD@C2
ASP38@OD2 — 0CD@C2 bond breaks
a proton from CWT757 is transferred to GLU/GLH245
```

Because `ROH` was removed and the catalytic water was renamed to `CWT`, the production must use:

```text
step3_input_noROH_CWT757.parm7
```

not the original:

```text
step3_input.parm7
```

For the first production run after equilibration, use:

```bash
mpirun -np 30 sander.MPI -O \
  -i step5_production.mdin \
  -p step3_input_noROH_CWT757.parm7 \
  -c step4.1_equilibration.rst7 \
  -o step5_production.mdout \
  -r step5_production.rst7 \
  -inf step5_production.mdinfo \
  -x step5_production.nc
```

For restarting/continuing production, use:

```bash
mpirun -np 30 sander.MPI -O \
  -i step5_production.mdin \
  -p step3_input_noROH_CWT757.parm7 \
  -c step5_production.rst7 \
  -o step5_production_2.mdout \
  -r step5_production_2.rst7 \
  -inf step5_production_2.mdinfo \
  -x step5_production_2.nc
```

QM/MM production should be performed with `sander` or `sander.MPI`, because QM/MM calculations are not performed with `pmemd.cuda`.

---

## 23. Important warning about QM atom indices after removing ROH

After removing `ROH`, the atom indices change. Therefore, if the QM region in `step5_production.mdin` is defined using `iqmatoms`, the atom numbers must be updated.

For the corrected system:

```text
0CD/fructose = residue 601
0CD atoms    = 9029–9050
CWT757 water = residue 757
CWT757 atoms = 9206–9208
```

The QM region should include:

```text
ASP38 side chain
ASP168 side chain
GLU/GLH245 side chain
0CD residue 601
CWT757 residue 757
```

A representative QM atom list for the corrected topology is:

```text
ASP38 side chain:      575–580
ASP168 side chain:     2468–2473
GLU/GLH245 side chain: 3704–3712
0CD/fructose:          9029–9050
CWT757 water:          9206–9208
```

Example `&qmmm` block:

```text
 &qmmm
  iqmatoms=575,576,577,578,579,580,
           2468,2469,2470,2471,2472,2473,
           3704,3705,3706,3707,3708,3709,3710,3711,3712,
           9029,9030,9031,9032,9033,9034,9035,9036,9037,9038,9039,9040,
           9041,9042,9043,9044,9045,9046,9047,9048,9049,9050,
           9206,9207,9208,
  qmcharge=-2,
  qm_theory='DFTB3',
  qmcut=12.0,
  qmshake=0,
  dftb_slko_path='#cur_folder/parameters/3ob-3-1',
  adjust_q=1
 /
```

The total QM charge must be checked carefully according to the protonation states and selected atoms.

---

## 24. Mechanistic distances for the second reaction step

Create:

```bash
nano mechanism_distances_step2.in
```

Paste:

```cpptraj
parm step3_input_noROH_CWT757.parm7
trajin step5_production.nc

# Covalent bond to be broken
distance d_Asp38_OD2_C2_0CD :38@OD2 :601@C2 out step2_d_Asp38_OD2_C2_0CD.dat

# Water attack distance
distance d_WAT757_O_C2_0CD :757@O :601@C2 out step2_d_WAT757_O_C2_0CD.dat

# Attack angle: O_water — C2_fructose — OD2_Asp38
angle ang_WAT757_O_C2_Asp38_OD2 :757@O :601@C2 :38@OD2 out step2_angle_WAT757_O_C2_Asp38_OD2.dat

# Interaction between catalytic water and GLU/GLH245
distance d_WAT757_O_GLU245_OE1 :757@O :245@OE1 out step2_d_WAT757_O_GLU245_OE1.dat
distance d_WAT757_O_GLU245_OE2 :757@O :245@OE2 out step2_d_WAT757_O_GLU245_OE2.dat

# Proton transfer candidates from catalytic water to GLU/GLH245
distance d_WAT757_H1_GLU245_OE1 :757@H1 :245@OE1 out step2_d_WAT757_H1_GLU245_OE1.dat
distance d_WAT757_H1_GLU245_OE2 :757@H1 :245@OE2 out step2_d_WAT757_H1_GLU245_OE2.dat
distance d_WAT757_H2_GLU245_OE1 :757@H2 :245@OE1 out step2_d_WAT757_H2_GLU245_OE1.dat
distance d_WAT757_H2_GLU245_OE2 :757@H2 :245@OE2 out step2_d_WAT757_H2_GLU245_OE2.dat

# Stabilization by ASP168
distance d_Asp168_OD2_C2_0CD :168@OD2 :601@C2 out step2_d_Asp168_OD2_C2_0CD.dat
distance d_Asp168_OD1_C2_0CD :168@OD1 :601@C2 out step2_d_Asp168_OD1_C2_0CD.dat

run
quit
```

Run:

```bash
cpptraj -i mechanism_distances_step2.in
```

Expected interpretation:

```text
Hydrolysis of covalent intermediate:
Asp38 OD2 — 0CD C2 distance increases

Nucleophilic water attack:
CWT757 O — 0CD C2 distance decreases

Water attack geometry:
CWT757 O — 0CD C2 — Asp38 OD2 angle remains high, preferably >130 degrees

Proton transfer:
one of the CWT757 hydrogens approaches GLU/GLH245 OE1/OE2

Transition-state stabilization:
ASP168 remains close to the fructose reactive center
```

---

## 25. Summary of the second reaction step workflow

```bash
# 1. Remove artificial ROH introduced by CHARMM-GUI
./remove_roh_cpptraj.sh

# 2. Rename catalytic water 757 from WAT to CWT in the topology
python3 rename_wat757_to_cwt.py

# 3. Check catalytic water geometry before minimization
cpptraj -i check_wat757_before_min.in

# 4. Run classical minimization with restraints on 0CD and catalytic water
pmemd.cuda -O \
  -i step4.0_minimization.mdin \
  -p step3_input_noROH_CWT757.parm7 \
  -c step3_input_noROH.rst7 \
  -o step4.0_minimization.mdout \
  -r step4.0_minimization.rst7 \
  -inf step4.0_minimization.mdinfo \
  -ref step3_input_noROH.rst7

# 5. Check geometry after minimization
cpptraj -p step3_input_noROH_CWT757.parm7 << EOF
trajin step4.0_minimization.rst7
distance d_wat757_C2 :757@O :601@C2 out after_min_dist_wat757_C2.dat
angle ang_attack :757@O :601@C2 :38@OD2 out after_min_angle_wat757_C2_Asp38.dat
run
quit
EOF

# 6. Run classical equilibration with restraints on 0CD and catalytic water
pmemd.cuda -O \
  -i step4.1_equilibration.mdin \
  -p step3_input_noROH_CWT757.parm7 \
  -c step4.0_minimization.rst7 \
  -o step4.1_equilibration.mdout \
  -r step4.1_equilibration.rst7 \
  -inf step4.1_equilibration.mdinfo \
  -ref step3_input_noROH.rst7 \
  -x step4.1_equilibration.nc

# 7. Check geometry during equilibration
cpptraj -i check_wat757_equilibration.in

# 8. Run QM/MM production with sander.MPI
mpirun -np 30 sander.MPI -O \
  -i step5_production.mdin \
  -p step3_input_noROH_CWT757.parm7 \
  -c step4.1_equilibration.rst7 \
  -o step5_production.mdout \
  -r step5_production.rst7 \
  -inf step5_production.mdinfo \
  -x step5_production.nc

# 9. Restart/continue QM/MM production
mpirun -np 30 sander.MPI -O \
  -i step5_production.mdin \
  -p step3_input_noROH_CWT757.parm7 \
  -c step5_production.rst7 \
  -o step5_production_2.mdout \
  -r step5_production_2.rst7 \
  -inf step5_production_2.mdinfo \
  -x step5_production_2.nc

# 10. Analyze the mechanism of the second reaction step
cpptraj -i mechanism_distances_step2.in
```

---

## 26. Final corrected files for the second reaction step

The corrected system for the second reaction step consists of:

```text
step3_input_noROH_CWT757.parm7
step3_input_noROH.rst7
step4.0_minimization.mdin
step4.1_equilibration.mdin
step5_production.mdin
```

The original `step3_input.parm7` should not be used for the second step after ROH correction, because it still contains the artificial CHARMM-GUI covalent-ligand modification.

The correct production topology is:

```text
step3_input_noROH_CWT757.parm7
```

The correct production starting coordinates are:

```text
step4.1_equilibration.rst7
```

for the first production run, or:

```text
step5_production.rst7
```

for a restarted production run.
