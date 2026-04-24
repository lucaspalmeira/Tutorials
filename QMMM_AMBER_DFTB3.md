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

# 1. Extract representative cluster frame

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

# 2. Measure mechanism distances

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

# 3. Measure distances only in representative frame

File: rep_frame_measure.in

```cpptraj
parm step3_input.parm7
trajin step5_centered.nc 10 10

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

# 4. Python plotting script

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

# 5. Expected outputs

cluster2_rep_frame10.pdb
mechanism_distances.agr
mechanism_distances.png

---

# 6. Interpretation guidelines

Covalent intermediate formation:

Asp38 OD2 – C2 distance ≈ 1.4–1.7 Å

Bond breaking:

C2 – O1 distance increases

Proton transfer:

GLH245 HE2 approaches O1

Transition-state stabilization:

Asp168 interacts with C2 or O1

---

# 7. Full workflow

```bash
cpptraj -i extract_cluster_rep.in
cpptraj -i mechanism_distancesd.in
cpptraj -i rep_frame_measure.in
python plot_mechanism_distances.py
```
