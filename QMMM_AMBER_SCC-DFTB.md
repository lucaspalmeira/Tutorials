# QM/MM Reaction Analysis – CPPTRAJ + Python Workflow

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
