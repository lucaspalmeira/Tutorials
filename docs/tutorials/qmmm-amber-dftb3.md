# DFTB3-Based QM/MM MD in AMBER

This tutorial describes a DFTB3-based QM/MM molecular dynamics workflow in AMBER for sucrose hydrolysis by invertase.

## Introduction

DFTB3 is an extension of SCC-DFTB, also known as DFTB2, based on a third-order expansion of the total DFT energy. This improves the description of charge variation, electrostatic interactions, and hydrogen bonds in biomolecular systems. Unlike SCC-DFTB, where the Hubbard parameter is fixed for each element, DFTB3 allows chemical hardness to depend on atomic charge state, making it more robust for reactions involving electronic redistribution.

In this tutorial, the system is invertase complexed with sucrose, prepared with CHARMM-GUI at pH 4.0. The AMBER simulation uses the FF19SB force field for the classical region, a temperature of 323.15 K, and a total simulation time of 1 ns. This setup provides a biologically consistent starting point with protonation states adjusted to an acidic environment compatible with invertase activity.

The reaction of interest is sucrose hydrolysis, characterized by cleavage of the beta-glycosidic beta(2->1) bond between glucose and fructose. Because this process involves bond breaking and bond formation, a purely classical treatment is not sufficient. DFTB3 within a QM/MM approach treats the reactive region quantum mechanically while the protein, solvent, and ions remain under molecular mechanics.

The workflow analyzes:

- nucleophilic attack of Asp38 on fructose C2
- proton donation from GLH245
- transition-state stabilization by Asp168

The reactive fructose residue is `0CU`, and the reactive atom is `0CU@C2`.

Target bond formation:

```text
Asp38@OD2 -> 0CU@C2
```

## 1. Run Production

```bash
mpirun -np 20 sander.MPI -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5_production.mdout -r step5_production.rst7 -inf step5_production.mdinfo -x step5_production.nc
```

Restart production:

```bash
mpirun -np 20 sander.MPI -O -i step5_production.mdin -p step3_input.parm7 -c step5_production.rst7 -o step5_production_2.mdout -r step5_production_2.rst7 -inf step5_production_2.mdinfo -x step5_production_2.nc
```

## 2. Extract Representative Cluster Frames

File: `extract_cluster_rep.in`

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

## 3. Measure Mechanism Distances

File: `mechanism_distances.in`

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

## 4. Measure Distances in the Representative Frame

File: `rep_frame_measure.in`

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

## 5. Plot Mechanism Distances

File: `plot_mechanism_distances.py`

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
    x, y = [], []
    for line in open(filepath):
        if line.startswith("#") or not line.strip():
            continue
        frame, value = line.split()[:2]
        x.append(float(frame))
        y.append(float(value))
    return np.array(x), np.array(y)

plt.figure(figsize=(10, 6))

for filename, label in FILES_AND_LABELS.items():
    if Path(filename).exists():
        x, y = read_cpptraj_dat(filename)
        plt.plot(x, y, label=label)

plt.axvline(REP_FRAME, linestyle="--")

plt.xlabel("Frame")
plt.ylabel("Distance (Angstrom)")
plt.legend()
plt.tight_layout()

plt.savefig(OUTPUT_FIG, dpi=600)
plt.show()
```

Run:

```bash
python plot_mechanism_distances.py
```

## 6. Expected Outputs

- `cluster2_rep_frame10.pdb`
- `mechanism_distances.agr`
- `mechanism_distances.png`

## 7. Interpretation Guidelines

| Event | Expected behavior |
| --- | --- |
| Covalent intermediate formation | Asp38 OD2 - C2 distance around 1.4-1.7 angstrom |
| Bond breaking | C2 - O1 distance increases |
| Proton transfer | GLH245 HE2 approaches O1 |
| Transition-state stabilization | Asp168 interacts with C2 or O1 |

## 8. Full Workflow

```bash
cpptraj -i extract_cluster_rep.in
cpptraj -i mechanism_distances.in
cpptraj -i rep_frame_measure.in
python plot_mechanism_distances.py
```
