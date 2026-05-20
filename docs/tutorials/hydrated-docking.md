# Molecular Docking with Explicit Waters

Using **AutoDock Vina 1.2.x**, **Meeko**, and **AutoGrid4**.

**Objective:** perform molecular docking while explicitly considering water molecules in the binding site.

**Recommended method:** AutoDock4 scoring with `--scoring ad4` and explicit waters attached to the ligand workflow.

## 1. Prerequisites and Installation

### 1.1 Install AutoDock Vina 1.2.x

```bash
# 1. Install compiler and dependencies
# Ubuntu/Debian
sudo apt update
sudo apt install build-essential libboost-all-dev swig

# 2. Download and compile Vina
git clone https://github.com/ccsb-scripps/AutoDock-Vina.git
cd AutoDock-Vina/build/linux/release
make

# 3. Put the executable in PATH
sudo cp vina /usr/local/bin/
# or add it to PATH:
# export PATH=$PATH:/path/to/AutoDock-Vina/build/linux/release
```

### 1.2 Required Scripts

Download the Meeko, molscrub, and AutoDock helper scripts:

```text
https://github.com/forlilab/Meeko/release/meeko/cli/mk_prepare_receptor.py
https://github.com/forlilab/Meeko/release/meeko/cli/mk_prepare_ligand.py
https://github.com/forlilab/molscrub/blob/develop/scripts/scrub.py
https://github.com/ccsb-scripps/AutoDock-Vina/develop/example/autodock_scripts/mapwater.py
https://github.com/ccsb-scripps/AutoDock-Vina/blob/develop/example/autodock_scripts/dry.py
```

### 1.3 Other Required Programs

- PDB2PQR for pH/protonation preparation
- PyMOL for inspection and conversion
- AutoGrid4, usually available with MGLTools or AutoDockTools

## 2. Receptor Preparation

```bash
# 1. Fix and protonate the receptor at the desired pH, here pH 5.0
pdb2pqr --ff=AMBER --ph-calc-method=propka --with-ph=5.0 --keep-chain --pdb-output receptor_pH5.0.pdb receptor_original.pdb receptor_pqr.pdb

# 2. Prepare the receptor for AutoDock with Meeko
python3 mk_prepare_receptor.py -i receptor_pH5.0.pdb -o receptor --box_center 12.0 9.0 -1.0 --box_size 25.5 31.5 25.5 -p -g
```

Generated files:

- `receptor.pdbqt`
- `receptor.gpf`
- `receptor.box.pdb`, useful for visual inspection

Open `receptor.box.pdb` in PyMOL to confirm that the docking box is well positioned.

## 3. Ligand Preparation

In PyMOL:

```pymol
load ligand.pdb
save ligand.sdf
```

Then run:

```bash
scrub.py ligand.sdf -o ligandH.sdf
```

Prepare the hydrated ligand:

```bash
python3 mk_prepare_ligand.py -i ligandH.sdf -o ligand_hydrated.pdbqt -w
```

## 4. Generate Affinity Maps with AutoGrid4

```bash
# Check ligand atom types
grep "^ATOM\|^HETATM" ligand_hydrated.pdbqt | cut -c 78-79 | sort | uniq

# Edit receptor.gpf and keep only the atom types that exist in the ligand.
# Common example for organic molecules:
# ligand_types C OA HD H

# Run AutoGrid4
autogrid4 -p receptor.gpf -l receptor.glg
```

## 5. Create the Special Water Map

```bash
python3 mapwater.py -r receptor.pdbqt -s receptor.W.map
```

## 6. Run Docking

For hydrated docking, use AutoDock4 scoring:

```bash
vina --ligand ligand_hydrated.pdbqt --maps receptor --scoring ad4 --exhaustiveness 32 --out ligand_out_ad4.pdbqt
```

## 7. Post-process and Classify Waters

```bash
python3 dry.py -r receptor.pdbqt -m receptor.W.map -i ligand_out_ad4.pdbqt
```

The output file has the suffix `_DRY_SCORED.pdbqt` and classifies waters as:

- `STRONG`: well-retained water
- `WEAK`: marginal water
- `DISPLC`: likely displaced water

## 8. Final Visualization

```bash
pymol receptor_pH5.0.pdb receptor.box.pdb ligand_out_ad4_DRY_SCORED.pdbqt
```

## Important Files

| Step | Main generated file | Purpose |
| --- | --- | --- |
| Protonation | `receptor_pH5.0.pdb`, `receptor_pqr.pdb` | Corrected and protonated receptor |
| Receptor preparation | `receptor.pdbqt`, `receptor.gpf` | AutoDock receptor files |
| Ligand preparation | `ligand_hydrated.pdbqt` | Ligand with explicit waters |
| Maps | `receptor.*.map`, `receptor.W.map` | Affinity maps and water map |
| Docking | `ligand_out_ad4.pdbqt` | Docking results |
| Post-processing | `ligand_out_ad4_DRY_SCORED.pdbqt` | Final result with classified waters |

## Final Notes

- Always use `--scoring ad4` for this hydrated docking workflow.
- Increase `--exhaustiveness` to 32-64 for fragments or flexible ligands.
- This method is not recommended for virtual screening of thousands of compounds.
