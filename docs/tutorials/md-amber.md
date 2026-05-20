# Classical Molecular Dynamics with AMBER

This tutorial describes a classical protein-ligand molecular dynamics workflow with **AMBER**.

It assumes that the protein-ligand complex was prepared with **CHARMM-GUI Solution Builder** and that AMBER input files are available.

## 1. Generate `dihe.restraint` for Glycans

This section explains how to identify a glycan ligand prepared by CHARMM-GUI and manually generate a `dihe.restraint` file for AMBER molecular dynamics.

This is especially useful when CHARMM-GUI does not automatically create dihedral restraints for carbohydrates or polymers, which is common for fructans and inulin. If CHARMM-GUI already created the file, skip this section and continue to minimization.

The example uses a beta-(2->1) fructose oligomer containing residues `571-576`.

### Context

- system prepared with CHARMM-GUI
- MD engine: AMBER
- ligand: inulin, a beta-fructan
- ligand residues: `1CU` and `0CU`
- main files:
  - `amber/step3_input.parm7`
  - `amber/step3_input.rst7`
  - `residues.dat`
  - `calc_inulin_dihedrals.in`
  - `write_dihe_restraint.py`

### Identify Ligand Residues

Inside the system directory, generate the residue list:

```bash
cpptraj step3_input.parm7 << EOF > residues.dat
resinfo :*
EOF
```

In this example, the ligand was identified as:

```text
  571 1CU    8794   8814     21   571     2
  572 1CU    8815   8835     21   572     2
  573 1CU    8836   8856     21   573     2
  574 1CU    8857   8877     21   574     2
  575 1CU    8878   8898     21   575     2
  576 0CU    8899   8920     22   576     2
```

This indicates a chain of six fructose units, with `0CU` as the terminal residue.

### Inspect Atom Names and Atom Indices

Open `cpptraj`:

```bash
cpptraj step3_input.parm7
```

List atoms for one ligand residue:

```cpptraj
atominfo :571
```

Relevant atoms for glycosidic dihedrals:

- `O5`
- `C2`
- `O1`
- `C1`

Repeat the inspection for the next residue:

```cpptraj
atominfo :572
```

### Correct Dihedrals for beta-(2->1) Fructan

For each glycosidic bond between residue `i` and residue `i+1`, restrain only the glycosidic dihedrals.

Phi:

```text
O5(i) - C2(i) - O1(i+1) - C1(i+1)
```

Psi:

```text
C2(i) - O1(i+1) - C1(i+1) - C2(i+1)
```

Do not restrain internal ring dihedrals.

For beta-(2->6) fructan, such as levan, the anomeric carbon `C2` of the terminal residue `0CU` binds to `O6` of the previous residue. The glycosidic axis is:

```text
O5(i) - C6(i) - O6(i) - C2(i+1)
```

Phi:

```text
O5(i) - C6(i) - O6(i) - C2(i+1)
```

Psi:

```text
C6(i) - O6(i) - C2(i+1) - O5(i+1)
```

## 2. Calculate Dihedrals with `cpptraj`

Create `calc_inulin_dihedrals.in`:

```cpptraj
parm step3_input.parm7
trajin step3_input.rst7 1 1

# Bond 571-572
dihedral phi_571_572 :571@O5 :571@C2 :572@O1 :572@C1 out phi_571_572.dat
dihedral psi_571_572 :571@C2 :572@O1 :572@C1 :572@O5 out psi_571_572.dat

# Bond 572-573
dihedral phi_572_573 :572@O5 :572@C2 :573@O1 :573@C1 out phi_572_573.dat
dihedral psi_572_573 :572@C2 :573@O1 :573@C1 :573@O5 out psi_572_573.dat

# Bond 573-574
dihedral phi_573_574 :573@O5 :573@C2 :574@O1 :574@C1 out phi_573_574.dat
dihedral psi_573_574 :573@C2 :574@O1 :574@C1 :574@O5 out psi_573_574.dat

# Bond 574-575
dihedral phi_574_575 :574@O5 :574@C2 :575@O1 :575@C1 out phi_574_575.dat
dihedral psi_574_575 :574@C2 :575@O1 :575@C1 :575@O5 out psi_574_575.dat

# Bond 575-576
dihedral phi_575_576 :575@O5 :575@C2 :576@O1 :576@C1 out phi_575_576.dat
dihedral psi_575_576 :575@C2 :576@O1 :576@C1 :576@O5 out psi_575_576.dat

run
```

Run:

```bash
cpptraj -i calc_inulin_dihedrals.in
```

Check the generated values:

```bash
head phi_571_572.dat
head psi_571_572.dat
```

## 3. Create `dihe.restraint`

Create `write_dihe_restraint.py`:

```python
import glob

RK = 20.0
DELTA = 10.0
OUTFILE = "dihe.restraint"

DIHEDRALS = {
    # PHI = O5(i) - C2(i) - O1(i+1) - C1(i+1)
    "phi_571_572": [8795, 8794, 8835, 8832],
    "phi_572_573": [8816, 8815, 8856, 8853],
    "phi_573_574": [8837, 8836, 8877, 8874],
    "phi_574_575": [8858, 8857, 8898, 8895],
    "phi_575_576": [8879, 8878, 8920, 8916],

    # PSI = C2(i) - O1(i+1) - C1(i+1) - O5(i+1)
    "psi_571_572": [8794, 8835, 8832, 8816],
    "psi_572_573": [8815, 8856, 8853, 8837],
    "psi_573_574": [8836, 8877, 8874, 8858],
    "psi_574_575": [8857, 8898, 8895, 8879],
    "psi_575_576": [8878, 8920, 8916, 8900],
}

def read_dihedral_value(filename):
    with open(filename) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                return float(line.split()[1])
    raise RuntimeError(f"Value not found in {filename}")

def build_restraint_window(angle, delta):
    r1 = -180.0
    r4 = 180.0
    r2 = max(angle - delta, -180.0)
    r3 = min(angle + delta, 180.0)

    if r2 >= r3:
        raise ValueError(
            f"Invalid window for dihedral {angle:.2f} "
            f"(r2={r2:.2f}, r3={r3:.2f})"
        )

    return r1, r2, r3, r4

with open(OUTFILE, "w") as out:
    dat_files = sorted(glob.glob("phi_*.dat") + glob.glob("psi_*.dat"))
    for dat in dat_files:
        key = dat.replace(".dat", "")
        if key not in DIHEDRALS:
            continue

        angle = read_dihedral_value(dat)
        r1, r2, r3, r4 = build_restraint_window(angle, DELTA)
        iat = DIHEDRALS[key]

        out.write("&rst\n")
        out.write(f" iat={iat[0]},{iat[1]},{iat[2]},{iat[3]},\n")
        out.write(f" r1={r1:.1f}, r2={r2:.1f}, r3={r3:.1f}, r4={r4:.1f},\n")
        out.write(f" rk2={RK:.1f}, rk3={RK:.1f},\n")
        out.write("/\n\n")

print("dihe.restraint generated.")
```

Run:

```bash
python write_dihe_restraint.py
```

Expected output:

```text
dihe.restraint
```

## 4. Use `dihe.restraint` in AMBER Inputs

In `step4.0_minimization.mdin` and `step4.1_equilibration.mdin`, enable NMR-style restraints:

```ini
&cntrl
  nmropt=1,
/
&wt type='END' /
DISANG=dihe.restraint
```

Recommended practice:

- use `dihe.restraint` only until the end of equilibration
- remove it completely for production
- use `rk2 = rk3 = 10-20` for typical glycan stabilization
- never restrain internal ring dihedrals

This helps preserve the initial glycan conformation and avoid nonphysical collapse at the beginning of the simulation.

## 5. GPU Selection for AMBER

To select GPU device 1, use `CUDA_VISIBLE_DEVICES`:

```bash
export CUDA_VISIBLE_DEVICES=1
```

or inline:

```bash
CUDA_VISIBLE_DEVICES=1 pmemd.cuda ...
```

For GPU-accelerated molecular dynamics, replace `pmemd` with `pmemd.cuda`:

```bash
pmemd.cuda -O -i mdin -o mdout -p prmtop -c inpcrd -r restrt -x mdcrd
```

If several GPUs are available, list candidate GPUs with commas:

```bash
export CUDA_VISIBLE_DEVICES=1,3
pmemd.cuda -O -i mdin -o mdout -p prmtop -c inpcrd -r restrt -x mdcrd
```

## 6. Minimization, Equilibration, and Production

### 6.1 Minimization

```bash
sed -e "s/FC/1.0/g" dihe.restraint > step4.0_minimization.rest
```

```bash
pmemd.cuda -O -i step4.0_minimization.mdin -p step3_input.parm7 -c step3_input.rst7 -o step4.0_minimization.mdout -r step4.0_minimization.rst7 -inf step4.0_minimization.mdinfo -ref step3_input.rst7
```

### 6.2 Equilibration

```bash
sed -e "s/FC/1.0/g" dihe.restraint > step4.1_equilibration.rest
```

```bash
pmemd.cuda -O -i step4.1_equilibration.mdin -p step3_input.parm7 -c step4.0_minimization.rst7 -o step4.1_equilibration.mdout -r step4.1_equilibration.rst7 -inf step4.1_equilibration.mdinfo -ref step3_input.rst7 -x step4.1_equilibration.nc
```

### 6.3 Production

```bash
pmemd.cuda -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5_production.mdout -r step5_production.rst7 -inf step5_production.mdinfo -x step5_production.nc
```

Alternative with `sander.MPI`:

```bash
mpirun -np 4 sander.MPI -O -i step5_production.mdin -p step3_input.parm7 -c step4.1_equilibration.rst7 -o step5_production.mdout -r step5_production.rst7 -inf step5_production.mdinfo -x step5_production.nc
```

Here, `-np 4` is the number of CPU cores used.

## 7. Production Replicas

Add the `ig` flag in `step5_production.mdin` to generate different trajectories:

```text
ig = -1
```

Run three replicas:

```bash
for i in 1 2 3
do
    echo "Running replica $i"

    pmemd.cuda -O \
    -i step5_production.mdin \
    -p step3_input.parm7 \
    -c step4.1_equilibration.rst7 \
    -o step5_production_${i}.mdout \
    -r step5_production_${i}.rst7 \
    -inf step5_production_${i}.mdinfo \
    -x step5_production_${i}.nc

done
```

## 8. AMBER Post-processing with `cpptraj`

Main files:

- topology: `step3_input.parm7`
- production trajectory: `step5_production.nc`
- reference coordinates: `step5_production.rst7`

The ligand corresponds to residues `1CU` and `0CU`.

Replace `XXX` with the last protein residue, excluding solvent and ligand.

### 8.1 Center the Trajectory

```bash
cpptraj -p step3_input.parm7 -y step5_production.nc << EOF
autoimage anchor :1-XXX
center :1-XXX mass origin
image origin center
rms first :1-XXX@CA
trajout step5_centered.nc

EOF
```

### 8.2 Ligand RMSD

```bash
cpptraj -p step3_input.parm7 -y step5_centered.nc << EOF
rms Protein first :1-XXX@CA
rms Ligand first :1CU,0CU out rmsd_ligand.dat
EOF
```

### 8.3 Residue RMSF

```bash
cpptraj -p step3_input.parm7 -y step5_centered.nc << EOF
rms first :1-XXX@CA
atomicfluct out rmsf_ca.dat :1-XXX@CA byres
EOF
```

### 8.4 Radius of Gyration

```bash
cpptraj -p step3_input.parm7 -y step5_centered.nc << EOF
radgyr :1-XXX out rg_protein.dat
EOF
```

### 8.5 Ligand-Protein Hydrogen Bonds

```bash
cpptraj -p step3_input.parm7 -y step5_centered.nc << EOF
hbond HB out hbond_lig_prot.dat \
donormask :1CU,0CU \
acceptormask :1-XXX
EOF
```

To include solvent-ligand and solvent-protein hydrogen bonds:

```bash
cpptraj -p step3_input.parm7 -y step5_centered.nc << EOF
hbond HB out hbond_lig_prot.dat \
solventdonor :WAT \
donormask :1CU,0CU \
acceptormask :1-XXX
EOF
```

### 8.6 Clustering

Cluster ligand conformations by RMSD:

```bash
cpptraj -p step3_input.parm7 -y step5_centered.nc << EOF
rms first :1CU,0CU
cluster hieragglo clusters 3 linkage average \
  summary cluster_summary.dat \
  repout cluster_rep repfmt pdb
EOF
```

Expected files:

- `cluster_summary.dat`
- `cluster_rep.c0.pdb`

The repository also includes scripts for repeated analyses:

```bash
bash run_cpptraj_replicatas.sh
python plot_md_replicates.py
bash cluster_3_replicates.sh
```

## 9. MM/PBSA and MM/GBSA with AMBER

MM/PBSA and MM/GBSA estimate receptor-ligand binding free energy from MD snapshots. This requires separate dry topologies for:

1. complex
2. receptor
3. ligand

### 9.1 Prepare Dry Topologies

Identify the receptor and ligand residue masks. For example, if the protein is residues `1-630` and the ligand is the remaining non-solvent residue:

```bash
ante-MMPBSA.py -p step3_input.parm7 \
               -c complex.prmtop \
               -r receptor.prmtop \
               -l ligand.prmtop \
               -m ':1-630' \
               -s ':WAT,Na+,Cl-' \
               --radii=mbondi2
```

The `-s` mask strips water and ions. Adjust residue masks for the actual system.

### 9.2 Create `mmpbsa.in`

```text
&general
   interval=1,
   verbose=1,
   keep_files=0,
   strip_mask=":WAT,Na+,Cl-"
/
&gb
   igb=5,
   saltcon=0.100
/
&pb
   istrng=0.100
/
&decomp
   idecomp=1,
   dec_verbose=1
/
```

### 9.3 Run Serial MMPBSA

```bash
MMPBSA.py -O -i mmpbsa.in \
          -o FINAL_RESULTS_MMPBSA.dat \
          -do FINAL_DECOMP_MMPBSA.dat \
          -sp step3_input.parm7 \
          -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop \
          -y step5_production.nc
```

### 9.4 Run Parallel MMPBSA

```bash
mpirun -np 16 MMPBSA.py.MPI -O -i mmpbsa.in \
          -o FINAL_RESULTS_MMPBSA.dat \
          -do FINAL_DECOMP_MMPBSA.dat \
          -sp step3_input.parm7 \
          -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop \
          -y step5_production.nc > mmpbsa_progress.log 2>&1
```

Do not use more MPI processes than trajectory frames. More processes also increase memory usage because each process loads part of the calculation context.

### 9.5 Interpret Results

`FINAL_RESULTS_MMPBSA.dat` reports average energy components for complex, receptor, ligand, and the difference:

```text
Complex - Receptor - Ligand
```

The resulting `TOTAL` value is the estimated binding free energy. Negative values indicate favorable binding. The calculation above does not include conformational entropy, so it is best used for relative comparisons unless an entropy term is added.

`FINAL_DECOMP_MMPBSA.dat` reports per-residue energetic contributions. Negative values favor binding, while positive values disfavor binding.
